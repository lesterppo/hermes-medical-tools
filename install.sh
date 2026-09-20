#!/usr/bin/env bash
# install.sh — Hermes Medical Tools installer (plugin-based)
#
# Installs the tools as a Hermes PLUGIN under $HERMES_HOME/plugins/ so that
# `hermes update` cannot wipe them. The old behaviour (copy into
# hermes-agent/tools/) is still available with --legacy, but it is documented
# as lossy: updates reset that directory and any toolsets.py entry with it.
#
# Usage:
#   ./install.sh              # install/refresh the plugin
#   ./install.sh --legacy     # copy into hermes-agent/tools/ (update-lossy)
#   ./install.sh --uninstall  # remove the plugin
#   ./install.sh --check      # report backend availability only

set -euo pipefail

HERMES_HOME="${HERMES_HOME:-$HOME/.hermes}"
PLUGIN_NAME="hermes_medical_tools"
PLUGIN_DIR="$HERMES_HOME/plugins/$PLUGIN_NAME"
legacy=0
uninstall=0
check_only=0

for arg in "$@"; do
    case "$arg" in
        --legacy) legacy=1 ;;
        --uninstall) uninstall=1 ;;
        --check) check_only=1 ;;
        -h|--help)
            sed -n '2,16p' "$0"
            exit 0
            ;;
        *)
            echo "Unknown option: $arg" >&2
            exit 2
            ;;
    esac
done

echo "Hermes Medical Tools — Installer"
echo "================================"

# ---------------------------------------------------------------------------
# 0. Uninstall
# ---------------------------------------------------------------------------
if [ "$uninstall" -eq 1 ]; then
    if [ -d "$PLUGIN_DIR" ]; then
        rm -rf "$PLUGIN_DIR"
        echo "✓ removed $PLUGIN_DIR"
        echo "  (restart Hermes / the gateway to unload the toolset)"
    else
        echo "nothing to remove at $PLUGIN_DIR"
    fi
    exit 0
fi

# ---------------------------------------------------------------------------
# 1. Backend detection
# ---------------------------------------------------------------------------
HAS_PSPP=0; HAS_R=0; HAS_SCIPY=0
command -v pspp >/dev/null 2>&1 && HAS_PSPP=1
command -v Rscript >/dev/null 2>&1 && HAS_R=1
python3 -c "import scipy" >/dev/null 2>&1 && HAS_SCIPY=1

echo
echo "Backends detected:"
if [ "$HAS_SCIPY" -eq 1 ]; then
    echo "  ✓ scipy        — med_stats, med_power"
else
    echo "  ✗ scipy        — install: pip install scipy numpy"
fi
if [ "$HAS_PSPP" -eq 1 ]; then
    echo "  ✓ pspp         — SPSS-syntax analyses (GNU PSPP)"
else
    echo "  ✗ pspp         — install: sudo apt install pspp"
fi
if [ "$HAS_R" -eq 1 ]; then
    echo "  ✓ Rscript      — jamovi-compatible analyses"
else
    echo "  ✗ Rscript      — install: sudo apt install r-base"
fi
echo "  ✓ stdlib tools — med_pubmed, med_trial, med_evidence (no backend needed)"

if [ "$check_only" -eq 1 ]; then
    exit 0
fi

# ---------------------------------------------------------------------------
# 2. Locate the Hermes install
# ---------------------------------------------------------------------------
if [ ! -d "$HERMES_HOME" ]; then
    echo
    echo "ERROR: Hermes home not found at $HERMES_HOME" >&2
    echo "       Install Hermes Agent first: https://github.com/NousResearch/hermes-agent" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 3. Legacy copy path (kept for compatibility, documented as lossy)
# ---------------------------------------------------------------------------
if [ "$legacy" -eq 1 ]; then
    HERMES_TOOLS="$HERMES_HOME/hermes-agent/tools"
    TOOLSETS_FILE="$HERMES_HOME/hermes-agent/toolsets.py"
    if [ ! -d "$HERMES_TOOLS" ]; then
        echo "ERROR: $HERMES_TOOLS not found" >&2
        exit 1
    fi
    echo
    echo "WARNING: --legacy copies into the Hermes git tree. The next"
    echo "         'hermes update' deletes these files and the toolsets.py"
    echo "         entry, leaving 'Unknown toolsets: medical'. Prefer the"
    echo "         plugin install (no flag)."
    for tool in medical_tools.py pspp_tool.py jmv_tool.py; do
        cp "tools/$tool" "$HERMES_TOOLS/$tool"
        echo "  ✓ $tool → $HERMES_TOOLS/$tool"
    done
    if ! grep -q '"medical"' "$TOOLSETS_FILE" 2>/dev/null; then
        echo
        echo "Add the toolset manually to $TOOLSETS_FILE:"
        echo '  "medical": {"description": "Medical research", '
        echo '              "tools": ["jmv","pspp","med_pubmed","med_trial","med_stats","med_power","med_evidence"],'
        echo '              "includes": []},'
    fi
    echo
    echo "Then: hermes tools enable medical && start a new session."
    exit 0
fi

# ---------------------------------------------------------------------------
# 4. Plugin install (recommended)
# ---------------------------------------------------------------------------
echo
echo "Installing plugin → $PLUGIN_DIR"

# Collision check: another plugin registering the same tool names would fight
# over the registry slot (last registration wins).
collisions=""
for other in "$HERMES_HOME"/plugins/*/; do
    [ -d "$other" ] || continue
    base="$(basename "$other")"
    [ "$base" = "$PLUGIN_NAME" ] && continue
    if grep -rqs "med_pubmed\|med_trial\|med_stats\|med_power\|med_evidence\|\"pspp\"\|'pspp'\|\"jmv\"\|'jmv'" "$other" 2>/dev/null; then
        collisions="$collisions $base"
    fi
done
if [ -n "$collisions" ]; then
    echo
    echo "  ⚠ These plugins already reference the same tool names:$collisions"
    echo "    Duplicate registrations overwrite each other. Disable one set:"
    for c in $collisions; do
        echo "      hermes plugins disable $(basename "$c")"
    done
fi

mkdir -p "$PLUGIN_DIR/tools"
cp plugin/plugin.yaml "$PLUGIN_DIR/plugin.yaml"
cp plugin/__init__.py "$PLUGIN_DIR/__init__.py"
cp plugin/tools/__init__.py "$PLUGIN_DIR/tools/__init__.py"
for tool in medical_tools.py pspp_tool.py jmv_tool.py; do
    cp "tools/$tool" "$PLUGIN_DIR/tools/$tool"
    echo "  ✓ tools/$tool"
done

echo
echo "Enabled toolset: medical"
echo
echo "Next steps:"
echo "  1. hermes plugins enable $PLUGIN_NAME      # if not auto-enabled"
echo "  2. restart Hermes (or the gateway) so the plugin loads"
echo "  3. the 'medical' toolset is auto-enabled for the platform"
echo
echo "Verify with:  hermes tools list | grep -E 'med_|pspp|jmv'"
