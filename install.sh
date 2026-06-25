#!/bin/bash
# install.sh — Hermes Medical Tools installer
# Copies tools into Hermes Agent and configures the medical toolset.
# Run from the repo root.

set -e

HERMES_HOME="${HERMES_HOME:-$HOME/.hermes}"
HERMES_TOOLS="${HERMES_HOME}/hermes-agent/tools"
TOOLSETS_FILE="${HERMES_HOME}/hermes-agent/toolsets.py"

echo "Hermes Medical Tools — Installer"
echo "================================="
echo ""

# 1. Verify Hermes Agent is installed
if [ ! -d "$HERMES_TOOLS" ]; then
    echo "ERROR: Hermes Agent tools directory not found at $HERMES_TOOLS"
    echo "       Install Hermes Agent first: https://github.com/NousResearch/hermes-agent"
    exit 1
fi
echo "✓ Hermes Agent found at $HERMES_TOOLS"

# 2. Check backend availability
echo ""
echo "Checking backends..."

HAS_PSPP=0
HAS_R=0
HAS_SCIPY=0

if command -v pspp &> /dev/null; then
    HAS_PSPP=1
    echo "  ✓ pspp: $(pspp --version 2>&1 | head -1)"
else
    echo "  ⚠ pspp: not found (install: sudo apt install pspp)"
fi

if command -v Rscript &> /dev/null; then
    HAS_R=1
    echo "  ✓ Rscript: $(Rscript --version 2>&1)"
else
    echo "  ⚠ Rscript: not found (install: sudo apt install r-base)"
fi

if python3 -c "import scipy" 2>/dev/null; then
    HAS_SCIPY=1
    echo "  ✓ scipy: available"
else
    echo "  ⚠ scipy: not found (install: pip install scipy numpy)"
fi

# 3. Copy tools
echo ""
echo "Installing tools..."

TOOLS=("medical_tools.py" "pspp_tool.py" "jmv_tool.py")
for tool in "${TOOLS[@]}"; do
    if [ -f "tools/$tool" ]; then
        cp "tools/$tool" "$HERMES_TOOLS/$tool"
        echo "  ✓ $tool → $HERMES_TOOLS/$tool"
    else
        echo "  ✗ $tool: not found in tools/"
        exit 1
    fi
done

# 4. Add medical toolset to toolsets.py
echo ""
echo "Configuring medical toolset..."

if grep -q '"medical"' "$TOOLSETS_FILE"; then
    echo "  ⚠ medical toolset already exists in toolsets.py — skipping"
else
    # Find the spotify toolset entry (a good insertion point before "Scenario-specific")
    if grep -q '"spotify"' "$TOOLSETS_FILE"; then
        python3 -c "
import re
path = '$TOOLSETS_FILE'
with open(path) as f:
    content = f.read()

med_block = '''
    \"medical\": {
        \"description\": \"Medical research — jmv/R, PSPP, PubMed, clinical trials, EBM\",
        \"tools\": [\"jmv\", \"pspp\", \"med_pubmed\", \"med_trial\", \"med_stats\", \"med_power\", \"med_evidence\"],
        \"includes\": [],
    },

    # Scenario-specific toolsets
'''

# Insert after spotify block
content = re.sub(
    r'(\\s*# Scenario-specific toolsets)',
    med_block,
    content,
    count=1
)

with open(path, 'w') as f:
    f.write(content)
print('added')
"
        echo "  ✓ medical toolset added to toolsets.py"
    else
        echo "  ⚠ Could not auto-insert toolset. Add manually:"
        echo ""
        echo '  "medical": {'
        echo '      "description": "Medical research — jmv/R, PSPP, PubMed, clinical trials, EBM",'
        echo '      "tools": ["jmv", "pspp", "med_pubmed", "med_trial", "med_stats", "med_power", "med_evidence"],'
        echo '      "includes": [],'
        echo '  },'
        echo ""
        echo "  to the TOOLSETS dict in $TOOLSETS_FILE"
    fi
fi

# 5. Summary
echo ""
echo "================================="
echo "Installation complete."
echo ""
echo "Available tools based on your system:"
[ $HAS_PSPP -eq 1 ] && echo "  ✓ pspp    — 135 SPSS commands (GNU PSPP)"
[ $HAS_PSPP -eq 0 ] && echo "  ✗ pspp    — install: sudo apt install pspp"
[ $HAS_R -eq 1 ] && echo "  ✓ jmv     — 9 statistical analyses (R base stats)"
[ $HAS_R -eq 0 ] && echo "  ✗ jmv     — install: sudo apt install r-base"
[ $HAS_SCIPY -eq 1 ] && echo "  ✓ med_stats, med_power — quick stats (scipy)"
[ $HAS_SCIPY -eq 0 ] && echo "  ✗ med_stats, med_power — install: pip install scipy numpy"
echo "  ✓ med_pubmed     — always available (stdlib)"
echo "  ✓ med_trial      — always available (stdlib)"
echo "  ✓ med_evidence   — always available (stdlib)"
echo ""
echo "Next step: hermes tools enable medical"
echo "           Then start a new Hermes session."
