# AGENTS.md — Hermes Medical Tools

For AI coding assistants and Hermes agents working with this repo.

## What This Is

7 native, token-efficient medical research tools for Hermes Agent. They register via `tools.registry` and are auto-discovered by Hermes's tool loader. Each tool is service-gated — it only appears in the agent's toolset when its backend is installed.

## File Structure

```
hermes-medical-tools/
├── tools/
│   ├── medical_tools.py    # med_pubmed, med_trial, med_stats, med_power, med_evidence
│   ├── pspp_tool.py         # PSPP subprocess wrapper
│   └── jmv_tool.py          # R subprocess wrapper (jamovi-compatible)
├── README.md                # User-facing install + usage
├── SKILL.md                 # Hermes skill file
├── AGENTS.md                # This file
└── install.sh               # Automated install script
```

## How Tools Register

Each tool file calls `registry.register()` at module level:

```python
from tools.registry import registry

registry.register(
    name="tool_name",
    toolset="medical",
    schema={...},
    handler=lambda args, **kw: tool_fn(...),
    check_fn=availability_check,
    emoji="...",
)
```

Hermes auto-discovers tools in `tools/*.py` by scanning for `registry.register()` calls. No manual import list needed.

## Adding to toolsets.py

The user must add the `medical` toolset to `~/.hermes/hermes-agent/toolsets.py`:

```python
"medical": {
    "description": "Medical research — jmv/R, PSPP, PubMed, clinical trials, EBM",
    "tools": ["jmv", "pspp", "med_pubmed", "med_trial", "med_stats", "med_power", "med_evidence"],
    "includes": [],
},
```

## Design Principles

- **Token efficiency**: Short param names, compact JSON output, terse descriptions
- **Service gating**: `check_fn` verifies backend availability — zero schema cost when unavailable
- **Single purpose**: Each tool does one thing well
- **Direct computation**: Tools do the work, return compact results — no MCP overhead
- **Privacy-safe**: No hardcoded paths, no personal data, no API keys

## Testing

Functional tests are inline in each tool file. Run:

```bash
python3 -c "
from tools.medical_tools import med_stats, med_evidence, med_pubmed
# Run a quick t-test
r = med_stats(test='ttest', a=[1,2,3,4,5], b=[6,7,8,9,10])
print(r)
"
```
