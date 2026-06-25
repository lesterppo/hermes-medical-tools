#!/usr/bin/env python3
"""
PSPP integration — native, token-efficient statistical analysis tool for Hermes.

PSPP is a free replacement for IBM SPSS, supporting 135 commands across
descriptive statistics, parametric/non-parametric tests, regression,
factor analysis, reliability, ROC curves, and more.

Design: single `pspp` tool. Model writes PSPP/SPSS syntax (already in
training data). Tool executes via subprocess, captures CSV output,
parses into compact JSON. Zero learning curve — model just writes
the SPSS syntax it already knows.

Output format:
- Each PSPP command produces a "Table: Name" section in the output.
- Tables are parsed from CSV into compact JSON: {cols:[...], rows:[...]}
- For raw text, pass out="text"
"""

from __future__ import annotations

import csv
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# ═══════════════════════════════════════════════════════════════════
# Dependency check
# ═══════════════════════════════════════════════════════════════════

def _pspp_available() -> bool:
    return shutil.which("pspp") is not None


# ═══════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════

def _ok(result: dict) -> str:
    return json.dumps(result, ensure_ascii=False, separators=(",", ":"))


def _err(msg: str) -> str:
    return _ok({"e": msg})


def _parse_csv_tables(output: str) -> List[Dict]:
    """Parse PSPP CSV output into structured tables.

    PSPP outputs CSV sections delimited by 'Table: <Name>' headers.
    Each table has a header row followed by data rows.
    """
    tables = []
    current_name = None
    current_lines = []

    for line in output.split("\n"):
        line = line.rstrip("\r")
        if line.startswith("Table:"):
            if current_name and current_lines:
                tables.append(_table_from_lines(current_name, current_lines))
                current_lines = []
            current_name = line[len("Table:"):].strip()
        elif current_name and line.strip():
            current_lines.append(line)

    if current_name and current_lines:
        tables.append(_table_from_lines(current_name, current_lines))

    return tables


def _table_from_lines(name: str, lines: List[str]) -> Dict:
    """Parse CSV lines into dict with cols and rows."""
    try:
        reader = csv.reader(io.StringIO("\n".join(lines)))
        rows = list(reader)
        if not rows:
            return {"n": name, "cols": [], "rows": []}

        # Filter empty trailing rows
        while rows and all(c.strip() == "" for c in rows[-1]):
            rows.pop()

        if not rows:
            return {"n": name, "cols": [], "rows": []}

        cols = rows[0]
        data_rows = rows[1:]

        # Truncate long rows for token efficiency
        if len(data_rows) > 50:
            data_rows = data_rows[:50]

        return {"n": name, "cols": cols, "rows": data_rows}
    except Exception:
        return {"n": name, "cols": [], "rows": [], "e": "parse error"}


def _compact_numbers(row: List[str]) -> List[Any]:
    """Convert numeric strings to numbers for compact JSON."""
    result = []
    for cell in row:
        stripped = cell.strip()
        if not stripped:
            result.append("")
            continue
        try:
            if "." in stripped:
                result.append(round(float(stripped), 4))
            else:
                result.append(int(stripped))
        except ValueError:
            result.append(stripped)
    return result


# ═══════════════════════════════════════════════════════════════════
# Tool: pspp — Statistical analysis via PSPP
# ═══════════════════════════════════════════════════════════════════

def pspp_run(
    syntax: str,
    data: Optional[str] = None,
    data_file: Optional[str] = None,
    out: str = "json",
    include_syntax: bool = False,
) -> str:
    """
    Run PSPP statistical analysis. Returns compact JSON tables by default.

    syntax:          PSPP/SPSS commands. Model writes the SPSS syntax it knows.
    data:            Inline data in PSPP format (DATA LIST ... BEGIN DATA ... END DATA).
                     If provided, prepended to syntax.
    data_file:       Path to a data file (.sav, .csv, .por, .ods, .gnumeric, .txt).
                     Use INSTEAD of data for large files or existing datasets.
                     For CSV: use "GET DATA /TYPE=TXT /FILE='path.csv' /DELIMITERS=','."
    out:             "json" (parsed tables) | "text" (raw output) | "csv" (raw CSV)
    include_syntax:  If true, include the syntax in output (for debugging)
    """
    try:
        # Build PSPP input
        lines = []

        # Load data file if specified
        if data_file:
            if not os.path.exists(data_file):
                return _err(f"File not found: {data_file}")
            ext = Path(data_file).suffix.lower()
            if ext == ".sav":
                lines.append(f"GET FILE='{data_file}'.")
            elif ext == ".csv":
                lines.append("GET DATA /TYPE=TXT")
                lines.append(f" /FILE='{data_file}'")
                lines.append(" /DELIMITERS=',' /QUALIFIER='\"' /FIRSTCASE=2.")
            elif ext in (".por", ".sps"):
                lines.append(f"GET FILE='{data_file}'.")
            elif ext in (".ods", ".gnumeric"):
                lines.append(f"GET DATA /TYPE={'ODS' if ext=='.ods' else 'GNM'}")
                lines.append(f" /FILE='{data_file}' /SHEET=1.")
            else:
                lines.append("GET DATA /TYPE=TXT")
                lines.append(f" /FILE='{data_file}' /DELIMITERS='\\t' /FIRSTCASE=2.")
            lines.append("")

        # Add inline data if provided
        if data:
            lines.append(data.strip())
            lines.append("")

        # Add syntax
        lines.append(syntax.strip())

        pspp_input = "\n".join(lines) + "\n"

        # Write syntax to temp file (PSPP needs a file for reliable output)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".sps", delete=False, encoding="utf-8"
        ) as f:
            f.write(pspp_input)
            sps_path = f.name

        try:
            # Run PSPP: suppress default ASCII driver, output CSV to stdout
            cmd = [
                "pspp",
                "--no-statrc",
                "--no-output",              # suppress default terminal driver
                "-o", "/dev/stdout",         # CSV output to stdout
                "-O", "format=csv",          # CSV format
                sps_path,
            ]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=60,
                env={**os.environ, "LC_ALL": "C.UTF-8"},
            )

            stderr = result.stderr.strip()
            stdout = result.stdout

            # Check for errors
            if result.returncode != 0 and not stdout:
                return _err(f"PSPP error (exit {result.returncode}): {stderr[:300]}")

            if out == "text":
                return _ok({"out": "text", "t": stdout[:8000]})

            if out == "csv":
                return _ok({"out": "csv", "t": stdout[:8000]})

            # Parse tables into compact JSON
            tables = _parse_csv_tables(stdout)

            if not tables:
                # Maybe output went to stderr only
                if not stdout.strip():
                    return _err(f"PSPP produced no output. Errors: {stderr[:300]}")
                return _ok({"out": "json", "tables": [], "raw": stdout[:500]})

            # Compact the tables
            compact = []
            for t in tables:
                entry = {"n": t["n"], "c": t["cols"], "r": []}
                for row in t["rows"]:
                    entry["r"].append(_compact_numbers(row))
                compact.append(entry)

            result_dict: Dict[str, Any] = {
                "out": "json",
                "tables": compact,
                "n_tables": len(compact),
            }

            # Include stderr warnings as a note
            if stderr:
                # Filter common non-error messages
                warnings = [l for l in stderr.split("\n")
                           if l.strip() and "warning" in l.lower()]
                if warnings:
                    result_dict["warnings"] = warnings[:5]

            if include_syntax:
                result_dict["syntax"] = pspp_input[:1000]

            return _ok(result_dict)

        finally:
            os.unlink(sps_path)

    except subprocess.TimeoutExpired:
        return _err("PSPP timed out (60s). Simplify syntax or reduce data.")
    except Exception as e:
        return _err(f"PSPP error: {e}")


# ═══════════════════════════════════════════════════════════════════
# Schema (token-efficient)
# ═══════════════════════════════════════════════════════════════════

PSPP_SCHEMA = {
    "name": "pspp",
    "description": (
        "Statistical analysis via GNU PSPP (free SPSS replacement). "
        "Write PSPP/SPSS syntax for DESCRIPTIVES, T-TEST, CROSSTABS, REGRESSION, "
        "FREQUENCIES, EXAMINE, CORRELATIONS, ONEWAY, GLM, FACTOR, RELIABILITY, "
        "ROC, NPAR TESTS, LOGISTIC REGRESSION, MEANS, RANK, GRAPH, etc. "
        "Returns compact JSON tables parsed from CSV output. "
        "For inline data, pass data='DATA LIST... BEGIN DATA... END DATA'. "
        "For files, pass data_file='/path/to/file.sav' or .csv. "
        "out='json' (parsed, default), 'text' (raw), or 'csv'."
    ),
    "parameters": {
        "type": "object",
        "properties": {
            "syntax": {
                "type": "string",
                "description": (
                    "PSPP/SPSS commands. E.g. 'DESCRIPTIVES var1 var2.' or "
                    "'T-TEST GROUPS=group(1 2) /VARIABLES=score.' or "
                    "'CROSSTABS var1 BY var2 /STATISTICS=CHISQ.' "
                    "Model knows SPSS syntax from training — write it naturally."
                ),
            },
            "data": {
                "type": "string",
                "description": (
                    "Inline data block: 'DATA LIST FREE /var1 var2. BEGIN DATA 1 2 3 4 END DATA.' "
                    "Use for quick analyses with small inline data."
                ),
            },
            "data_file": {
                "type": "string",
                "description": (
                    "Path to .sav, .csv, .por, .ods, .gnumeric file. "
                    "Use INSTEAD of data for large files. Auto-detects format by extension."
                ),
            },
            "out": {
                "type": "string",
                "description": "'json' (parsed tables, default), 'text' (raw), 'csv' (raw CSV)",
            },
            "include_syntax": {
                "type": "boolean",
                "description": "Echo back syntax in response (for debugging)",
            },
        },
        "required": ["syntax"],
    },
}


# ═══════════════════════════════════════════════════════════════════
# Registry
# ═══════════════════════════════════════════════════════════════════

from tools.registry import registry

registry.register(
    name="pspp",
    toolset="medical",
    schema=PSPP_SCHEMA,
    handler=lambda args, **kw: pspp_run(
        syntax=args.get("syntax", ""),
        data=args.get("data"),
        data_file=args.get("data_file"),
        out=args.get("out", "json"),
        include_syntax=args.get("include_syntax", False),
    ),
    check_fn=_pspp_available,
    emoji="📊",
)
