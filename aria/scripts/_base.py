"""
ARIA Analysis Script Base Contract
------------------------------------
Every script executed by EnvironmentManager MUST follow this pattern.

Interface:
    python script.py <input_json_path> <output_json_path>

Input JSON:  parameters dict written by the calling agent
Output JSON: results dict read by the calling agent

The output dict MUST always contain:
    {"status": "success" | "error", ...additional keys...}

On error, MUST contain:
    {"status": "error", "error_type": str, "details": str}

This contract ensures that EnvironmentManager can always parse the result
regardless of which script produced it.

Usage (in your analysis script):
    from aria.scripts._base import run_script

    def my_analysis(params: dict) -> dict:
        # do work
        return {"n_cells": 2643, "qc_passed": True}

    if __name__ == "__main__":
        run_script(my_analysis)
"""

from __future__ import annotations

import json
import sys
import traceback
import warnings
from pathlib import Path
from typing import Callable


def run_script(analysis_fn: Callable[[dict], dict]) -> None:
    """
    Standard entry point for all ARIA analysis scripts.

    Reads params from input JSON, calls analysis_fn, writes results
    to output JSON. Catches all exceptions and writes them as structured
    error responses so the agent always gets parseable JSON back.

    Args:
        analysis_fn: Function that takes params dict, returns results dict.
                     Must include {"status": "success"} in the return value
                     or the base will add it automatically.
    """
    # Silence common bioinformatics warnings that pollute stderr
    # These are safe to suppress for production runs
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning,
                            message=".*sparse.*")
    warnings.filterwarnings("ignore", category=RuntimeWarning,
                            message=".*log2.*")

    if len(sys.argv) != 3:
        _write_error(
            output_path=None,
            error_type="InvalidArguments",
            details=(
                f"Expected: python {sys.argv[0]} <input.json> <output.json>\n"
                f"Got: {' '.join(sys.argv)}"
            ),
        )
        sys.exit(1)

    input_path  = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    # Read input
    try:
        with open(input_path, "r") as f:
            params = json.load(f)
    except Exception as e:
        _write_error(output_path, "InputReadError", str(e))
        sys.exit(1)

    # Execute analysis
    try:
        result = analysis_fn(params)

        # Ensure status is always present
        if "status" not in result:
            result["status"] = "success"

        with open(output_path, "w") as f:
            json.dump(result, f, default=_json_serializer)

    except Exception as e:
        _write_error(
            output_path,
            error_type=type(e).__name__,
            details=f"{str(e)}\n\n{traceback.format_exc()}",
        )
        sys.exit(1)


def _write_error(output_path: Path | None,
                 error_type: str,
                 details: str) -> None:
    """Write a structured error response to output JSON."""
    payload = {
        "status":     "error",
        "error_type": error_type,
        "details":    details,
    }
    if output_path is not None:
        try:
            with open(output_path, "w") as f:
                json.dump(payload, f)
        except Exception:
            # Last resort: print to stderr
            print(json.dumps(payload), file=sys.stderr)
    else:
        print(json.dumps(payload), file=sys.stderr)


def _json_serializer(obj):
    """
    Handle numpy and other non-JSON-serializable types.
    Called by json.dump as the default= handler.
    """
    # numpy scalar types
    try:
        import numpy as np
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
    except ImportError:
        pass

    # pandas types
    try:
        import pandas as pd
        if isinstance(obj, pd.Series):
            return obj.tolist()
        if isinstance(obj, pd.DataFrame):
            return obj.to_dict(orient="list")
    except ImportError:
        pass

    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")
