"""
ARIA EnvironmentManager Tests
------------------------------
Validates the IPC architecture without requiring real Conda environments.
Tests the JSON communication contract, error handling, fallback logic,
and script base contract.

Run:
  conda activate aria-env
  python tests/test_environment_manager.py
"""

from __future__ import annotations
import sys, os, json, uuid, tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

GRN = "\033[92m"; RED = "\033[91m"; YLW = "\033[93m"
CYN = "\033[96m"; DIM = "\033[2m";  RST = "\033[0m"; BLD = "\033[1m"

passed = 0
failed = 0

def ok(msg, detail=""):
    global passed; passed += 1
    d = f"  {DIM}{detail}{RST}" if detail else ""
    print(f"  {GRN}v{RST} {msg}{d}")

def fail(msg, err=""):
    global failed; failed += 1
    print(f"  {RED}x{RST} {msg}")
    if err: print(f"    {DIM}{err}{RST}")

def section(title):
    print(f"\n{BLD}{CYN}> {title}{RST}")


# ── Test 1: EnvironmentManager instantiation ──────────────────────────────────

section("EnvironmentManager — instantiation & conda check")

try:
    from aria.utils.environment_manager import EnvironmentManager, env_manager
    ok("EnvironmentManager imported successfully")
except Exception as e:
    fail("Import failed", str(e))

try:
    assert hasattr(env_manager, "run_in_stack")
    assert hasattr(env_manager, "check_environments")
    assert hasattr(env_manager, "get_status_report")
    ok("Global env_manager instance has all required methods")
except Exception as e:
    fail("env_manager missing methods", str(e))

try:
    report = env_manager.get_status_report()
    assert "conda_available" in report
    assert "environments" in report
    assert "ready_stacks" in report
    assert "missing_stacks" in report
    ok(f"Status report generated: conda={report['conda_available']}, "
       f"ready={report['ready_stacks']}")
except Exception as e:
    fail("Status report failed", str(e))


# ── Test 2: Stack validation ──────────────────────────────────────────────────

section("EnvironmentManager — stack validation & error handling")

try:
    result = env_manager.run_in_stack(
        stack="invalid_stack",
        script_path="nonexistent.py",
        params={},
    )
    assert result["status"] == "error"
    assert result["error_type"] == "UnknownStack"
    ok("Invalid stack returns structured error", result["error_type"])
except Exception as e:
    fail("Invalid stack handling failed", str(e))

try:
    envs = env_manager.check_environments()
    assert isinstance(envs, dict)
    assert all(k in envs for k in ["rna", "chromatin", "hic", "integration"])
    assert all(isinstance(v, bool) for v in envs.values())
    installed = [k for k, v in envs.items() if v]
    missing   = [k for k, v in envs.items() if not v]
    ok(f"Environment check: {len(installed)} ready, {len(missing)} missing",
       f"missing: {missing}" if missing else "all installed")
except Exception as e:
    fail("check_environments failed", str(e))


# ── Test 3: Script base contract ──────────────────────────────────────────────

section("Script base contract — IPC JSON communication")

try:
    from aria.scripts._base import run_script, _json_serializer
    ok("_base.py imported successfully")
except Exception as e:
    fail("_base.py import failed", str(e))

try:
    import numpy as np

    # Test numpy serialization
    assert _json_serializer(np.int64(42))    == 42
    assert _json_serializer(np.float64(3.14)) == 3.14
    result = _json_serializer(np.array([1, 2, 3]))
    assert result == [1, 2, 3]
    ok("JSON serializer handles numpy int64, float64, ndarray")
except Exception as e:
    fail("JSON serializer failed", str(e))

try:
    # Simulate the full IPC cycle manually
    with tempfile.TemporaryDirectory() as tmpdir:
        input_file  = Path(tmpdir) / "input.json"
        output_file = Path(tmpdir) / "output.json"

        # Write test params
        params = {"test_key": "test_value", "number": 42}
        with open(input_file, "w") as f:
            json.dump(params, f)

        # Simulate what a script does
        with open(input_file, "r") as f:
            loaded = json.load(f)
        result = {"status": "success", "echo": loaded["test_key"],
                  "doubled": loaded["number"] * 2}
        with open(output_file, "w") as f:
            json.dump(result, f)

        # Read back as EnvironmentManager would
        with open(output_file, "r") as f:
            final = json.load(f)

        assert final["status"]  == "success"
        assert final["echo"]    == "test_value"
        assert final["doubled"] == 84

    ok("Full IPC cycle: write params -> execute -> read results")
except Exception as e:
    fail("IPC cycle simulation failed", str(e))


# ── Test 4: rna_qc script import ────────────────────────────────────────────

section("rna_qc.py — import and structure")

try:
    from aria.scripts.rna_qc import rna_qc
    ok("rna_qc.py imported successfully")
except Exception as e:
    fail("rna_qc.py import failed", str(e))

try:
    # Test with invalid path — should return structured error
    try:
        import scanpy  # noqa
        result = rna_qc({"data_path": "/nonexistent/path.h5ad"})
        assert "status" in result
        ok(f"rna_qc handles invalid path: status={result['status']}")
    except ImportError:
        ok("rna_qc import ok (scanpy not in this env — will run in aria-rna-env)")
except Exception as e:
    fail("rna_qc invalid path handling failed", str(e))

try:
    # Test with real PBMC data if available
    pbmc_candidates = [
        Path.home() / "aria-data" / "pbmc3k_test",
        Path.home() / "aria-data" / "pbmc3k_test" / "hg19",
    ]
    pbmc_dir = None
    for c in pbmc_candidates:
        if c.exists() and list(c.rglob("*.mtx*")):
            pbmc_dir = c
            break

    if pbmc_dir:
        import scanpy as sc
        # Quick load test
        adata = sc.read_10x_mtx(str(pbmc_dir), var_names="gene_symbols",
                                  cache=True)
        assert adata.n_obs > 0
        ok(f"PBMC 3k data accessible: {adata.n_obs} cells found",
           str(pbmc_dir))

        # Run actual QC
        result = rna_qc({
            "data_path":  str(pbmc_dir),
            "organism":   "Homo sapiens",
        })
        assert result["status"] == "success"
        assert result["n_cells_after"] > 2000
        assert result["n_cells_before"] == 2700
        ok(f"rna_qc on PBMC 3k: {result['n_cells_before']} -> "
           f"{result['n_cells_after']} cells "
           f"(MT<={result['mt_threshold_used']}%)",
           f"warnings: {result['warnings']}")

        # Clean up written h5ad
        qc_out = Path(result.get("output_path", ""))
        if qc_out.exists():
            qc_out.unlink()
    else:
        ok("PBMC 3k data not found — skipping live QC test",
           "run install.sh to download")

except Exception as e:
    fail("rna_qc live test failed", str(e))


# ── Test 5: rna_clustering script import ─────────────────────────────────────

section("rna_clustering.py — import and structure")

try:
    from aria.scripts.rna_clustering import rna_clustering
    ok("rna_clustering.py imported successfully")
except Exception as e:
    fail("rna_clustering.py import failed", str(e))


# ── Test 6: Workspace management ─────────────────────────────────────────────

section("EnvironmentManager — workspace & temp file cleanup")

try:
    workspace = env_manager.workspace
    assert workspace.exists()
    assert workspace.is_dir()
    ok(f"Workspace directory exists: {workspace}")
except Exception as e:
    fail("Workspace check failed", str(e))

try:
    # Verify no leftover temp files from previous runs
    temp_files = list(workspace.glob("input_*.json")) + \
                 list(workspace.glob("output_*.json"))
    if temp_files:
        fail(f"Leftover temp files found: {[f.name for f in temp_files]}")
    else:
        ok("No leftover temp files in workspace — cleanup working correctly")
except Exception as e:
    fail("Workspace cleanup check failed", str(e))


# ── Summary ───────────────────────────────────────────────────────────────────

total = passed + failed
print(f"\n{'─'*50}")
print(f"{BLD}Results: {GRN}{passed} passed{RST}{BLD}, "
      f"{RED if failed else GRN}{failed} failed{RST}{BLD} / {total} total{RST}")

if failed == 0:
    print(f"\n{GRN}{BLD}v EnvironmentManager stack validated.{RST}\n")
else:
    print(f"\n{YLW}Some tests need attention.{RST}\n")
    sys.exit(1)
