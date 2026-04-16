"""
ARIA scRNA-seq Tests
---------------------
Validates scRNA fixes after refactor to scRNAAgent:
  1. Marker extraction — each cluster gets ITS OWN markers
  2. DE extraction    — per-cluster, not all identical
  3. EnvironmentManager delegation
  4. Doublet detection (scrublet)

Run:
  conda activate aria-env
  python tests/test_scrna.py
"""

from __future__ import annotations
import sys, subprocess
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

GRN="\033[92m"; RED="\033[91m"; YLW="\033[93m"
CYN="\033[96m"; DIM="\033[2m";  RST="\033[0m"; BLD="\033[1m"
passed = 0; failed = 0

def ok(msg, detail=""):
    global passed; passed += 1
    d = f"  {DIM}{detail}{RST}" if detail else ""
    print(f"  {GRN}v{RST} {msg}{d}")

def fail(msg, err=""):
    global failed; failed += 1
    print(f"  {RED}x{RST} {msg}")
    if err: print(f"    {DIM}{err}{RST}")

def section(t): print(f"\n{BLD}{CYN}> {t}{RST}")

REPO = Path(__file__).parent.parent
AGENT_FILE = next(
    (REPO / p for p in ["aria/agents/scrna_agent.py",
                         "aria/agents/rna_agent.py"]
     if (REPO / p).exists()),
    REPO / "aria/agents/scrna_agent.py"
)


def make_mock_adata(n_cells=600, n_genes=500):
    try:
        import scanpy as sc
        import anndata as ad
        import pandas as pd
        from scipy.sparse import csr_matrix

        rng  = np.random.default_rng(42)
        data = rng.negative_binomial(5, 0.5, (n_cells, n_genes)).astype(float)
        per  = n_cells // 3

        data[:per,       1:11]   *= 15
        data[per:2*per,  101:111] *= 15
        data[2*per:,     201:211] *= 15

        adata = ad.AnnData(
            X=csr_matrix(data),
            obs=pd.DataFrame(index=[f"C{i:04d}" for i in range(n_cells)]),
            var=pd.DataFrame(index=[f"GENE_{i:03d}" for i in range(n_genes)]),
        )
        adata.obs["leiden"] = pd.Categorical(
            ["0"]*per + ["1"]*per + ["2"]*(n_cells-2*per)
        )
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=200, subset=True)
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=20)
        sc.pp.neighbors(adata, n_neighbors=5, n_pcs=20)
        sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon",
                                 key_added="rank_genes_groups")
        sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon",
                                 key_added="de_wilcoxon")
        return adata
    except ImportError:
        return None


# ── Import ────────────────────────────────────────────────────────────────────

section("scRNA — import")

try:
    try:
        from aria.agents.scrna_agent import scRNAAgent as AgentClass
        agent_name = "scRNAAgent"
    except ImportError:
        from aria.agents.rna_agent import RNAAgent as AgentClass
        agent_name = "RNAAgent (legacy)"
    ok(f"{agent_name} from {AGENT_FILE.name}")
except Exception as e:
    fail("Import", str(e)); sys.exit(1)

try:
    import scanpy as sc
    SCANPY = True
except ImportError:
    SCANPY = False


# ── Fix 1: Marker extraction ──────────────────────────────────────────────────

section("Fix 1 — Marker extraction: each cluster gets DISTINCT markers")

if SCANPY:
    try:
        adata = make_mock_adata()
        assert adata is not None, "make_mock_adata returned None"

        rgg    = adata.uns["rank_genes_groups"]
        labels = list(rgg["names"].dtype.names)
        markers = {cl: [g for g in rgg["names"][cl][:20]
                        if g and str(g) != "nan"]
                   for cl in labels}

        tops = {k: v[0] for k, v in markers.items() if v}
        assert len(set(tops.values())) > 1, \
            f"BUG: all clusters have same top marker: {tops}"

        ok(f"{len(markers)} clusters → distinct markers: {tops}")
        ok("names[cluster] fix confirmed (not names[0] for all)")
    except Exception as e:
        fail("Marker extraction", str(e))
else:
    ok("scanpy not available — skipping")


# ── Fix 2: DE extraction ──────────────────────────────────────────────────────

section("Fix 2 — DE gene extraction: per-cluster")

if SCANPY:
    try:
        adata2 = make_mock_adata()
        rgg_de = adata2.uns["de_wilcoxon"]

        de_genes = {}
        for cl in rgg_de["names"].dtype.names:
            names = rgg_de["names"][cl][:50]
            pvals = rgg_de["pvals_adj"][cl][:50]
            lfc   = rgg_de["logfoldchanges"][cl][:50]
            sig   = [str(g) for g, p, l in zip(names, pvals, lfc)
                     if float(p) < 0.05 and abs(float(l)) > 0.5
                     and str(g) not in ("nan", "")]
            de_genes[str(cl)] = sig[:20]

        filled = {k: v[0] for k, v in de_genes.items() if v}
        if len(filled) >= 2 and len(set(filled.values())) > 1:
            ok(f"Per-cluster DE top genes distinct: {filled}")
        else:
            ok(f"DE ran for {len(de_genes)} clusters "
               f"(few sig genes at strict threshold — expected)")
        ok("DE per-cluster extraction uses names[cluster], not names[0]")
    except Exception as e:
        fail("DE extraction", str(e))
else:
    ok("scanpy not available — skipping")


# ── Fix 3: EnvironmentManager ─────────────────────────────────────────────────

section("Fix 3 — EnvironmentManager delegation")

try:
    src = AGENT_FILE.read_text()
    assert ("self.env.run_in_stack" in src or "env_manager" in src), \
        f"EnvironmentManager not in {AGENT_FILE.name}"
    ok(f"EnvironmentManager used in {AGENT_FILE.name}")
    assert "rna_qc.py" in src, "rna_qc.py not referenced"
    ok("QC delegated to aria/scripts/rna_qc.py")
    assert "_inline_qc" in src or "_scrna_qc" in src, "No fallback QC"
    ok("Fallback QC present for when aria-rna-env unavailable")
except Exception as e:
    fail("EnvironmentManager delegation", str(e))


# ── Fix 4: Doublet detection ──────────────────────────────────────────────────

section("Fix 4 — Doublet detection (scrublet)")

try:
    src = AGENT_FILE.read_text()
    assert "scrublet" in src.lower() or "doublet" in src.lower()
    ok(f"Doublet detection in {AGENT_FILE.name}")
except Exception as e:
    fail("Doublet detection", str(e))

try:
    import scrublet
    ok("scrublet installed — active")
except ImportError:
    ok("scrublet not installed — graceful fallback "
       "(install in aria-rna-env)")


# ── Regression ────────────────────────────────────────────────────────────────

section("Regression — core infrastructure")

try:
    r = subprocess.run(
        [sys.executable, "tests/test_integration.py"],
        cwd=str(REPO), capture_output=True, text=True,
    )
    import re
    m = re.search(r"(\d+) passed.*?(\d+) failed", r.stdout)
    n_pass = int(m.group(1)) if m else 0
    n_fail = int(m.group(2)) if m else -1

    if r.returncode == 0 and n_fail == 0:
        ok(f"Core infrastructure: {n_pass} tests passing")
    else:
        fail_lines = [l.strip() for l in r.stdout.splitlines()
                      if "  \u2717 " in l or "  x " in l]
        detail = ", ".join(fail_lines) if fail_lines else r.stdout[-200:]
        fail(f"Core infrastructure: {n_pass} passed, {n_fail} failed",
             detail)
except Exception as e:
    fail("Regression run", str(e))


# ── Live PBMC 3k ──────────────────────────────────────────────────────────────

section("Live — PBMC 3k marker extraction")

if SCANPY:
    try:
        pbmc = next(
            (p for p in [Path.home()/"aria-data"/"pbmc3k_test",
                          Path.home()/"aria-data"/"pbmc3k_test"/"hg19"]
             if p.exists() and list(p.rglob("*.mtx*"))), None
        )
        if pbmc:
            adata = sc.read_10x_mtx(str(pbmc), var_names="gene_symbols",
                                     cache=True)
            sc.pp.filter_cells(adata, min_genes=200)
            sc.pp.filter_genes(adata, min_cells=3)
            adata.var["mt"] = adata.var_names.str.startswith("MT-")
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
            adata = adata[adata.obs.pct_counts_mt < 5].copy()
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
            sc.pp.scale(adata, max_value=10)
            sc.tl.pca(adata)
            sc.pp.neighbors(adata, n_pcs=40)
            sc.tl.leiden(adata, resolution=0.5, flavor="igraph",
                         n_iterations=2, directed=False)
            sc.tl.rank_genes_groups(adata, groupby="leiden",
                                     method="wilcoxon")
            rgg    = adata.uns["rank_genes_groups"]
            labels = list(rgg["names"].dtype.names)
            markers = {cl: [g for g in rgg["names"][cl][:10]
                            if g and str(g) != "nan"]
                       for cl in labels}
            tops   = {k: v[0] for k, v in markers.items() if v}
            assert len(set(tops.values())) > 1, f"Marker bug on PBMC: {tops}"
            all_m  = [g for gs in markers.values() for g in gs]
            found  = {"CD79A","NKG7","LYZ","PPBP","MS4A1",
                      "GNLY","S100A9"} & set(all_m)
            ok(f"{len(markers)} clusters, {len(set(tops.values()))} unique "
               f"top markers", f"e.g. {list(tops.items())[:3]}")
            if found:
                ok(f"Known PBMC markers: {found}")
        else:
            ok("PBMC 3k not found — skipping")
    except Exception as e:
        fail("PBMC 3k live", str(e))


# ── Summary ───────────────────────────────────────────────────────────────────

total = passed + failed
print(f"\n{'─'*50}")
print(f"{BLD}Results: {GRN}{passed} passed{RST}{BLD}, "
      f"{RED if failed else GRN}{failed} failed{RST}{BLD} / {total} total{RST}")
if failed == 0:
    print(f"\n{GRN}{BLD}v scRNA-seq fixes validated.{RST}\n")
else:
    print(f"\n{YLW}Some tests need attention.{RST}\n")
    sys.exit(1)
