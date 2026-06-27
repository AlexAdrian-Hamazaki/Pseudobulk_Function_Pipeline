#!/usr/bin/env python3
"""
Excitatory neuron marker genes via scanpy t-test with variance overestimation.

More formal than Wilcoxon because:
  - Assumes log-normal distribution (appropriate for log1p-normalised data)
  - Overestimates within-group variance (Bayesian-style shrinkage toward a
    shared variance, similar in spirit to limma's moderated t-test) → more
    conservative / fewer false positives than the standard t-test
  - Reports actual mean-based log-fold-changes and per-group detection rates

Outputs
-------
  results/mast_markers_excitatory.csv   ranked marker table
  results/volcano_excitatory.pdf/.png   volcano plot
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse

# ── config ────────────────────────────────────────────────────────────────────
H5AD = (
    "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes"
    "/bin/preprocessing/preprocessSC_pipe/data/h5ad_datasets"
    "/processed/cpm/brain_sc_with_metadata_pc_cpm.h5ad"
)
GROUP         = "Excitatory neurons"
CELL_TYPE_COL = "Cell type"
OUT_CSV       = "results/mast_markers_excitatory.csv"
OUT_MARKERS   = "results/mast_markers_excitatory_sig.csv"
OUT_PDF       = "results/volcano_excitatory.pdf"
OUT_PNG       = "results/volcano_excitatory.png"

MIN_CELLS         = 100
MIN_FRAC_IN_GROUP = 0.1

# significance / effect-size thresholds for colouring the volcano
PADJ_THRESH = 0.01
LFC_THRESH  = 2.0


# ── load & preprocess ─────────────────────────────────────────────────────────
print("Loading data...")
adata = ad.read_h5ad(H5AD)
print(f"  {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ── gene filter (same thresholds as notebook) ─────────────────────────────────
print("Filtering genes...")
group_mask = adata.obs[CELL_TYPE_COL] == GROUP
X = adata.X

if sparse.issparse(X):
    n_expr   = np.asarray(X.astype(bool).sum(axis=0)).ravel()
    frac_grp = np.asarray(X[group_mask].astype(bool).sum(axis=0)).ravel() / group_mask.sum()
else:
    n_expr   = (X > 0).sum(axis=0)
    frac_grp = (X[group_mask] > 0).sum(axis=0) / group_mask.sum()

keep = (n_expr >= MIN_CELLS) | (frac_grp >= MIN_FRAC_IN_GROUP)
adata = adata[:, keep].copy()
print(f"  Kept {keep.sum():,} of {len(keep):,} genes")

# ── one-vs-all label ──────────────────────────────────────────────────────────
adata.obs["one_vs_all"] = adata.obs[CELL_TYPE_COL].eq(GROUP).astype(str)

# ── differential expression ───────────────────────────────────────────────────
print("Running t-test with variance overestimation (scanpy)...")
sc.tl.rank_genes_groups(
    adata,
    groupby="one_vs_all",
    method="t-test_overestim_var",   # more formal than Wilcoxon for log-normal data
    reference="rest",
    pts=True,                        # report % cells expressing per group
    corr_method="benjamini-hochberg",
)

deg = sc.get.rank_genes_groups_df(adata, group="True")
deg.rename(columns={"pct_nz_group": "pct_group", "pct_nz_reference": "pct_rest"},
           inplace=True, errors="ignore")
deg.to_csv(OUT_CSV, index=False)
print(f"Saved full results → {OUT_CSV}")

sig = deg[(deg["pvals_adj"] <= PADJ_THRESH) & (deg["logfoldchanges"] > LFC_THRESH)]
sig.to_csv(OUT_MARKERS, index=False)
print(f"Significant markers (adj p ≤ {PADJ_THRESH}, logFC > {LFC_THRESH}): {len(sig):,}")
print(f"Saved markers → {OUT_MARKERS}")
print(sig.head(10).to_string(index=False))


# ── volcano plot ──────────────────────────────────────────────────────────────
print("Plotting volcano...")

x = deg["logfoldchanges"]
y = -np.log10(deg["pvals_adj"].clip(lower=1e-300))

# colour by significance category
def _category(row):
    if row["pvals_adj"] <= PADJ_THRESH and row["logfoldchanges"] > LFC_THRESH:
        return "up"
    if row["pvals_adj"] <= PADJ_THRESH and row["logfoldchanges"] < -LFC_THRESH:
        return "down"
    return "ns"

deg["category"] = deg.apply(_category, axis=1)
palette = {"up": "#C0392B", "down": "#2980B9", "ns": "#BDC3C7"}

sns.set_theme(style="white")
sns.set_context("paper")
fig, ax = plt.subplots(figsize=(4, 4))

for cat, colour in palette.items():
    mask = deg["category"] == cat
    ax.scatter(x[mask], y[mask], s=4, color=colour, alpha=0.6, linewidths=0, label=cat)

# threshold lines
ax.axhline(-np.log10(PADJ_THRESH), color="black", lw=0.8, ls="--")
ax.axvline( LFC_THRESH,            color="black", lw=0.8, ls="--")
ax.axvline(-LFC_THRESH,            color="black", lw=0.8, ls="--")

# annotate top 10 up-regulated genes
top = deg[deg["category"] == "up"].head(10)
for _, row in top.iterrows():
    xi = row["logfoldchanges"]
    yi = -np.log10(max(row["pvals_adj"], 1e-300))
    ax.text(xi + 0.05, yi, row["names"], fontsize=5, va="center")

counts = deg["category"].value_counts()
ax.legend(
    title="",
    labels=[f"up ({counts.get('up', 0):,})",
            f"ns ({counts.get('ns', 0):,})",
            f"down ({counts.get('down', 0):,})"],
    markerscale=2, fontsize=7,
)

ax.set_xlabel("Log2 fold change (Excitatory vs rest)", fontsize=9)
ax.set_ylabel(r"$-\log_{10}$(adj p-value)",            fontsize=9)
ax.set_title("Excitatory neurons — t-test overestim_var", fontsize=9)
sns.despine(ax=ax)
plt.tight_layout()

fig.savefig(OUT_PDF, transparent=True)
fig.savefig(OUT_PNG, dpi=300, transparent=True)
print(f"Volcano saved → {OUT_PDF} / {OUT_PNG}")
