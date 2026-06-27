#!/usr/bin/env python3
"""
Compare learnability (AUROC) of simulated excitatory-neuron-label GO terms
between simulated bulk data (all cells, including excitatory neurons) and
real bulk data.

Sim data  : EGADSimulatedGOTerms/data/simgo_allcells/master_melted_df.csv.gz
Real data : bulkEGADPipeline/results/realbulk_simexilabels/.../Brain_split.csv_melted_EGADs.csv.gz

Output: results/corr_sim_vs_real_exi.{pdf,png}, results/merged_aucs.csv
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# ── paths ──────────────────────────────────────────────────────────────────────
BASE = "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes"

SIM_DATA_PATH = (
    f"{BASE}/bin/EGADSimulatedGOTerms/data/simgo_allcells/master_melted_df.csv.gz"
)
REAL_DATA_PATH = (
    f"{BASE}/bin/bulkEGADPipeline/results/realbulk_simexilabels"
    "/EGAD/melted_dfs/Brain_split.csv_melted_EGADs.csv.gz"
)
OUT_DIR = "results"
os.makedirs(OUT_DIR, exist_ok=True)


# ── loaders ────────────────────────────────────────────────────────────────────
def load_sim(path):
    """Load EGADSimulatedGOTerms output, filter to experimental conditions."""
    df = pd.read_csv(path, index_col=0)
    if "organism_part" in df.columns:
        df = df[df["organism_part"].str.split("_").str[0] == "exp"]
    df["Q"] = df["index"].str.split("_").str[-1]
    df = df[df["Q"].isin(["Q2", "Q3"])]
    agg = df.groupby(["index", "Q"])["auc"].mean().reset_index()
    agg.columns = ["go_id", "Q", "auc_sim"]
    return agg


def load_real(path):
    """Load bulkEGADPipeline output for simulated excitatory-neuron labels."""
    df = pd.read_csv(path, index_col=0)
    df["Q"] = df["index"].str.split("_").str[-1]
    df = df[df["Q"].isin(["Q2", "Q3"])]
    agg = df.groupby(["index", "Q"])["auc"].mean().reset_index()
    agg.columns = ["go_id", "Q", "auc_real"]
    return agg


# ── load & merge ───────────────────────────────────────────────────────────────
print("Loading sim data (all cells, excitatory neurons included) ...")
sim_agg = load_sim(SIM_DATA_PATH)
print(f"  {sim_agg['go_id'].nunique()} SIMGO terms, {len(sim_agg)} rows")

print("Loading real data ...")
real_agg = load_real(REAL_DATA_PATH)
print(f"  {real_agg['go_id'].nunique()} SIMGO terms, {len(real_agg)} rows")

merged = pd.merge(sim_agg, real_agg, on=["go_id", "Q"])
print(f"  {len(merged)} matched rows after merge")

merged.to_csv(f"{OUT_DIR}/merged_aucs.csv", index=False)
print(f"Saved → {OUT_DIR}/merged_aucs.csv")


# ── Pearson correlation (per Q group) ─────────────────────────────────────────
corr = {}
for q, grp in merged.groupby("Q"):
    r, p = pearsonr(grp["auc_sim"], grp["auc_real"])
    corr[q] = (r, p)
    print(f"  {q}: r = {r:.4f},  p = {p:.3g},  n = {len(grp)}")


# ── scatter plots (one panel per Q group) ─────────────────────────────────────
sns.set_theme(style="white")
sns.set_context("paper")

palette  = {"Q2": "#66C2A5", "Q3": "#FC8D62"}
q_labels = {"Q2": "Non-enriched (Q2)", "Q3": "Exci-enriched (Q3)"}

fig, axes = plt.subplots(1, 2, figsize=(5.5, 2.8))

for ax, q in zip(axes, ["Q2", "Q3"]):
    grp = merged[merged["Q"] == q]
    r, p = corr[q]

    ax.scatter(
        grp["auc_sim"], grp["auc_real"],
        color=palette[q], s=8, alpha=0.5, linewidths=0,
    )

    m, b = np.polyfit(grp["auc_sim"], grp["auc_real"], 1)
    x_line = np.linspace(grp["auc_sim"].min(), grp["auc_sim"].max(), 300)
    ax.plot(x_line, m * x_line + b, color="black", lw=1.5, zorder=5)

    ax.text(
        0.05, 0.97,
        f"r = {r:.3f}\np = {p:.2g}",
        transform=ax.transAxes, va="top", fontsize=7,
    )

    ax.set_xlabel("AUC  (Simulated Bulk)", fontsize=8)
    ax.set_ylabel("AUC  (Real Bulk)", fontsize=8)
    ax.set_title(q_labels[q], fontsize=8)
    ax.tick_params(labelsize=7)
    sns.despine(ax=ax)

fig.suptitle("Sim Excitatory-Neuron Labels", fontsize=9, y=1.01)
plt.tight_layout()

plt.savefig(f"{OUT_DIR}/corr_sim_vs_real_exi.pdf", transparent=True)
plt.savefig(f"{OUT_DIR}/corr_sim_vs_real_exi.png", dpi=300, transparent=True)
print(f"Saved → {OUT_DIR}/corr_sim_vs_real_exi.png")
plt.show()
