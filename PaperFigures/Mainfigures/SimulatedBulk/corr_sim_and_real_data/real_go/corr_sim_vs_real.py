#!/usr/bin/env python3
"""
Compare learnability (AUROC) of real GO terms between simulated bulk data
(all cell types) and real bulk data.

Sim data  : EGADSimulatedGOTerms/data/realgo_allcells/master_melted_df.csv.gz
Real data : bulkEGADPipeline/results/data/splitOPs1/.../Brain_split.csv_melted_EGADs.csv.gz

Output: results/corr_sim_vs_real_realgo.{pdf,png}
        results/corr_sim_vs_real_realgo_auc06.{pdf,png}
        results/merged_aucs.csv
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

# ── paths ──────────────────────────────────────────────────────────────────────
BASE = "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes"

SIM_DATA_PATH = (
    f"{BASE}/bin/EGADSimulatedGOTerms/data/realgo_allcells/master_melted_df.csv.gz"
)
REAL_DATA_PATH = (
    f"{BASE}/bin/bulkEGADPipeline/results/data/splitOPs1"
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
    # Use median to match go_term_exci_expression.py convention
    agg = df.groupby("index")["auc"].median().reset_index()
    agg.columns = ["go_id", "auc_sim"]
    return agg


def load_real(path):
    """Load bulkEGADPipeline output for real GO labels."""
    df = pd.read_csv(path, index_col=0)
    agg = df.groupby("index")["auc"].mean().reset_index()
    agg.columns = ["go_id", "auc_real"]
    return agg


# ── scatter helper ─────────────────────────────────────────────────────────────
def scatter_plot(df, title, out_stem):
    r, p = pearsonr(df["auc_sim"], df["auc_real"])
    print(f"  {title}: r = {r:.4f},  p = {p:.3g},  n = {len(df)}")

    fig, ax = plt.subplots(figsize=(3, 3))
    ax.scatter(df["auc_sim"], df["auc_real"],
               color="#8DA0CB", s=10, alpha=0.5, linewidths=0)

    m, b = np.polyfit(df["auc_sim"], df["auc_real"], 1)
    x_line = np.linspace(df["auc_sim"].min(), df["auc_sim"].max(), 300)
    ax.plot(x_line, m * x_line + b, color="black", lw=1.5, zorder=5)

    ax.text(0.05, 0.97, f"r = {r:.3f}\np = {p:.2g}",
            transform=ax.transAxes, va="top", fontsize=7)
    ax.set_xlabel("AUC  (Simulated Bulk)", fontsize=8)
    ax.set_ylabel("AUC  (Real Bulk)", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    sns.despine(ax=ax)
    plt.tight_layout()

    plt.savefig(f"{OUT_DIR}/{out_stem}.pdf", transparent=True)
    plt.savefig(f"{OUT_DIR}/{out_stem}.png", dpi=300, transparent=True)
    print(f"  Saved → {OUT_DIR}/{out_stem}.png")
    plt.close()


# ── load & merge ───────────────────────────────────────────────────────────────
print("Loading sim data (real GO terms, all cells) ...")
sim_agg = load_sim(SIM_DATA_PATH)
print(f"  {len(sim_agg)} GO terms (sim data)")

print("Loading real data ...")
real_agg = load_real(REAL_DATA_PATH)
print(f"  {len(real_agg)} GO terms (real data)")

merged = pd.merge(sim_agg, real_agg, on="go_id")
print(f"  {len(merged)} matched GO terms after merge")

merged.to_csv(f"{OUT_DIR}/merged_aucs.csv", index=False)
print(f"Saved → {OUT_DIR}/merged_aucs.csv")


# ── correlations ───────────────────────────────────────────────────────────────
sns.set_theme(style="white")
sns.set_context("paper")

print("\n── All GO terms ──")
scatter_plot(merged, "Real GO Labels", "corr_sim_vs_real_realgo")

high_auc = merged[merged["auc_sim"] >= 0.6]
print(f"\n── Sim AUC ≥ 0.  ({len(high_auc)} terms) ──")
scatter_plot(high_auc, "Real GO Labels  (Sim AUC ≥ 0.6)", "corr_sim_vs_real_realgo_auc06")
