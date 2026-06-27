#!/usr/bin/env python3
"""
Compare learnability (AUROC) of simulated brain-label GO terms
between simulated bulk data and real bulk data.

Sim data  : EGADSimulatedGOTerms/data/sim_go_apr23_brain/master_melted_df.csv.gz
Real data : bulkEGADPipeline/results/realbulk_simbrainlabels/.../Brain_split.csv_melted_EGADs.csv.gz
MGC data  : analysis/PaperFigs/SimulatedBulk/sim_functions_sim_data/brain/data/simulated_brain_GOs.csv

Output: results/corr_sim_vs_real_brain.{pdf,png}
        results/corr_sim_vs_real_brain_mgcbins.{pdf,png}
        results/merged_aucs.csv
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
    f"{BASE}/bin/EGADSimulatedGOTerms/data/sim_go_apr23_brain/master_melted_df.csv.gz"
)
REAL_DATA_PATH = (
    f"{BASE}/bin/bulkEGADPipeline/results/realbulk_simbrainlabels"
    "/EGAD/melted_dfs/Brain_split.csv_melted_EGADs.csv.gz"
)
SIMGO_LABELS_PATH = (
    f"{BASE}/analysis/PaperFigs/SimulatedBulk/sim_functions_sim_data"
    "/brain/data/simulated_brain_GOs.csv"
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
    """Load bulkEGADPipeline output for simulated brain labels."""
    df = pd.read_csv(path, index_col=0)
    df["Q"] = df["index"].str.split("_").str[-1]
    df = df[df["Q"].isin(["Q2", "Q3"])]
    agg = df.groupby(["index", "Q"])["auc"].mean().reset_index()
    agg.columns = ["go_id", "Q", "auc_real"]
    return agg


def load_mgc_per_simgo(path):
    """Compute mean brain MGC score per SIMGO GO.ID from the gene-level labels file."""
    df = pd.read_csv(path)
    mgc = df.groupby("GO.ID")["score"].mean().reset_index()
    mgc.columns = ["go_id", "mean_mgc"]
    return mgc


def assign_mgc_bins(df, mgc_col="mean_mgc"):
    """Label each row as Bottom 5%, Center 90%, or Top 5% by MGC score."""
    p5  = df[mgc_col].quantile(0.05)
    p95 = df[mgc_col].quantile(0.95)
    conditions = [df[mgc_col] <= p5, df[mgc_col] >= p95]
    choices    = ["Bottom 5%", "Top 5%"]
    df = df.copy()
    df["mgc_bin"] = np.select(conditions, choices, default="Center 90%")
    df["mgc_bin"] = pd.Categorical(
        df["mgc_bin"],
        categories=["Bottom 5%", "Center 90%", "Top 5%"],
        ordered=True,
    )
    return df


# ── load & merge ───────────────────────────────────────────────────────────────
print("Loading sim data ...")
sim_agg = load_sim(SIM_DATA_PATH)
print(f"  {sim_agg['go_id'].nunique()} SIMGO terms, {len(sim_agg)} rows")

print("Loading real data ...")
real_agg = load_real(REAL_DATA_PATH)
print(f"  {real_agg['go_id'].nunique()} SIMGO terms, {len(real_agg)} rows")

merged = pd.merge(sim_agg, real_agg, on=["go_id", "Q"])
print(f"  {len(merged)} matched rows after merge")

print("Loading MGC scores ...")
mgc_per_go = load_mgc_per_simgo(SIMGO_LABELS_PATH)
print(f"  {len(mgc_per_go)} GO.ID entries")

merged = merged.merge(mgc_per_go, on="go_id", how="left")
n_missing = merged["mean_mgc"].isna().sum()
if n_missing:
    print(f"  WARNING: {n_missing} rows had no MGC score; dropping them")
    merged = merged.dropna(subset=["mean_mgc"])

merged = assign_mgc_bins(merged)
print("MGC bin counts:")
print(merged["mgc_bin"].value_counts().sort_index())

merged.to_csv(f"{OUT_DIR}/merged_aucs.csv", index=False)
print(f"Saved → {OUT_DIR}/merged_aucs.csv")


# ── Pearson correlation (per Q group) ─────────────────────────────────────────
print("\nOverall correlations per Q:")
corr = {}
for q, grp in merged.groupby("Q"):
    r, p = pearsonr(grp["auc_sim"], grp["auc_real"])
    corr[q] = (r, p)
    print(f"  {q}: r = {r:.4f},  p = {p:.3g},  n = {len(grp)}")


# ── Pearson correlation per (Q group × MGC bin) ───────────────────────────────
print("\nCorrelations per Q × MGC bin:")
corr_bins = {}
for (q, b), grp in merged.groupby(["Q", "mgc_bin"], observed=True):
    if len(grp) < 3:
        corr_bins[(q, b)] = (float("nan"), float("nan"), len(grp))
        continue
    r, p = pearsonr(grp["auc_sim"], grp["auc_real"])
    corr_bins[(q, b)] = (r, p, len(grp))
    print(f"  {q} | {b}: r = {r:.4f},  p = {p:.3g},  n = {len(grp)}")


# ── scatter plots: overall (one panel per Q) ──────────────────────────────────
sns.set_theme(style="white")
sns.set_context("paper")

palette  = {"Q2": "#66C2A5", "Q3": "#FC8D62"}
q_labels = {"Q2": "Non-enriched (Q2)", "Q3": "Brain-enriched (Q3)"}

fig, axes = plt.subplots(1, 2, figsize=(5.5, 2.8))
for ax, q in zip(axes, ["Q2", "Q3"]):
    grp = merged[merged["Q"] == q]
    r, p = corr[q]
    ax.scatter(grp["auc_sim"], grp["auc_real"],
               color=palette[q], s=8, alpha=0.5, linewidths=0)
    m, b = np.polyfit(grp["auc_sim"], grp["auc_real"], 1)
    x_line = np.linspace(grp["auc_sim"].min(), grp["auc_sim"].max(), 300)
    ax.plot(x_line, m * x_line + b, color="black", lw=1.5, zorder=5)
    ax.text(0.05, 0.97, f"r = {r:.3f}\np = {p:.2g}",
            transform=ax.transAxes, va="top", fontsize=7)
    ax.set_xlabel("AUC  (Simulated Bulk)", fontsize=8)
    ax.set_ylabel("AUC  (Real Bulk)", fontsize=8)
    ax.set_title(q_labels[q], fontsize=8)
    ax.tick_params(labelsize=7)
    sns.despine(ax=ax)

fig.suptitle("Sim Brain Labels", fontsize=9, y=1.01)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/corr_sim_vs_real_brain.pdf", transparent=True)
plt.savefig(f"{OUT_DIR}/corr_sim_vs_real_brain.png", dpi=300, transparent=True)
print(f"Saved → {OUT_DIR}/corr_sim_vs_real_brain.png")
plt.close()


# ── scatter plots: MGC bins (2 Q rows × 3 bin cols) ───────────────────────────
BINS        = ["Bottom 5%", "Center 90%", "Top 5%"]
bin_palette = {"Bottom 5%": "#7B9EC8", "Center 90%": "#B2B2B2", "Top 5%": "#C8796E"}

fig, axes = plt.subplots(2, 3, figsize=(8.5, 5.5))

for row_i, q in enumerate(["Q2", "Q3"]):
    for col_i, b in enumerate(BINS):
        ax  = axes[row_i, col_i]
        grp = merged[(merged["Q"] == q) & (merged["mgc_bin"] == b)]
        r, p, n = corr_bins.get((q, b), (float("nan"), float("nan"), 0))

        ax.scatter(grp["auc_sim"], grp["auc_real"],
                   color=bin_palette[b], s=8, alpha=0.5, linewidths=0)

        if len(grp) >= 3 and not np.isnan(r):
            m, intercept = np.polyfit(grp["auc_sim"], grp["auc_real"], 1)
            x_line = np.linspace(grp["auc_sim"].min(), grp["auc_sim"].max(), 300)
            ax.plot(x_line, m * x_line + intercept, color="black", lw=1.2, zorder=5)

        stat_str = (f"r = {r:.3f}\np = {p:.2g}\nn = {n}"
                    if not np.isnan(r) else f"n = {n}\n(too few)")
        ax.text(0.05, 0.97, stat_str,
                transform=ax.transAxes, va="top", fontsize=6.5)

        ax.set_xlabel("AUC  (Simulated Bulk)", fontsize=7.5)
        ax.set_ylabel("AUC  (Real Bulk)", fontsize=7.5)
        ax.set_title(f"{q_labels[q]}\n{b}", fontsize=7.5)
        ax.tick_params(labelsize=6.5)
        sns.despine(ax=ax)

fig.suptitle("Sim Brain Labels — MGC Bins (bottom 5% / center 90% / top 5%)",
             fontsize=9, y=1.01)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/corr_sim_vs_real_brain_mgcbins.pdf", transparent=True)
plt.savefig(f"{OUT_DIR}/corr_sim_vs_real_brain_mgcbins.png", dpi=300, transparent=True)
print(f"Saved → {OUT_DIR}/corr_sim_vs_real_brain_mgcbins.png")
plt.show()
