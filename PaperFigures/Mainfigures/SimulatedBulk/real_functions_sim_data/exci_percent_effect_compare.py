#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

# ── paths ─────────────────────────────────────────────────────────────────────
DATASETS = {
    "all":   "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/EGADSimulatedGOTerms/data/real_go_feb2026_allcells/master_melted_df.csv.gz",
    "noexi": "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/EGADSimulatedGOTerms/data/realgo_noneur/master_melted_df.csv.gz"
}
# "noexi": "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/EGADSimulatedGOTerms/data/real_go_apr2026_noexi/master_melted_df.csv.gz",
MARKERS_PATH = "results/mast_markers_excitatory_sig.csv"
GO_ANNOT_PATH = (
    "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes"
    "/bin/preprocessing/preprocessGO_pipe/data/2024_march"
    "/data/processing/bp_annotations_withGeneData_qc_annotations.csv"
)
OUT_PDF  = "results/auc_vs_percent_compare.pdf"
OUT_PNG  = "results/auc_vs_percent_compare.png"
OUT_HIST_PDF  = "results/auc_diff_histogram.pdf"
OUT_HIST_PNG  = "results/auc_diff_histogram.png"
OUT_DROP_PDF  = "results/auc_high_vs_drop.pdf"
OUT_DROP_PNG  = "results/auc_high_vs_drop.png"

# ── marker genes & GO percent ─────────────────────────────────────────────────
dea = pd.read_csv(MARKERS_PATH)
go_annot = pd.read_csv(GO_ANNOT_PATH, usecols=["GO ID", "DB_Object_Symbol", "ensembl_gene_id"])

go_sizes = go_annot.groupby("GO ID").size().reset_index(name="size")
n_in = (
    pd.merge(go_annot, dea, left_on="ensembl_gene_id", right_on="names")
    .groupby("GO ID").size().reset_index(name="size_in")
)
percent = pd.merge(go_sizes, n_in, on="GO ID")
percent["percent"] = percent["size_in"] / percent["size"]

# ── load and stack both AUC datasets ─────────────────────────────────────────
frames = []
for label, path in DATASETS.items():
    df = pd.read_csv(path, index_col=0)
    df["exp"] = df["organism_part"].str.split("_").str[0]
    df = df[df["exp"] == "exp"]
    agg = df.groupby("index")["auc"].median().reset_index()
    agg.columns = ["GO ID", "auc"]
    agg["dataset"] = label
    frames.append(agg)

df_plot = pd.concat(frames, ignore_index=True)
df_plot = df_plot.merge(percent, on="GO ID", how="left")
df_plot = df_plot[df_plot["size_in"] >= 5]
df_plot["percent"] = df_plot["percent"].fillna(0)

# ── model: auc ~ percent + size_in + dataset ─────────────────────────────────
# "all" is the reference level; coefficient on dataset[T.noexi] captures the
# mean AUC shift between the two conditions after controlling for percent and size
model = smf.ols("auc ~ percent * C(dataset, Treatment('all')) + size_in", data=df_plot).fit()
print(model.summary())

r_squared    = model.rsquared
p_percent    = model.pvalues["percent"]
p_dataset    = model.pvalues["C(dataset, Treatment('all'))[T.noexi]"]
p_interaction = model.pvalues["percent:C(dataset, Treatment('all'))[T.noexi]"]

# ── plot ──────────────────────────────────────────────────────────────────────
_dark2   = plt.cm.Dark2.colors
palette  = {"all": _dark2[0], "noexi": _dark2[1]}
DISPLAY  = {"all": "With neurons", "noexi": "No excitatory neurons"}

sns.set_theme(style="white")
sns.set_context("paper")

fig, ax = plt.subplots(figsize=(3, 2.5))

for label, grp in df_plot.groupby("dataset"):
    ax.scatter(grp["percent"], grp["auc"],
               color=palette[label], s=10, alpha=0.6, label=DISPLAY[label])

# marginal regression lines per dataset (size_in fixed at overall mean)
x_range   = np.linspace(df_plot["percent"].min(), df_plot["percent"].max(), 200)
mean_size = df_plot["size_in"].mean()
for label in DATASETS:
    pred = model.predict(pd.DataFrame({
        "percent": x_range,
        "size_in": mean_size,
        "dataset": label,
    }))
    ax.plot(x_range, pred, color=palette[label], lw=2)

ax.text(
    0.05, 0.95,
    f"$R^2$ = {r_squared:.3f}\n$p_{{\\%}}$ = {p_percent:.3g}\n$p_{{dataset}}$ = {p_dataset:.3g}\n$p_{{interaction}}$ = {p_interaction:.3g}",
    transform=ax.transAxes,
    verticalalignment="top",
    fontsize=7,
)

ax.set_xlabel("% Excitatory Neuron Genes in Term", fontsize=8)
ax.set_ylabel("Learnability (AUROC)", fontsize=8)
ax.tick_params(labelsize=7)
ax.legend(title="Condition", fontsize=7, title_fontsize=7, frameon=False)
sns.despine(ax=ax)

plt.tight_layout()
plt.savefig(OUT_PDF, transparent=True)
plt.savefig(OUT_PNG, dpi=300, transparent=True)
print(f"Saved → {OUT_PDF} / {OUT_PNG}")

# ── histogram of per-GO AUC differences (all − noexi) ────────────────────────
wide = df_plot.pivot(index="GO ID", columns="dataset", values="auc").dropna()
wide["diff"] = wide["all"] - wide["noexi"]

from scipy.stats import ttest_1samp
t_stat, p_val = ttest_1samp(wide["diff"], popmean=0)

print("\nTop 10 highest-performing GO terms in all-cells dataset:")
print(wide.sort_values("all", ascending=False).head(10)[["all", "noexi", "diff"]].to_string())
median_diff = wide["diff"].median()

fig2, ax2 = plt.subplots(figsize=(3, 2.5))
ax2.hist(wide["diff"], bins=40, color="#5A7D9A", edgecolor="white", linewidth=0.3)
ax2.axvline(0,           color="black", lw=1,   ls="--")
ax2.axvline(median_diff, color="#D5723F", lw=1.5, ls="-", label=f"median = {median_diff:.3f}")

ax2.text(
    0.97, 0.95,
    f"t = {t_stat:.2f}\np = {p_val:.3g}\nn = {len(wide)}",
    transform=ax2.transAxes,
    ha="right", va="top",
    fontsize=7,
)

ax2.set_xlabel("AUC (all cells) − AUC (no excitatory)", fontsize=8)
ax2.set_ylabel("GO terms", fontsize=8)
ax2.tick_params(labelsize=7)
ax2.legend(fontsize=7)
sns.despine(ax=ax2)

fig2.tight_layout()
fig2.savefig(OUT_HIST_PDF, transparent=True)
fig2.savefig(OUT_HIST_PNG, dpi=300, transparent=True)
print(f"Saved → {OUT_HIST_PDF} / {OUT_HIST_PNG}")
print(f"Median diff = {median_diff:.4f},  t = {t_stat:.3f},  p = {p_val:.3g}")

# ── do high-AUC terms (all cells) drop the most when excitatory are removed? ──
from scipy.stats import pearsonr, spearmanr

# merge percent back so we can colour points by excitatory marker content
wide_pct = wide.merge(percent[["GO ID", "percent"]], on="GO ID", how="left").fillna(0)

r_pearson, p_pearson = pearsonr(wide_pct["all"], wide_pct["diff"])
r_spearman, p_spearman = spearmanr(wide_pct["all"], wide_pct["diff"])
print(f"Pearson  r = {r_pearson:.3f},  p = {p_pearson:.3g}")
print(f"Spearman r = {r_spearman:.3f},  p = {p_spearman:.3g}")

fig3, ax3 = plt.subplots(figsize=(3, 2.5))

sc3 = ax3.scatter(
    wide_pct["all"], wide_pct["diff"],
    c=wide_pct["percent"], cmap="YlOrRd",
    s=12, alpha=0.7, linewidths=0,
)
plt.colorbar(sc3, ax=ax3, label="% excitatory markers", shrink=0.85)

# regression line
m, b = np.polyfit(wide_pct["all"], wide_pct["diff"], 1)
x_line = np.linspace(wide_pct["all"].min(), wide_pct["all"].max(), 200)
ax3.plot(x_line, m * x_line + b, color="#D5723F", lw=1.5)

ax3.axhline(0, color="black", lw=0.8, ls="--")

ax3.text(
    0.05, 0.95,
    f"Pearson $r$ = {r_pearson:.3f}\n$p$ = {p_pearson:.3g}\n"
    f"Spearman $\\rho$ = {r_spearman:.3f}",
    transform=ax3.transAxes,
    va="top", fontsize=7,
)

ax3.set_xlabel("AUC (all cells)", fontsize=8)
ax3.set_ylabel("AUC drop (all − no excitatory)", fontsize=8)
ax3.tick_params(labelsize=7)
sns.despine(ax=ax3)

fig3.tight_layout()
fig3.savefig(OUT_DROP_PDF, transparent=True)
fig3.savefig(OUT_DROP_PNG, dpi=300, transparent=True)
print(f"Saved → {OUT_DROP_PDF} / {OUT_DROP_PNG}")
