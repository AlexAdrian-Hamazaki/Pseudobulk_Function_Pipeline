#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

# ── paths ─────────────────────────────────────────────────────────────────────
AUC_PATH = (
    "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes"
    "/bin/archive/bulkSimulationOneProfile/results"
    "/fakebulk_reallabels_attempt2/master_melted_df.csv.gz"
)
AUC_PATH = (
    "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes"
    "/bin/archive/bulkSimulationOneProfile/results"
    "/fakebulk_reallabels/master_melted_df.csv.gz"
)
# AUC_PATH = "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/EGADSimulatedGOTerms/data/real_go_feb2026_allcells/master_melted_df.csv.gz"
# AUC_PATH = "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/EGADSimulatedGOTerms/data/real_go_apr2026_noexi/master_melted_df.csv.gz"
MARKERS_PATH = "results/mast_markers_excitatory_sig.csv"
GO_ANNOT_PATH = (
    "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes"
    "/bin/preprocessing/preprocessGO_pipe/data/2024_march"
    "/data/processing/bp_annotations_withGeneData_qc_annotations.csv"
)
OUT_PDF = "results/auc_vs_percent.pdf"
OUT_PNG = "results/auc_vs_percent.png"

# ── AUC data ──────────────────────────────────────────────────────────────────
df = pd.read_csv(AUC_PATH, index_col=0)
df["exp"] = df["organism_part"].str.split("_").str[0]
df = df[df["exp"] == "exp"]

# ── marker genes ──────────────────────────────────────────────────────────────
dea = pd.read_csv(MARKERS_PATH)

# ── GO annotations ────────────────────────────────────────────────────────────
go_annot = pd.read_csv(GO_ANNOT_PATH, usecols=["GO ID", "DB_Object_Symbol", "ensembl_gene_id"])

go_sizes = go_annot.groupby("GO ID").size().reset_index(name="size")

n_in = (
    pd.merge(go_annot, dea, left_on="ensembl_gene_id", right_on="names")
    .groupby("GO ID")
    .size()
    .reset_index(name="size_in")
)

percent = pd.merge(go_sizes, n_in, on="GO ID")
percent["percent"] = percent["size_in"] / percent["size"]

# ── merge with AUC, filter ────────────────────────────────────────────────────
df_plot = df.groupby(["index", "exp"])["auc"].median().reset_index()
df_plot.columns = ["GO ID", "exp", "auc"]
df_plot = df_plot.merge(percent, on="GO ID", how="left")
# df_plot = df_plot[df_plot["size_in"] >= df_plot["size_in"].quantile(0.9)]
df_plot = df_plot[df_plot["size_in"] >= 5]

df_plot["percent"] = df_plot["percent"].fillna(0)

# ── multiple linear regression: auc ~ percent + size_in ──────────────────────
model = smf.ols("auc ~ percent + size_in", data=df_plot).fit()
print(model.summary())

p_percent = model.pvalues["percent"]
r_squared = model.rsquared

# ── scatter + regression line (percent marginalised over size_in mean) ────────
sns.set_theme(style="white", palette="muted")
sns.set_context("paper")

fig, ax = plt.subplots(figsize=(2.5, 2.5))

# draw scatter coloured by size_in
sc = ax.scatter(
    df_plot["percent"], df_plot["auc"],
    c=df_plot["size_in"], cmap="YlOrRd", s=15, alpha=0.8,
)
plt.colorbar(sc, ax=ax, label="size_in", shrink=0.8)

# marginal regression line: vary percent, fix size_in at its mean
x_range = np.linspace(df_plot["percent"].min(), df_plot["percent"].max(), 200)
pred = model.predict(pd.DataFrame({
    "percent": x_range,
    "size_in": df_plot["size_in"].mean(),
}))
ax.plot(x_range, pred, color="#D5723F", lw=2)

ax.text(
    0.05, 0.95,
    f"$R^2$ = {r_squared:.3f}\n$p_{{\\%}}$ = {p_percent:.3g}",
    transform=ax.transAxes,
    verticalalignment="top",
    fontsize=8,
)

ax.set_xlabel("% Excitatory Neuron Genes in Term", fontsize=8)
ax.set_ylabel("Learnability (AUROC)", fontsize=8)
ax.tick_params(axis="x", labelsize=7)
ax.tick_params(axis="y", labelsize=7)
sns.despine(ax=ax)

plt.tight_layout()
plt.savefig(OUT_PDF, transparent=True)
plt.savefig(OUT_PNG, dpi=300, transparent=True)
print(f"Saved → {OUT_PDF} / {OUT_PNG}")
