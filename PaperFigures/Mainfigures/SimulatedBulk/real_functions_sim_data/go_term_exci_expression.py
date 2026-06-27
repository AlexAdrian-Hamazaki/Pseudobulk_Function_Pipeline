#!/usr/bin/env python3
"""
For each GO term, identify what fraction of its genes have peak expression
in Excitatory neurons (using cell-type mean profiles), then:
  - scatter that fraction against learnability (AUROC) for both datasets
  - draw expression heatmaps for the top N GO terms by AUC
"""

import os
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from scipy.stats import pearsonr

# ── paths ─────────────────────────────────────────────────────────────────────
AUC_PATH_ALL = "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/EGADSimulatedGOTerms/data/realgo_allcells/master_melted_df.csv.gz"
AUC_PATH_NONEURONS = "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/EGADSimulatedGOTerms/data/realgo_noneur/master_melted_df.csv.gz"
CT_PROFILES_PATH = (
    "data/cell_type_profiles.csv"
)
GO_ANNOT_PATH = (
    "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes"
    "/bin/preprocessing/preprocessGO_pipe/data/2024/data/final"
    "/bp_annotations_withGeneData_qc_annotations_no_dependance.csv"
)
OUT_SCATTER_DIR = "results/profile_dominance_scatterplots"
OUT_HEATMAP_DIR = "results/top_heatmaps"
OUT_CORR_DIR    = "results/correlations"
OUT_SUMMARY_CSV = "results/go_exci_expression_summary.csv"

EXCITATORY     = "Excitatory neurons"
INHIBITORY     = "Inhibitory neurons"
NEURONS        = "Neurons"
TOP_N_HEATMAPS = 6
MIN_GENES      = 3

_dark2  = plt.cm.Dark2.colors
PALETTE = {"all cells": _dark2[4], "no neurons": _dark2[4]}
CMAP    = LinearSegmentedColormap.from_list(
    "exci_heat", ["#F7FBFF", "#6BAED6", "#08306B"]
)

os.makedirs(OUT_HEATMAP_DIR, exist_ok=True)
os.makedirs(OUT_SCATTER_DIR, exist_ok=True)
os.makedirs(OUT_CORR_DIR,    exist_ok=True)

# ── load data ─────────────────────────────────────────────────────────────────
print("Loading CT profiles...")
profiles = pd.read_csv(CT_PROFILES_PATH, index_col=0)
print(f"  {profiles.shape[0]} cell types × {profiles.shape[1]:,} genes")

print("Loading GO annotations...")
go_annot = pd.read_csv(GO_ANNOT_PATH, usecols=["GO ID", "ensembl_gene_id"])

print("Loading AUC data...")
def _load_auc(path):
    df = pd.read_csv(path, index_col=0)
    df["exp"] = df["organism_part"].str.split("_").str[0]
    df = df[df["exp"] == "exp"]
    return (df.groupby("index")["auc"].median()
              .reset_index()
              .rename(columns={"index": "GO ID"}))

def _load_auc_raw(path):
    df = pd.read_csv(path, index_col=0)
    df["exp"] = df["organism_part"].str.split("_").str[0]
    df = df[df["exp"] == "exp"]
    return df.rename(columns={"index": "GO ID"})[["GO ID", "auc"]]

auc_all       = _load_auc(AUC_PATH_ALL).rename(columns={"auc": "auc_all"})
auc_noneurons = _load_auc(AUC_PATH_NONEURONS).rename(columns={"auc": "auc_noneurons"})
auc_wide      = auc_all.merge(auc_noneurons, on="GO ID", how="inner")

# ── per-gene peak cell type (original profiles, used for heatmaps) ───────────
peak_ct = profiles.idxmax(axis=0)

# ── combined neuron profile for scatter analysis ──────────────────────────────
# Merge Excitatory + Inhibitory into a single "Neurons" row using element-wise
# max, not sum. Sum inflates the Neurons row ~2× above other cell types, causing
# idxmax to assign ~50% of all genes to Neurons regardless of GO term content.
# Max asks: "does this gene beat all other cell types in either neuron subtype?"
profiles_neurons = profiles.copy()
profiles_neurons.loc[NEURONS] = profiles.loc[[EXCITATORY, INHIBITORY]].max(axis=0)
profiles_neurons = profiles_neurons.drop(index=[EXCITATORY, INHIBITORY])
cell_types_scatter = profiles_neurons.index.tolist()
peak_ct_neurons    = profiles_neurons.idxmax(axis=0)

# ── per-GO term summary (neuron-peak fraction, used for heatmap selection) ────
print("Computing neuron-peak fraction per GO term...")
rows = []
for go_id, grp in go_annot.groupby("GO ID"):
    genes = [g for g in grp["ensembl_gene_id"].values if g in profiles_neurons.columns]
    if len(genes) < MIN_GENES:
        continue
    n_neurons = (peak_ct_neurons[genes] == NEURONS).sum()
    rows.append({
        "GO ID"           : go_id,
        "n_genes"         : len(genes),
        "n_neurons_peak"  : n_neurons,
        "pct_neurons_peak": n_neurons / len(genes),
    })

summary = pd.DataFrame(rows).merge(auc_wide, on="GO ID", how="inner")
summary["auc"] = summary["auc_all"]
summary.to_csv(OUT_SUMMARY_CSV, index=False)
print(f"  {len(summary):,} GO terms  |  saved → {OUT_SUMMARY_CSV}")

# ── per-GO fractions for every cell type (combined neuron profile) ────────────
print("Computing per-cell-type peak fractions per GO term...")
ct_rows = []
for go_id, grp in go_annot.groupby("GO ID"):
    genes = [g for g in grp["ensembl_gene_id"].values
             if g in profiles_neurons.columns]
    if len(genes) < MIN_GENES:
        continue
    peaks = peak_ct_neurons[genes]
    row   = {"GO ID": go_id, "n_genes": len(genes)}
    for ct in cell_types_scatter:
        row[f"pct_{ct}"] = (peaks == ct).sum() / len(genes)
    ct_rows.append(row)

ct_summary = pd.DataFrame(ct_rows).merge(auc_wide, on="GO ID", how="inner")

# Add Excitatory and Inhibitory fractions separately (from original peak_ct)
for _ct in [EXCITATORY, INHIBITORY]:
    _pct_col = f"pct_{_ct}"
    _fracs = []
    for _go_id, _grp in go_annot.groupby("GO ID"):
        _genes = [g for g in _grp["ensembl_gene_id"].values if g in profiles.columns]
        if len(_genes) < MIN_GENES:
            continue
        _fracs.append({"GO ID": _go_id, _pct_col: (peak_ct[_genes] == _ct).sum() / len(_genes)})
    ct_summary = ct_summary.merge(pd.DataFrame(_fracs), on="GO ID", how="left")

cell_types_scatter = cell_types_scatter + [EXCITATORY, INHIBITORY]
print(f"  {len(ct_summary):,} GO terms with per-cell-type fractions")

# ── per-cell-type scatter plots ───────────────────────────────────────────────
# Every plot shows both "all cells" AUC and "no neurons" AUC.
print("\nSaving per-cell-type scatter plots...")
sns.set_theme(style="white")
sns.set_context("paper")

for ct in cell_types_scatter:
    pct_col = f"pct_{ct}"
    sub     = ct_summary[ct_summary[pct_col] > 0].copy()
    if len(sub) < 3:
        continue

    safe_ct  = ct.lower().replace(" ", "_").replace("/", "_")
    sub2     = sub.dropna(subset=["auc_noneurons"])
    panels   = [
        ("All cells",  sub,  "auc_all",        PALETTE["all cells"]),
        ("No neurons", sub2, "auc_noneurons",   PALETTE["no neurons"]),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(3, 2), sharey=True)

    for ax, (label, data, auc_col, color) in zip(axes, panels):
        if len(data) < 3:
            ax.set_visible(False)
            continue
        x = data[pct_col].values * 100
        y = data[auc_col].values
        ax.scatter(x, y, color='steelblue', s=12, alpha=0.8, linewidths=0)
        m, b   = np.polyfit(x, y, 1)
        x_line = np.linspace(x.min(), x.max(), 200)
        ax.plot(x_line, m * x_line + b, color='black', lw=1.5)
        ax.set_title(label, fontsize=8, pad=3)
        # ax.set_ylabel("Performance (AUROC)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.xaxis.set_major_locator(plt.MultipleLocator(20))
        sns.despine(ax=ax)

    fig.text(0.5, 0, "% CT-specific genes", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.18)
    fig.savefig(f"{OUT_SCATTER_DIR}/{safe_ct}.pdf", transparent=True)
    fig.savefig(f"{OUT_SCATTER_DIR}/{safe_ct}.png", dpi=300, transparent=True)
    plt.close(fig)
    print(f"  Saved → {OUT_SCATTER_DIR}/{safe_ct}.png")

# ── per-cell-type histograms: distribution of % CT-specific genes per GO term ──
OUT_HIST_DIR = "results/cell_type_histograms"
os.makedirs(OUT_HIST_DIR, exist_ok=True)
print("\nSaving per-cell-type histograms...")

ct_hist = [ct for ct in cell_types_scatter if ct not in (EXCITATORY, INHIBITORY)]
n_ct  = len(ct_hist)
ncols = 3
nrows = int(np.ceil(n_ct / ncols))
_col_w = 180 / 25.4 / ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(180 / 25.4, nrows * _col_w * 0.85), sharex=True)
axes_flat = axes.flatten()

for ax, ct in zip(axes_flat, ct_hist):
    pct_col = f"pct_{ct}"
    vals = ct_summary[pct_col].dropna().values * 100
    ax.hist(vals, bins=20, color='steelblue', edgecolor="none")
    ax.set_title(ct, fontsize=8, pad=3)
    ax.tick_params(labelsize=7)
    ax.set_xlabel("% CT-specific genes", fontsize=7)
    ax.set_ylabel("GO terms", fontsize=7)
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    sns.despine(ax=ax)

for ax in axes_flat[n_ct:]:
    ax.set_visible(False)

fig.tight_layout()
fig.savefig(f"{OUT_HIST_DIR}/ct_pct_histograms.pdf", transparent=True)
fig.savefig(f"{OUT_HIST_DIR}/ct_pct_histograms.png", dpi=300, transparent=True)
plt.close(fig)
print(f"  Saved → {OUT_HIST_DIR}/ct_pct_histograms.png")

# ── OLS model (single model per cell type) ───────────────────────────────────
# auc ~ pct_ct * is_all_cells
# "No neurons" is baseline (is_all_cells = 0).
# Coefficients:
#   pct_ct            → slope in No neurons
#   is_all_cells      → intercept shift for All cells
#   pct_ct:is_all_cells → change in slope for All cells (the interaction)
# Marginal All-cells slope = pct_ct + pct_ct:is_all_cells (tested via Wald contrast)
from statsmodels.formula.api import mixedlm
import matplotlib as mpl
import matplotlib.ticker as ticker
from plotnine import (
    ggplot, aes, geom_tile,
    scale_fill_cmap, scale_fill_gradient2, scale_color_manual,
    labs, theme_bw, theme, element_text, element_blank,
    guides, guide_colorbar,
)

_slope_cmap = mpl.colors.LinearSegmentedColormap.from_list(
    "slopes_div", ["#2471A3", "white", "#C0392B"]
)
mpl.colormaps.register(_slope_cmap, force=True)

print("\nLoading raw AUC data for mixed effects models...")
auc_all_raw   = _load_auc_raw(AUC_PATH_ALL)
auc_noexi_raw = _load_auc_raw(AUC_PATH_NONEURONS)

df_all        = auc_all_raw.copy();   df_all["is_all_cells"]   = 1
df_no         = auc_noexi_raw.copy(); df_no["is_all_cells"]    = 0
auc_combined  = pd.concat([df_all, df_no], ignore_index=True)

print("Fitting mixed effects models (this may take a minute)...")
model_rows = []
slope_rows = []
for ct in cell_types_scatter:
    if ct in (EXCITATORY, INHIBITORY):
        continue
    pct_col  = f"pct_{ct}"
    pct_data = ct_summary[["GO ID", pct_col]].dropna()

    merged = auc_combined.merge(pct_data, on="GO ID", how="inner")
    merged = merged[merged[pct_col] > 0].copy()
    merged = merged.rename(columns={pct_col: "pct_ct"})
    merged["pct_ct"] = merged["pct_ct"] * 10   # fraction → units of 10%

    if merged["GO ID"].nunique() < 5:
        continue

    try:
        # Random intercept per GO term accounts for repeated measures across
        # the two conditions (all-cells vs no-neurons) for the same GO term.
        result = mixedlm(
            "auc ~ pct_ct * is_all_cells", merged,
            groups=merged["GO ID"],
        ).fit(reml=False)

        # fixed-effect coefficients
        for term in result.fe_params.index:
            model_rows.append({
                "cell_type": ct,
                "term"     : term,
                "coef"     : result.fe_params[term],
                "pval"     : result.pvalues[term],
                "n_terms"  : merged["GO ID"].nunique(),
                "n_obs"    : len(merged),
            })

        # random intercept variance (per-GO-term)
        model_rows.append({
            "cell_type": ct,
            "term"     : "RE_var_intercept",
            "coef"     : float(result.cov_re.values[0, 0]),
            "pval"     : np.nan,
            "n_terms"  : merged["GO ID"].nunique(),
            "n_obs"    : len(merged),
        })

        # No neurons marginal slope = pct_ct fixed effect
        slope_no = result.fe_params["pct_ct"]
        pval_no  = result.pvalues["pct_ct"]

        # All cells marginal slope = pct_ct + pct_ct:is_all_cells
        # p-value is the interaction term: tests whether slope differs from No Neurons
        slope_all = float(result.fe_params["pct_ct"] + result.fe_params["pct_ct:is_all_cells"])
        pval_all  = result.pvalues["pct_ct:is_all_cells"]

        for cond, slope, pval in [
            ("No Neurons", slope_no,  pval_no),
            ("All cell types",  slope_all, pval_all),
        ]:
            slope_rows.append({
                "cell_type": ct,
                "condition": cond,
                "slope"    : slope,
                "pval"     : pval,
            })

        print(f"  {ct}: slope_noneur={slope_no:.3f}(p={pval_no:.3g})  slope_all={slope_all:.3f}(p={pval_all:.3g})")
    except Exception as e:
        print(f"  FAILED {ct}: {e}")

model_df = pd.DataFrame(model_rows)
model_df.to_csv(f"{OUT_SCATTER_DIR}/lme_results.csv", index=False)
print(f"Saved → {OUT_SCATTER_DIR}/lme_results.csv")

# ── effect-size heatmap ───────────────────────────────────────────────────────
_fe_terms = ["Intercept", "pct_ct", "is_all_cells", "pct_ct:is_all_cells"]
plot_df = model_df[model_df["term"].isin(_fe_terms)].copy()
plot_df["sig"] = (plot_df["pval"] < 0.05).map({True: "p<0.05", False: "ns"})
plot_df["term"] = pd.Categorical(plot_df["term"], categories=_fe_terms)

p_lme = (
    ggplot(plot_df, aes(x="cell_type", y="term", fill="coef", color="sig"))
    + geom_tile(size=1.2)
    + scale_fill_gradient2(
        low="#2471A3", mid="white", high="#C0392B",
        midpoint=0,
        name="ΔAUC per 10 pp increase in\ncell-type gene dominance",
    )
    + scale_color_manual(values={"p<0.05": "#4D4D4D", "ns": "white"})
    + guides(
        fill=guide_colorbar(title_position="top", title_hjust=0.5),
        color=False,
    )
    + labs(x="", y="")
    + theme_bw()
    + theme(
        figure_size=(0.45 * len(cell_types_scatter) + 1, 3),
        axis_text=element_text(size=8),
        axis_text_x=element_text(angle=30, ha="right"),
        panel_grid=element_blank(),
    )
)
p_lme.save(f"{OUT_SCATTER_DIR}/lme_effect_heatmap.pdf", verbose=False)
p_lme.save(f"{OUT_SCATTER_DIR}/lme_effect_heatmap.png", dpi=300, verbose=False)
print(f"LME heatmap saved → {OUT_SCATTER_DIR}/lme_effect_heatmap.png")

# ── slope heatmap: AUC gain per unit increase in cell-type gene fraction ───────
slope_df = pd.DataFrame(slope_rows)
slope_df.to_csv(f"{OUT_SCATTER_DIR}/slopes.csv", index=False)

slope_df["sig"] = (slope_df["pval"] < 0.05).map({True: "p<0.05", False: "ns"})
slope_df["condition"] = pd.Categorical(
    slope_df["condition"], categories=["No Neurons", "All cell types"]
)

_slope_abs_max = max(abs(slope_df["slope"].min()), abs(slope_df["slope"].max()))
_slope_mag      = 10 ** np.floor(np.log10(_slope_abs_max))
_slope_nice_max = float(np.ceil(_slope_abs_max / _slope_mag) * _slope_mag)
_slope_step     = float(_slope_mag)
_slope_decimals = max(0, int(-np.floor(np.log10(_slope_step))))
_slope_ticks    = np.arange(-_slope_nice_max, _slope_nice_max + _slope_step / 2, _slope_step)

p_slopes = (
    ggplot(slope_df, aes(x="cell_type", y="condition", fill="slope", color="sig"))
    + geom_tile(size=1.2)
    + scale_fill_cmap("slopes_div", limits=(-_slope_nice_max, _slope_nice_max))
    + scale_color_manual(values={"p<0.05": "#4D4D4D", "ns": "white"})
    + guides(color=False)
    + labs(x="", y="")
    + theme_bw()
    + theme(
        figure_size=(2, 1),
        axis_text=element_text(size=8),
        axis_text_x=element_text(angle=30, ha="right"),
        panel_grid=element_blank(),
        legend_position="none",
    )
)
p_slopes.save(f"{OUT_SCATTER_DIR}/slope_heatmap.pdf", verbose=False)
p_slopes.save(f"{OUT_SCATTER_DIR}/slope_heatmap.png", dpi=300, verbose=False)
print(f"Slope heatmap saved → {OUT_SCATTER_DIR}/slope_heatmap.png")

# ── slope heatmap legend (separate figure) ────────────────────────────────────
_cb_norm = mpl.colors.Normalize(vmin=-_slope_nice_max, vmax=_slope_nice_max)

fig_cb = plt.figure(figsize=(1.2, 1.0))
ax_cb = fig_cb.add_axes([0.05, 0.08, 0.09, 0.84])
sm = plt.cm.ScalarMappable(cmap=_slope_cmap, norm=_cb_norm)
sm.set_array([])
cb = fig_cb.colorbar(sm, cax=ax_cb, orientation="vertical")
cb.set_label(
    "ΔAUC per 10% increase in cell-type gene dominance",
    fontsize=7, rotation=90, labelpad=6, va="center",
)
cb.set_ticks(_slope_ticks)
cb.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter(f"%.{_slope_decimals}f"))
cb.ax.tick_params(labelsize=6)
fig_cb.savefig(f"{OUT_SCATTER_DIR}/slope_heatmap_legend.pdf", transparent=True, bbox_inches="tight")
fig_cb.savefig(f"{OUT_SCATTER_DIR}/slope_heatmap_legend.png", dpi=300, transparent=True, bbox_inches="tight")
plt.close(fig_cb)
print(f"Slope legend saved → {OUT_SCATTER_DIR}/slope_heatmap_legend.png")

# ── heatmap helpers ───────────────────────────────────────────────────────────
profiles_glial = profiles.drop(index=[EXCITATORY, INHIBITORY])

def _norm(prof, gene_order):
    e = prof[gene_order]
    return e.div(e.max(axis=0).replace(0, np.nan)).fillna(0)

def _genes_for_go(go_id):
    genes = [g for g in
             go_annot.loc[go_annot["GO ID"] == go_id, "ensembl_gene_id"].values
             if g in profiles.columns]
    gene_peak = peak_ct[genes]
    is_exci   = gene_peak == EXCITATORY
    order     = list(gene_peak[is_exci].index) + list(gene_peak[~is_exci].index)
    n_exci    = int(is_exci.sum())
    return order, n_exci

def _draw_heatmap(ax, expr_norm, n_exci, title):
    data = expr_norm.values
    im = ax.imshow(data, aspect="auto", cmap=CMAP, vmin=0, vmax=1,
                   interpolation="nearest")
    if 0 < n_exci < data.shape[1]:
        ax.axvline(n_exci - 0.5, color="white", lw=1.2, ls="--")
    ax.set_yticks(range(len(expr_norm.index)))
    ax.set_yticklabels(expr_norm.index, fontsize=6)
    ax.set_xticks([])
    ax.set_xlabel("")
    ax.set_title(title, fontsize=7, pad=3)
    return im

# ── top GO terms ──────────────────────────────────────────────────────────────
top_terms = summary.nlargest(TOP_N_HEATMAPS, "auc")["GO ID"].tolist()

# ── Pearson r between excitatory and inhibitory profiles per top GO term ──────
print("\nComputing Exci vs Inh Pearson correlations for top GO terms...")
corr_rows = []
for go_id in top_terms:
    genes = [g for g in
             go_annot.loc[go_annot["GO ID"] == go_id, "ensembl_gene_id"].values
             if g in profiles.columns]
    exci_vec = profiles.loc[EXCITATORY, genes].values.astype(float)
    inh_vec  = profiles.loc[INHIBITORY,  genes].values.astype(float)
    r, p = pearsonr(exci_vec, inh_vec)
    auc_all_val   = summary.loc[summary["GO ID"] == go_id, "auc_all"].values[0]
    auc_noneuro   = summary.loc[summary["GO ID"] == go_id, "auc_noneurons"].values[0]
    n_neurons_peak   = summary.loc[summary["GO ID"] == go_id, "n_neurons_peak"].values[0]
    corr_rows.append({
        "GO ID"              : go_id,
        "n_genes"            : len(genes),
        "n_neurons_peak"        : n_neurons_peak,
        "pearson_r_exci_inh" : r,
        "pearson_p_exci_inh" : p,
        "auc_all"            : auc_all_val,
        "auc_noneurons"      : auc_noneuro,
    })
    print(f"  {go_id}  r = {r:.3f}  p = {p:.3g}")

corr_df = pd.DataFrame(corr_rows)
out_corr_csv = f"{OUT_HEATMAP_DIR}/top_go_exci_inh_correlation.csv"
corr_df.to_csv(out_corr_csv, index=False)
print(f"Saved → {out_corr_csv}")

# ── pairwise cell-type Pearson correlation matrices per top GO term ───────────
print(f"\nSaving cell-type pairwise Pearson correlation matrices → {OUT_CORR_DIR}/")
for go_id in top_terms:
    genes = [g for g in
             go_annot.loc[go_annot["GO ID"] == go_id, "ensembl_gene_id"].values
             if g in profiles.columns]
    ct_corr = profiles[genes].T.corr()
    safe_id = go_id.replace(":", "_")
    out_path = f"{OUT_CORR_DIR}/{safe_id}_cell_type_pearson.csv"
    ct_corr.to_csv(out_path)
    print(f"  {go_id} ({len(genes)} genes) → {out_path}")

# ── individual heatmaps: one file per (GO term × condition) ──────────────────
print(f"\nSaving individual heatmaps for top {TOP_N_HEATMAPS} GO terms...")
for go_id in top_terms:
    order, n_exci = _genes_for_go(go_id)
    auc_all_val   = summary.loc[summary["GO ID"] == go_id, "auc_all"].values[0]
    auc_noneuro   = summary.loc[summary["GO ID"] == go_id, "auc_noneurons"].values[0]
    safe_id       = go_id.replace(":", "_")

    for cond_label, prof, auc_val in [
        ("all_cells",  profiles,       auc_all_val),
        ("no_neurons", profiles_glial, auc_noneuro),
    ]:
        expr_norm = _norm(prof, order)
        n_ct      = len(prof.index)
        fig, ax   = plt.subplots(figsize=(5, max(1.8, n_ct * 0.24 + 0.4)))
        im = _draw_heatmap(
            ax, expr_norm, n_exci,
            f"{go_id}  AUC = {auc_val:.3f}",
        )
        cbar = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
        cbar.set_label("Scaled Expression", fontsize=7)
        cbar.ax.tick_params(labelsize=6)
        fig.tight_layout()
        fig.savefig(
            f"{OUT_HEATMAP_DIR}/{safe_id}_{cond_label}.pdf",
            transparent=True, bbox_inches="tight",
        )
        fig.savefig(
            f"{OUT_HEATMAP_DIR}/{safe_id}_{cond_label}.png",
            dpi=200, transparent=True, bbox_inches="tight",
        )
        plt.close(fig)
    print(f"  {go_id}")

# ── summary heatmaps: all GO terms stacked, one per condition ─────────────────
for cond_label, prof, auc_col, out_pdf, out_png in [
    (
        "All cells", profiles, "auc_all",
        f"{OUT_HEATMAP_DIR}/summary_all_cells.pdf",
        f"{OUT_HEATMAP_DIR}/summary_all_cells.png",
    ),
    (
        "No neurons", profiles_glial, "auc_noneurons",
        f"{OUT_HEATMAP_DIR}/summary_no_neurons.pdf",
        f"{OUT_HEATMAP_DIR}/summary_no_neurons.png",
    ),
]:
    n_ct  = len(prof.index)
    fig, axes = plt.subplots(
        TOP_N_HEATMAPS, 1,
        figsize=(6, TOP_N_HEATMAPS * (n_ct * 0.24 + 0.7)),
    )
    if TOP_N_HEATMAPS == 1:
        axes = [axes]

    last_im = None
    for ax, go_id in zip(axes, top_terms):
        order, n_exci = _genes_for_go(go_id)
        auc_val       = summary.loc[summary["GO ID"] == go_id, auc_col].values[0]
        expr_norm     = _norm(prof, order)
        last_im = _draw_heatmap(
            ax, expr_norm, n_exci,
            f"{go_id}  AUC = {auc_val:.3f}",
        )

    cbar = fig.colorbar(last_im, ax=list(axes), shrink=0.4, pad=0.02)
    cbar.set_label("Scaled Expression", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    fig.suptitle(f"Top {TOP_N_HEATMAPS} GO terms — {cond_label}", fontsize=10, y=1.01)
    fig.tight_layout()
    fig.savefig(out_pdf, transparent=True, bbox_inches="tight")
    fig.savefig(out_png, dpi=150, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"Summary heatmap saved → {out_pdf}")
