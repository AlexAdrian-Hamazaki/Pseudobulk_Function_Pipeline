# PaperFigures

Jupyter notebooks that generate all main and supplemental figures for the paper. Each notebook reads pre-computed pipeline outputs (EGAD AUROCs, cell type profiles, simulation results) and produces publication-ready figures.

---

## Structure

```
PaperFigures/
├── Mainfigures/
│   ├── SimulatedBulk/               Figures from the simulated bulk analyses
│   │   ├── sim_functions_sim_data/  Simulated GO terms + simulated bulk
│   │   ├── sim_functions_real_data/ Simulated GO terms + real GTEx bulk
│   │   ├── real_functions_real_data/ Real GO terms + real GTEx bulk
│   │   ├── real_functions_sim_data/ Real GO terms + simulated bulk
│   │   └── corr_sim_and_real_data/  Correlation: simulated vs real bulk EGAD performance
│   ├── CellTypeProfileAnalysis/     Figures from CT profile-based analyses
│   │   ├── expansion/               EGAD performance vs number of CT profiles
│   │   ├── panel_learnability/      CT profile correlation with bulk performance
│   │   └── comparison_with_compcoexpression/ CT profiles vs composite co-expression
│   └── MGC/                         Multi-gene co-expression / brain-specific gene analyses
└── Supplemental/                    Supplemental figures
    ├── GO_curation.ipynb
    └── sc_processing.ipynb
```

---

## Main Figure Notebooks

### SimulatedBulk

| Notebook / Script | Output | Description |
|-------------------|--------|-------------|
| `sim_functions_sim_data/excitatory_neurons/simulate_go_terms_for_ct.ipynb` | `exi_learnability.pdf/svg` | Learnability of excitatory neuron-enriched GO terms as a function of excitatory neuron proportion in simulated bulk |
| `sim_functions_sim_data/excitatory_neurons/analysis_exci_functions.ipynb` | `exi_function_dist.svg` | Distribution of simulated GO term gene set sizes |
| `sim_functions_sim_data/brain/simulate_go_terms_for_brain.ipynb` | `brain_performance_by_brain_mgc.svg`, `brain_ginis.svg` | Brain-enriched GO term learnability vs brain cell type composition |
| `sim_functions_real_data/exci_percent_effect.ipynb` | `profile_dominance_scatterplots/`, `top_heatmaps/` | Effect of excitatory neuron percent on real GO term AUROC in simulated bulk |
| `sim_functions_real_data/no_exci_percent_effect.ipynb` | — | Same analysis with excitatory neurons removed |
| `real_functions_real_data/AUCvsMGC.ipynb` | `AUC_MGC_Brain_bulk.pdf/svg` | EGAD AUROC vs MGC score for real GTEx brain bulk |
| `corr_sim_and_real_data/real_go/corr_sim_vs_real.py` | `corr_sim_vs_real_realgo.pdf` | Correlation of AUROC between simulated and real bulk (real GO terms) |
| `corr_sim_and_real_data/sim_brain_labels/corr_sim_vs_real.py` | `corr_sim_vs_real_brain.pdf` | Correlation using brain-simulated GO term labels |
| `corr_sim_and_real_data/sim_exi_labels/corr_sim_vs_real.py` | `corr_sim_vs_real_exi.pdf` | Correlation using excitatory neuron-simulated GO term labels |

### CellTypeProfileAnalysis

| Notebook | Output | Description |
|----------|--------|-------------|
| `expansion/Expansion.ipynb` | `expansion_model_predictions.pdf/svg`, `performance_vs_ct_breadth.svg` | EGAD AUROC as CT profile count increases; linear mixed model fit |
| `expansion/Removal.ipynb` | — | EGAD AUROC when CT profiles are removed |
| `panel_learnability/panel.ipynb` | `ct_corr_bulk.svg` | Correlation between CT profile co-expression and bulk EGAD performance |
| `comparison_with_compcoexpression/real_GOs/ct_pro_realGOs.ipynb` | `scatter_ct_vs_comp.pdf`, `scatter_ct_vs_mgc.pdf` | CT profile AUROC vs composite co-expression AUROC vs MGC (real GO terms) |
| `comparison_with_compcoexpression/brain_simGOs/ct_pro_sim_brainGOs.ipynb` | `scatter_brain_enriched.pdf`, `bootstrap_learnability_boxplot.pdf` | CT profile performance on brain-simulated GO terms; t-test results |

### MGC

| Notebook | Output | Description |
|----------|--------|-------------|
| `make_MGC_scores.ipynb` | — | Computes MGC (multi-gene co-expression) scores for brain genes |
| `MGC_enrichment.ipynb` | `brain_mgcs.pdf/svg`, `percent_exi_genes.pdf` | MGC enrichment analysis; percent excitatory neuron marker genes in top MGC gene sets |

---

## Supplemental

| Notebook | Description |
|----------|-------------|
| `GO_curation.ipynb` | Documents GO term curation decisions (size filtering, dependency removal) |
| `sc_processing.ipynb` | Documents scRNA-seq preprocessing steps (normalization, cell type annotation) |

---

## Notes

- All notebooks read from pre-computed pipeline output paths hardcoded in each notebook. Update these paths if data has moved.
- Figure outputs (`.pdf`, `.svg`) are committed alongside the notebooks in `results/` or `figs/` subdirectories.
- Most analyses use `statsmodels` for mixed linear effects models (LME) and `scipy.stats` for t-tests.
