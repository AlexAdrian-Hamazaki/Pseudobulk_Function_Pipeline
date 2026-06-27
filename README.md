# Pseudobulk Function Prediction Pipeline

**Author:** Alex Adrian-Hamazaki

---

## Overview

This repository contains the Nextflow pipelines and analysis notebooks for a thesis investigating how **cellular composition of bulk RNA-seq data affects gene function prediction** using the EGAD (Extending Gene set Analysis with Diagnostics) framework.

The core question: when bulk RNA-seq is used to build co-expression networks for gene function prediction (via neighbor voting / AUROC), how much of the predictive signal comes from the tissue's **cell type composition** rather than intrinsic gene biology?

---

## Repository Structure

```
.
├── function_prediction_pipelines/   # Nextflow pipelines
│   ├── bulkEGADPipeline/            # EGAD on real GTEx bulk RNA-seq
│   ├── ctProfilePipeline/           # EGAD on scRNA-seq cell type profiles
│   ├── EGADSimulatedGOTerms/        # EGAD on simulated bulk with controlled composition
│   └── expansionSimulation/         # EGAD while varying number of cell type profiles
└── PaperFigures/                    # Notebooks generating all paper figures
    ├── Mainfigures/
    │   ├── SimulatedBulk/           # Figures from simulated bulk analyses
    │   ├── CellTypeProfileAnalysis/ # Figures from CT profile analyses
    │   └── MGC/                     # Multi-gene co-expression figures
    └── Supplemental/                # Supplemental figure notebooks
```

---

## Methods Summary

### Gene Function Prediction with EGAD

All pipelines use **EGAD's neighbor voting** algorithm:
1. Build a Pearson co-expression network from an expression matrix (samples × genes)
2. For each GO term, predict membership using network neighborhood labels
3. Report AUROC per GO term — values > 0.5 indicate learnable signal

### Data Sources

| Source | Type | Use |
|--------|------|-----|
| GTEx | Real bulk RNA-seq (multi-tissue) | Baseline function prediction performance |
| Brain scRNA-seq | Single-cell h5ad | Cell type profiles + bulk simulation ground truth |

### Pipelines at a Glance

| Pipeline | Input | Purpose |
|----------|-------|---------|
| `bulkEGADPipeline` | GTEx h5ad | Benchmark function prediction on real bulk data |
| `ctProfilePipeline` | scRNA-seq h5ad | Function prediction from aggregated cell type profiles |
| `EGADSimulatedGOTerms` | scRNA-seq h5ad + composition JSON | Simulate bulk at controlled compositions; test real & simulated GO terms |
| `expansionSimulation` | scRNA-seq h5ad | Test how adding/removing cell type profiles changes prediction |

---

## Dependencies

All pipelines are managed with **Nextflow** and **Conda/Mamba** environments.

- **Nextflow** ≥ 23.x
- **Python** (via conda env `main_env`): `anndata`, `scanpy`, `pandas`, `numpy`
- **R**: `EGAD`, `tidyverse`

---

## Quickstart

Each pipeline has its own `nextflow.config` and `bootstrap_pipe.sh`. See the README in each subdirectory under `function_prediction_pipelines/` for pipeline-specific instructions.

```bash
# Example: run the simulated GO terms pipeline
cd function_prediction_pipelines/EGADSimulatedGOTerms
nextflow run sim_pipe.nf -entry bootstrap
```

---

## Analysis Notebooks

Paper figures are generated from Jupyter notebooks in `PaperFigures/`. See [PaperFigures/README.md](PaperFigures/README.md) for a breakdown of what each notebook produces.
