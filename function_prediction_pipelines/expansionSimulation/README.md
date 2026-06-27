# expansionSimulation

Tests how the **number of cell type profiles** included in the expression matrix affects EGAD gene function prediction performance. Starting from all available cell type profiles, the pipeline iteratively subsamples smaller sets and measures AUROC at each step.

This directly addresses the question: does adding more cell type diversity improve or degrade function prediction?

---

## What It Does

1. **Compute cell type profiles** — average expression per cell type from scRNA-seq h5ad
2. **Subsample iteratively** — create a series of matrices, each containing a different number of cell type profiles (stepping by `ct_profile_step`)
3. **Run EGAD** — build a co-expression network from each subsampled matrix and predict GO term membership

---

## Pipeline Steps

```
scRNA-seq h5ad
    └─ calc_ct_profiles       Average expression per cell type → all_human_CT_profiles.csv
         └─ subsample_ct_profiles   Generate subsampled CT profile matrices at each step
              └─ EGAD              Build co-expression network + EGAD (Egad.R) for each matrix
```

---

## Key Scripts

| Script | Language | Description |
|--------|----------|-------------|
| `main.nf` | Nextflow | Main pipeline definition |
| `bin/makeCellTypeProfiles.py` | Python | Computes mean expression per cell type from scRNA-seq h5ad |
| `bin/downsampleCTProfiles.py` | Python | Iteratively samples subsets of CT profiles, stepping by `ct_profile_step` |
| `bin/Egad.R` | R | Builds Pearson co-expression network and runs EGAD neighbor voting |

---

## Subsampling Strategy (`downsampleCTProfiles.py`)

Starting from the full set of cell type profiles:
- Begin with `ct_profile_step` profiles sampled randomly
- Increase by `ct_profile_step` each iteration until all profiles are included
- Each subset is saved as `{n}_ct_pros.csv`

This produces a ladder of datasets from small (few cell types) to full (all cell types), allowing EGAD performance to be plotted as a function of cell type breadth.

---

## Configuration (`nextflow.config`)

| Parameter | Description |
|-----------|-------------|
| `params.path_to_sc_adata` | Path to scRNA-seq h5ad file |
| `params.go_annotations` | Path to GO BP annotation CSV |
| `params.min_ct_profile` | Minimum number of CT profiles to start from |
| `params.ct_profile_step` | Step size for adding CT profiles each iteration |
| `params.publish` | Output directory |
| `params.python3_9` | Path to conda environment |

---

## Running

```bash
nextflow run main.nf
```

---

## Outputs

```
<publish>/
├── all_human_CT_profiles.csv    Full cell type profile matrix
├── ct_profiles/
│   └── {n}_ct_pros.csv          Subsampled CT profile matrices (one per step)
└── EGAD/
    └── *_EGAD.csv               EGAD AUROC results per subsampled matrix
```

The EGAD outputs can then be joined on the number of CT profiles included to plot AUROC vs. cell type breadth.
