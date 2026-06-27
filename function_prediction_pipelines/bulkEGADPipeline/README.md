# bulkEGADPipeline

Runs EGAD gene function prediction on **real bulk RNA-seq data** from GTEx. This pipeline establishes the baseline performance of function prediction before any compositional manipulation.

---

## What It Does

1. **Loads GTEx** — reads a merged multi-tissue h5ad expression file
2. **Splits by tissue** — extracts each organism part (e.g., Brain, Blood) into its own CSV
3. **Bootstrap subsampling** — for each tissue, repeatedly samples N bulk samples with replacement
4. **Runs EGAD** — builds a Pearson co-expression network per bootstrap and predicts GO term membership via neighbor voting (AUROC)
5. **Merges results** — collects all bootstrap EGAD outputs into one melted dataframe per tissue

---

## Pipeline Steps

```
GTEx h5ad
    └─ splitMerged            Split into per-tissue CSVs
         └─ sampleWithReplacement    Bootstrap subsample N samples
              └─ sim_bulk_EGAD       Build co-expression network + EGAD (bulkEgad.R)
                   └─ makeMeltedMergedDF   Collect all bootstrap EGADs → melted CSV
```

---

## Key Scripts

| Script | Language | Description |
|--------|----------|-------------|
| `bulkPipe.nf` | Nextflow | Main pipeline definition |
| `bin/bulkEgad.R` | R | Builds Pearson co-expression network, runs EGAD neighbor voting, outputs AUROC CSV |
| `bin/splitadata.py` | Python | Splits merged h5ad by organism part into CSV |
| `bin/subsampleBulk.py` | Python | Subsamples N rows (bulk samples) with replacement |
| `bin/adata_to_csv.py` | Python | Converts h5ad to CSV for EGAD input |
| `bin/make_merged_EGAD_of_splits.py` | Python | Merges per-bootstrap EGAD CSVs into one melted dataframe |
| `bootstrap_pipe.sh` | Bash | Wrapper to run multiple bootstrap iterations in parallel |

---

## Configuration (`nextflow.config`)

| Parameter | Description |
|-----------|-------------|
| `params.bulk_merged` | Path to GTEx merged h5ad file |
| `params.organism_parts` | List of tissues to process (e.g., `["Brain", "Blood"]`) |
| `params.go_annotations` | Path to GO BP annotation CSV |
| `params.gene_column` | Gene name column in GO file (`DB_Object_Symbol`) |
| `params.num_bootstrap` | Number of bootstrap iterations |
| `params.bulk_size` | Number of bulk samples per bootstrap subsample |
| `params.publish` | Output directory |

---

## Running

### Bootstrap mode (recommended)

```bash
nextflow run bulkPipe.nf -entry bootstrap
```

This runs `num_bootstrap` iterations, each subsampling `bulk_size` samples per tissue.

### Bootstrap shell wrapper (parallel runs)

```bash
./bootstrap_pipe.sh <num_bootstraps> <publishDir> [max_parallel]
# e.g.
./bootstrap_pipe.sh 100 data/bootstrapped 4
```

### Run on full GTEx (all tissues, no split)

```bash
nextflow run bulkPipe.nf -entry run_on_all_Gtex
```

---

## Outputs

```
results/
└── <publishDir>/
    ├── splits/
    │   ├── <tissue>_split.csv.gz            Per-tissue expression matrices
    │   └── subsamples/
    │       └── <tissue>_<bootstrap>.csv.gz  Bootstrapped subsamples
    └── EGAD/
        ├── *_EGAD.csv                       Per-bootstrap EGAD AUROCs
        └── melted_dfs/
            └── <tissue>_melted_EGADs.csv.gz Merged results for plotting
```

---

## EGAD Method Notes

`bulkEgad.R` does the following:
- Loads the expression CSV (samples × genes)
- Computes Pearson correlations across all gene pairs → co-expression network
- Filters GO terms to those with ≥ 20 measured genes
- Runs `neighbor_voting()` with 3-fold cross-validation
- Outputs AUROC per GO term
