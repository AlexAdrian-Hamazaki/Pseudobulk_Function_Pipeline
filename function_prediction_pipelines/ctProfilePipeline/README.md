# ctProfilePipeline

Runs EGAD gene function prediction directly on **cell type profiles** derived from scRNA-seq data. Each cell type profile is the average gene expression across all cells of that type, producing a compact representation of each cell type's transcriptional signature.

---

## What It Does

1. **Computes cell type profiles** — reads a single-cell h5ad and averages gene expression per cell type, producing one row per cell type
2. **Runs EGAD** — builds a Pearson co-expression network from the cell type profile matrix and predicts GO term membership via neighbor voting (AUROC)

The key insight this pipeline tests: can a matrix with **only ~8–15 rows** (one per cell type) still produce meaningful co-expression signal for function prediction?

---

## Pipeline Steps

```
scRNA-seq h5ad
    └─ calc_ct_profiles    Average expression per cell type → brain_profiles.csv.gz
         └─ EGAD           Build co-expression network + EGAD (Egad.R)
```

---

## Key Scripts

| Script | Language | Description |
|--------|----------|-------------|
| `main.nf` | Nextflow | Main pipeline definition |
| `bin/makeCellTypeProfiles.py` | Python | Groups cells by type, computes mean expression per cell type |
| `bin/Egad.R` | R | Builds Pearson co-expression network, runs EGAD neighbor voting, outputs AUROC CSV |

---

## Configuration (`nextflow.config`)

| Parameter | Description |
|-----------|-------------|
| `params.path_to_sc_adata` | Path to scRNA-seq h5ad file (cells × genes, with `obs["Cell type"]` annotation) |
| `params.go_annotations` | Path to GO BP annotation CSV |
| `params.gene_column` | Gene name column in GO file |
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
├── brain_profiles.csv.gz    Cell type profile matrix (cell types × genes)
└── EGAD/
    └── *_EGAD.csv           EGAD AUROC results per GO term
```

---

## Cell Type Profile Construction

`makeCellTypeProfiles.py` groups cells by the `obs["Cell type"]` column and computes the mean expression vector for each group. NaN values are filled with 0 before averaging. The resulting matrix has shape `(n_cell_types, n_genes)`.

This matrix is then treated as if each cell type were a "sample" — so the co-expression network captures how genes co-vary **across cell types** rather than across individual cells or bulk samples.
