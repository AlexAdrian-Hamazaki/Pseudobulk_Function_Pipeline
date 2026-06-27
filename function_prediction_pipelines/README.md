# Function Prediction Pipelines

This directory contains four Nextflow pipelines that together form the gene function prediction analysis. Each pipeline is independent and addresses a distinct experimental question.

---

## Pipelines

| Directory | Entry Point | Purpose |
|-----------|-------------|---------|
| [`bulkEGADPipeline/`](bulkEGADPipeline/) | `bulkPipe.nf -entry bootstrap` | Run EGAD on real GTEx bulk RNA-seq, bootstrapped |
| [`ctProfilePipeline/`](ctProfilePipeline/) | `main.nf` | Run EGAD directly on scRNA-seq cell type profiles |
| [`EGADSimulatedGOTerms/`](EGADSimulatedGOTerms/) | `sim_pipe.nf -entry bootstrap` | Simulate pseudo-bulk from cell type profiles at controlled compositions; run EGAD |
| [`expansionSimulation/`](expansionSimulation/) | `main.nf` | Systematically add/remove cell type profiles and measure EGAD performance |

---

## Shared Design Patterns

All pipelines follow the same general pattern:

1. **Expression matrix in** — either real bulk (GTEx) or simulated from scRNA-seq cell type profiles
2. **EGAD neighbor voting** — builds a Pearson co-expression network and predicts GO term membership via AUROC
3. **GO annotation filter** — only GO Biological Process terms with ≥ 20 measured genes are used
4. **Results out** — EGAD AUROC CSVs, melted for downstream plotting

### Shared Tools

| Script | Language | Role |
|--------|----------|------|
| `Egad.R` / `bulkEgad.R` / `simBulkEgad.R` | R | Build co-expression network and run EGAD neighbor voting |
| `makeCellTypeProfiles.py` | Python | Average gene expression per cell type from an AnnData h5ad |
| `simBulk.py` | Python | Simulate bulk samples by scaling cell type profiles with stochastic compositions |

### Common Parameters

| Parameter | Description |
|-----------|-------------|
| `params.go_annotations` | Path to GO Biological Process annotation CSV |
| `params.gene_column` | Gene identifier column (`DB_Object_Symbol` or `ensembl_gene_id`) |
| `params.publish` | Output directory for results |
| `params.python3_9` | Path to conda environment for Python scripts |

---

## Running Pipelines

All pipelines assume:
- **Nextflow** is installed and on `PATH`
- **Mamba/Conda** is available (set `conda.enabled = true` and `conda.useMamba = true` in config)
- Input data paths in `nextflow.config` are updated to your local paths

```bash
# General pattern
cd <pipeline_dir>
nextflow run <main_script>.nf [options]
```

See each subdirectory's README for pipeline-specific commands and parameters.
