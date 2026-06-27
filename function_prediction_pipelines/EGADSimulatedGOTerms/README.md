# EGADSimulatedGOTerms

The core simulation pipeline. Generates **pseudo-bulk RNA-seq datasets** by mixing scRNA-seq cell type profiles at stochastically controlled cell type compositions, then runs EGAD to measure how composition affects GO term learnability.

This pipeline can be run with either **real GO terms** (Biological Process annotations) or **simulated GO terms** (gene sets constructed to be enriched in specific cell types), making it possible to directly test whether compositional bias drives prediction performance.

---

## What It Does

1. **Compute cell type profiles** — average expression per cell type from scRNA-seq h5ad
2. **Simulate pseudo-bulk samples** — mix cell type profiles at proportions drawn from a distribution parameterized by a baseline composition and a variance factor
3. **Run EGAD** — build co-expression networks from simulated bulk matrices and predict GO term membership
4. **Collect results** — melt all bootstrap × variance × tissue EGAD outputs into one master dataframe

---

## Pipeline Steps

```
scRNA-seq h5ad
    └─ makeCTProfiles          Compute cell type profiles (one row per cell type)
         └─ simBulk            Simulate N pseudo-bulk samples at K variance levels
              ├─ getSimulationStats    Summarize composition distributions (QC)
              ├─ graphAverageComposition   Plot average cell type proportions (QC)
              └─ sim_bulk_EGAD      Build co-expression network + EGAD (simBulkEgad.R)
                   └─ makeMeltedDF       Melt per-bootstrap EGAD results
                        └─ makeALLMeltedDF   Concatenate all bootstraps → master_melted_df.csv.gz
```

---

## Key Scripts

| Script | Language | Description |
|--------|----------|-------------|
| `sim_pipe.nf` | Nextflow | Main pipeline definition |
| `bin/makeCellTypeProfiles.py` | Python | Computes mean expression per cell type from scRNA-seq h5ad |
| `bin/simBulk.py` | Python | Simulates pseudo-bulk by scaling cell type profiles by stochastic proportions |
| `bin/simBulkEgad.R` | R | Builds co-expression network and runs EGAD neighbor voting |
| `bin/makeMeltedDf.py` | Python | Melts per-variance EGAD results into a long-format CSV |
| `bin/makeMasterMeltedDf.py` | Python | Concatenates all per-bootstrap melted DFs into the master output |
| `bin/getSimulationStats.py` | Python | Computes average cell type proportions across simulations |
| `bin/graphSparsity.py` | Python | QC: plots expression sparsity of cell type profiles |
| `bin/graphSimulationComposition.py` | Python | QC: plots per-variance composition consistency |
| `bin/graphAverageComposition.py` | Python | QC: plots average composition across variance levels |
| `bin/graphEGAD_avg.py` | Python | Plots mean EGAD AUROC by variance level |
| `bin/graphEGAD_splitCTAffiliation.py` | Python | Plots AUROC stratified by GO term cell type affiliation |

---

## Bulk Simulation Method (`simBulk.py`)

For each simulated bulk sample:
1. Draw per-cell-type proportions from a **normal distribution** centered on the baseline composition from `cell_type_proportions.json`, with standard deviation scaled by `variance_factor`
2. Clip negative proportions to 0; rescale to sum to `totalSampleSize` (default: 1000 cells)
3. Multiply each cell type profile by its assigned count → uncollapsed simulated bulk
4. Sum across cell types → one collapsed bulk sample (a single gene expression vector)
5. Repeat N times → a simulated bulk dataset (samples × genes)

---

## Configuration (`nextflow.config`)

| Parameter | Description |
|-----------|-------------|
| `params.tissue_dir` | Path to scRNA-seq h5ad file |
| `params.cell_type_proportions` | JSON with baseline mean/stdev per cell type |
| `params.compositional_variance` | List of variance factors (e.g., `[0.05, 0.1, 0.5, 1]`) |
| `params.go_annotations` | Path to GO BP annotation CSV |
| `params.gene_column` | Gene name column in GO file |
| `params.num_sims` | Number of simulated bulk samples per dataset |
| `params.num_bootstrap` | Number of full bootstrap iterations |
| `params.publish` | Output directory |

### Cell Type Proportion JSON Format

```json
{
  "brain_sc_with_metadata": {
    "Excitatory neurons": [0.45, 0.08],
    "Inhibitory neurons": [0.12, 0.03],
    "Astrocytes": [0.10, 0.02],
    ...
  }
}
```

Each entry is `[mean_proportion, stdev]`. The tissue key is matched by prefix to the h5ad filename.

### Variant Configs

| Config file | Description |
|-------------|-------------|
| `allcells_realgo.config` | All cell types, real GO BP terms |
| `allcells_simgo.config` | All cell types, simulated GO terms (excitatory neuron-enriched) |
| `allcells_simgobrain.config` | All cell types, simulated GO terms (brain-enriched) |
| `noexi_simgo.config` | Excitatory neurons removed, simulated GO terms |
| `noneur_realgo.config` | All neurons removed, real GO terms |

---

## Running

### Bootstrap mode

```bash
nextflow run sim_pipe.nf -entry bootstrap
```

### With a specific config variant

```bash
nextflow run sim_pipe.nf -entry bootstrap -c allcells_simgo.config
```

---

## Outputs

```
<publish>/
├── master_melted_df.csv.gz          All results merged (used for plotting)
└── <bootstrap>/
    ├── CTProfiles/
    │   └── exp_*_cell_type_profiles.csv
    ├── simulations/<variance>/<profile>/
    │   ├── *_.csv.gz                     Simulated bulk dataset
    │   └── *_n_sim_*_profiles.csv        Per-simulation cell type composition
    ├── EGAD/<profile>/
    │   └── *_EGAD.csv                    EGAD AUROC results
    ├── EGAD/melted_dfs/
    │   └── *_melted_df.csv.gz
    ├── stats/avg_composition/
    └── graphs/
```
