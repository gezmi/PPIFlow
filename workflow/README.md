# PPIFlow Automated Pipeline

Fully automated Snakemake pipeline: one command in, final ranking out. Six modes:

- **`target_affinity`** -- Affinity maturation of an existing nanobody-receptor complex
- **`target_full`** -- Same as affinity but with AF3Score filtering before relax
- **`target_denovo`** -- De novo nanobody design from a receptor PDB
- **`target_rescore`** -- Rescore AF3 outputs with AF3Score + Rosetta relax
- **`target_nanobody_full`** -- Two-round de novo nanobody design (paper B.2)
- **`target_binder_full`** -- Two-round de novo binder design (paper B.1)

AF3 refolding is integrated directly into the pipeline (no manual steps).

## Quick Start

### Affinity maturation (most common)

**On a compute node** (GPU available):
```bash
snakemake --snakefile workflow/Snakefile \
  --configfile workflow/config_default.yaml \
  --config input_file=/path/to/complex.cif output_dir=/path/to/output \
    af3_home=/path/to/af3 \
  target_affinity
```

**From a login node** (LSF job submission):
```bash
snakemake --snakefile workflow/Snakefile \
  --configfile workflow/config_default.yaml \
  --config input_file=/path/to/complex.cif output_dir=/path/to/output \
    af3_home=/path/to/af3 \
  --cluster "bash workflow/lsf_submit.sh {resources.gpu} {resources.mem_mb} {threads} {rule}" \
  --default-resources gpu=0 mem_mb=8000 \
  --latency-wait 120 --jobs 20 \
  target_affinity
```

### De novo nanobody design

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile workflow/config_default.yaml \
  --config input_file=/path/to/receptor.pdb output_dir=/path/to/output \
    mode=denovo af3_home=/path/to/af3 \
    antigen_chain=C specified_hotspots="C100,C120" \
  target_denovo
```

### Full pipeline (with AF3Score)

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile workflow/config_default.yaml \
  --config input_file=/path/to/complex.cif output_dir=/path/to/output \
    mode=full af3_home=/path/to/af3 af3score_dir=/path/to/af3score \
  target_full
```

### Rescore AF3 outputs

Takes a directory of AF3 output directories, converts CIF to PDB, runs AF3Score + Rosetta relax in parallel, and produces a ranked CSV. No DockQ.

**From a login node** (LSF job submission):
```bash
snakemake -s workflow/Snakefile --configfile workflow/config_default.yaml \
  --config input_dir=/path/to/af3_output output_dir=/path/to/rescore_output \
    mode=rescore af3score_dir=/path/to/af3score \
  --cluster "bash workflow/lsf_submit.sh {resources.gpu} {resources.mem_mb} {threads} {rule}" \
  --default-resources gpu=0 mem_mb=8000 \
  --latency-wait 120 --jobs 20 \
  -- target_rescore
```

By default, converts the top model per design (`{name}_model.cif`). Use `rescore_samples=all` to process all seed-samples, or `rescore_samples=seed-10_sample-0` for a specific one.

### Two-round nanobody design (paper B.2)

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile workflow/config_nanobody_full.yaml \
  --config input_file=/path/to/receptor.pdb output_dir=/path/to/output \
    af3score_dir=/path/to/af3score af3_home=/path/to/af3 \
  -- target_nanobody_full
```

### Two-round binder design (paper B.1)

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile workflow/config_binder_full.yaml \
  --config input_file=/path/to/receptor.pdb output_dir=/path/to/output \
    af3score_dir=/path/to/af3score af3_home=/path/to/af3 \
  -- target_binder_full
```

## Pipeline DAGs

### Two-round nanobody/binder mode

```
00_input (receptor PDB) →
  01_denovo [GPU] (backbone generation) →
  01b_filtered_backbones (CDR interface ratio filter) →
  02_position_analysis_r1 (CDR/framework extraction) →
  03_seqdesign_r1 [GPU] (AbMPNN/ProteinMPNN) →
  04_sidechain_packing_r1 [GPU] (FlowPacker) →
  05_af3score_r1 [GPU] (ipTM > 0.2 filter) →
  06_enrichment (interface energy → key contacts) →
  07_partial_flow [GPU×N] (per-backbone refinement) →
  08_position_analysis_r2 →
  09_seqdesign_r2 [GPU] →
  10_sidechain_packing_r2 [GPU] (FlowPacker) →
  11_af3score_r2 [GPU] (ipTM > 0.5, pTM > 0.8) →
  12_rosetta_relax (interface or full) →
  13_af3_refold [GPU×N] →
  14_dockq →
  final_ranking.csv
```

FlowPacker requires `flowpacker_dir` in config. Falls back to PyRosetta if not set.

### Affinity mode

```
00_input (CIF→PDB, chain rename) →
  01_rosetta_interface (per-residue energy) →
  02_cdr_analysis (CDR extraction, fixed positions) →
  03_partial_flow [GPU] (backbone generation) →
  04_abmpnn [GPU] (sequence design) →
  05_sidechain_packing (thread + repack) →
  06_rosetta_relax (energy minimization) ──────────────────┐
  07_af3_refold (prepare JSONs → refold [GPU×N] → CIF→PDB) →
  08_dockq (structural comparison) →
  final_ranking.csv
```

### De novo mode

```
00_input (receptor PDB) →
  01_denovo [GPU] (CDR generation) →
  02_cdr_analysis (framework indices) →
  03_abmpnn [GPU] (sequence design) →
  04_sidechain_packing (thread + repack) →
  05_rosetta_relax (energy minimization) ──────────────────┐
  05_af3_refold (prepare JSONs → refold [GPU×N] → CIF→PDB) →
  06_dockq (structural comparison) →
  final_ranking.csv
```

### Full mode (with AF3Score)

Same as affinity, but inserts `06_af3score` filtering after sidechain_packing. Relax shifts to `07_`, AF3 refold to `08_`, DockQ to `09_`.

### Rescore mode

```
00_pdbs (CIF→PDB conversion) ─┬─→ 01_af3score (prepare → jax [GPU] →
                               │     inference [GPU] → metrics)
                               │
                               └─→ 02_rosetta_relax [20 CPUs]
                                              ↓
                                rescore_ranking.csv  ← merges both
```

Rosetta relax and AF3Score run in parallel (both only need CIF conversion done). Ranking = `iptm * 100 - interface_score`. Filters: iptm > 0.5, ptm_A > 0.8.

## Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input_file` | *(required)* | Input CIF or PDB file |
| `output_dir` | *(required)* | Output directory |
| `af3_home` | *(required)* | Path to AF3 installation (conda env) |
| `mode` | `affinity` | Pipeline mode: `denovo`, `affinity`, or `full` |
| `nanobody_chain` | `B` | Nanobody chain ID in input (affinity/full) |
| `receptor_chain` | `A` | Receptor chain ID in input (affinity/full) |
| `energy_threshold` | `-5.0` | Key contact energy cutoff (REU) |
| `start_t` | `0.6` | Partial flow start time (0=full noise, 1=original) |
| `samples_per_target` | `8` | Backbone samples to generate |
| `num_seq_per_target` | `4` | AbMPNN sequences per backbone |
| `sampling_temp` | `0.1` | AbMPNN sampling temperature |
| `af3_refold_seeds` | `20` | AF3 model seeds for refolding |
| `af3_model_dir` | | AF3 model weights (defaults to AF3 internal) |

### De novo mode only

| Parameter | Default | Description |
|-----------|---------|-------------|
| `framework_pdb` | `Framework/7xl0_nanobody_framework.pdb` | Framework PDB |
| `antigen_chain` | `C` | Antigen chain in input PDB |
| `heavy_chain` | `H` | Heavy chain in framework PDB |
| `cdr_length` | `CDRH1,5-12,...` | CDR length ranges |
| `specified_hotspots` | | Target residues (e.g. `C100,C120`) |

### Full mode only

| Parameter | Default | Description |
|-----------|---------|-------------|
| `af3score_dir` | | AF3Score installation path |
| `af3score_num_jobs` | `4` | AF3Score batch jobs |

### Two-round mode only (`nanobody_full` / `binder_full`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `flowpacker_dir` | | Path to FlowPacker installation (GPU sidechain packing) |
| `flowpacker_ckpt` | `cluster.pth` | FlowPacker checkpoint name |
| `cdr_interface_ratio_threshold` | `0.6` | Min CDR interface ratio for backbone filtering |
| `r1_samples_per_target` | `100` | R1 backbone samples (scale to 25K for production) |
| `r1_num_seq_per_target` | `8` | AbMPNN sequences per R1 backbone |
| `r1_sampling_temp` | `0.5` (nb) / `0.2` (binder) | R1 sequence design temperature |
| `r2_start_t` | `0.6` | Partial flow conservativeness |
| `r2_samples_per_backbone` | `8` | R2 samples per enriched backbone |
| `relax_method` | `interface` | `"interface"` (fast) or `"full"` (unconstrained) |
| `af3score_dir` | *(required)* | AF3Score installation path |
| `af3_home` | *(required)* | AF3 installation path |

### Rescore mode only

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input_dir` | *(required)* | Directory of AF3 output dirs or flat PDB/CIF files |
| `rescore_samples` | `top` | `"top"` (top model), `"all"` (every seed-sample), or specific e.g. `"seed-10_sample-0"` |
| `af3score_dir` | *(required)* | AF3Score installation path |

## Cluster Execution

The `workflow/lsf_submit.sh` wrapper submits each Snakemake rule as an LSF job:
- Rules declaring `resources: gpu=1` get a GPU on `short-gpu`
- CPU-only rules go to `short`
- AF3 refold jobs (`af3_refold_single`) each get their own GPU job

With `--jobs 20`, multiple AF3 refold jobs run in parallel on different GPUs.

## Filtering Criteria

**Affinity/De novo** (`final_ranking.csv`):
- DockQ > 0.49
- Score = DockQ * 100 - interface_score (higher is better)

**Full pipeline** (`final_ranking.csv`):
- DockQ > 0.49, iptm > 0.7, ptm_A > 0.8
- Score = iptm * 100 - interface_score

## Hotspot-Guided Patch Design (Surface Scanning)

To scan an entire target surface, use `pipeline_scripts/split_to_patches.py` to generate hotspot patches, then run one de novo pipeline per patch:

```bash
# Generate hotspot patches from target surface
python pipeline_scripts/split_to_patches.py -i target.pdb -c C -o patches.txt

# Each line in patches.txt is a hotspot string (e.g. "C100,C120,C145")
# Launch one pipeline per patch — see results/260317_rbd_hotspot_patches/launch_all.sh
```

This is useful for exploring which epitope regions are designable before committing to a full campaign.

## Cluster Execution Notes

When using `--cluster` with `lsf_submit.sh`, always include `--default-resources gpu=0` because some Snakemake rules do not define GPU resources. Without it, `{resources.gpu}` is undefined and the submit script fails:

```bash
--cluster "bash workflow/lsf_submit.sh {resources.gpu} {resources.mem_mb} {threads} {rule}" \
--default-resources gpu=0 mem_mb=8000
```

Also note that `--rerun-triggers` uses `nargs='+'`, which means it will consume the target rule name if placed right before it. Use `--` to separate flags from the target:

```bash
snakemake ... --rerun-triggers mtime -- target_affinity
```

## Directory Structure (two-round nanobody_full mode)

```
output_dir/
  00_input/                       # Receptor PDB (chain C)
  01_denovo/                      # Generated backbone PDBs + sample_metrics.csv [GPU]
  01b_filtered_backbones/         # PDBs passing CDR interface ratio filter
  02_position_analysis_r1/        # CDR/framework extraction
  03_seqdesign_r1/                # R1 AbMPNN sequences (.fa files)
  04_sidechain_packing_r1/        # FlowPacker output PDBs [GPU]
    flowpacker_input.csv          # Threading CSV
    flowpacker_raw/               # Raw FlowPacker output (run_1/)
    {backbone}_seq{N}.pdb         # Renamed for enrichment
  05_af3score_r1/                 # R1 AF3Score filtering
    filtered_links/               # PDBs passing ipTM threshold
  06_enrichment/                  # Key contact aggregation
    {backbone}.json               # Per-backbone enriched positions
    enrichment_manifest.csv       # Summary
  07_partial_flow/                # R2 backbone refinement [GPU×N]
    per_backbone/{name}/          # Per-backbone output
    gathered/                     # Flat collection of all R2 backbones
  08_position_analysis_r2/        # R2 fixed positions
  09_seqdesign_r2/                # R2 AbMPNN sequences
  10_sidechain_packing_r2/        # R2 FlowPacker output [GPU]
  11_af3score_r2/                 # R2 AF3Score filtering
    filtered_links/               # PDBs passing stricter thresholds
  12_rosetta_relax/               # Energy minimization
  13_af3_refold/                  # AF3 structural validation [GPU×N]
  14_dockq/                       # DockQ scoring
  final_ranking.csv               # Ranked designs
```

## Directory Structure (affinity mode)

```
output_dir/
  00_input/                 # Preprocessed PDB (A=nanobody, C=receptor)
  01_rosetta_interface/     # Per-residue interface energy
  02_cdr_analysis/          # CDR positions, fixed residues
  03_partial_flow/          # Generated backbone PDBs [GPU]
  04_abmpnn/                # Redesigned sequences [GPU]
  05_sidechain_packing/     # Threaded + repacked PDBs
  06_rosetta_relax/         # Energy-minimized structures
  07_af3_refold/
    json/                   # AF3 input JSONs (one per design)
    af3_raw/                # AF3 raw output (CIF)
    pdb/                    # Converted PDBs (one subdir per design)
  08_dockq/                 # DockQ scores
  final_ranking.csv         # Ranked designs
```
