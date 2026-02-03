# PPIFlow Automated Pipeline

Fully automated Snakemake pipeline: one command in, final ranking out. Three modes:

- **`target_affinity`** -- Affinity maturation of an existing nanobody-receptor complex
- **`target_full`** -- Same as affinity but with AF3Score filtering before relax
- **`target_denovo`** -- De novo nanobody design from a receptor PDB

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

## Pipeline DAGs

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
