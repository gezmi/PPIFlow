# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PPIFlow is a flow-matching-based framework for de novo generation of high-affinity protein binders targeting precise epitopes. It uses discrete-continuous flow matching for backbone generation and affinity maturation (sequence refinement). Supported tasks: PPI binders, nanobodies (VHH), antibodies, motif scaffolding, and unconditional monomer generation.

## Environment Setup

```bash
conda env create -f environment.yml
conda activate ppiflow
```

Key stack: Python 3.10, PyTorch 2.3.1+cu121, PyTorch Lightning 2.5.0, PyG 2.6.1, Hydra/OmegaConf for config, PyRosetta for structure analysis. The `env/` directory contains the full conda environment and should be excluded from searches.

Model checkpoints (downloaded separately from Google Drive) go in `checkpoints/`: `binder.ckpt`, `antibody.ckpt`, `nanobody.ckpt`, `monomer.ckpt`.

## Running Inference

Entry points are the `sample_*.py` scripts at the repo root:

```bash
# Binder design
python sample_binder.py --input_pdb target.pdb --target_chain B --binder_chain A \
    --config configs/inference_binder.yaml --model_weights checkpoints/binder.ckpt \
    --specified_hotspots "B119,B141,B200" --samples_per_target 5 --output_dir output/

# Antibody/Nanobody CDR design
python sample_antibody_nanobody.py --antigen_pdb antigen.pdb --framework_pdb framework.pdb \
    --antigen_chain C --heavy_chain A --light_chain B \
    --config configs/inference_nanobody.yaml --model_weights checkpoints/nanobody.ckpt \
    --cdr_length "CDRH1,8-8,CDRH2,8-8,CDRH3,10-20" --output_dir output/

# Partial flow (affinity maturation) - local redesign of existing structure
python sample_antibody_nanobody_partial.py --complex_pdb complex.pdb \
    --fixed_positions A97-111 --cdr_position A26-33,A51-58,A97-111 \
    --start_t 0.5 --config configs/inference_nanobody.yaml \
    --model_weights checkpoints/nanobody.ckpt --output_dir output/

# Unconditional monomer generation
python sample_monomer.py --config configs/inference_unconditional.yaml \
    --model_weights checkpoints/monomer.ckpt --length_subset "[50, 100]" --samples_num 5

# Motif scaffolding
python sample_monomer.py --config configs/inference_scaffolding.yaml \
    --model_weights checkpoints/monomer.ckpt --motif_csv motif.csv \
    --motif_names "['01_1LDB']" --samples_num 5
```

## Snakemake Automated Pipeline

The `workflow/` directory contains a Snakemake pipeline for end-to-end nanobody affinity maturation:

```bash
# Affinity maturation pipeline (no AF3Score)
snakemake --snakefile workflow/Snakefile \
    --configfile workflow/config_default.yaml \
    --config input_file=/path/to/complex.cif output_dir=/path/to/output \
    target_affinity

# Full pipeline with AF3Score
snakemake --snakefile workflow/Snakefile \
    --configfile workflow/config_default.yaml \
    --config input_file=/path/to/complex.cif output_dir=/path/to/output \
    af3score_dir=/path/to/af3score/ target_full
```

Pipeline stages (output subdirectories numbered 00-08):
1. `prepare_input` - CIF-to-PDB conversion, chain renaming (nanobody→A, receptor→C)
2. `interface_energy` - PyRosetta per-residue interface energy analysis
3. `cdr_analysis` - CDR extraction (IMGT numbering), fixed position computation
4. `partial_flow` - Backbone generation via partial flow matching (GPU)
5. `abmpnn` - Sequence design via AbMPNN with fixed framework+key contacts (GPU)
6. `sidechain_packing` - Thread sequences and repack sidechains (PyRosetta)
7. `af3score` (full pipeline only) - Structure quality assessment
8. `rosetta_relax` - Energy minimization
9. Phase 2 (`score_affinity`/`score_full`) - DockQ ranking after manual AF3 refold step

Key pipeline parameters in config YAML: `energy_threshold` (-5.0 REU), `start_t` (0.6, higher=more conservative), `samples_per_target` (8), `num_seq_per_target` (4), `sampling_temp` (0.1).

## Architecture

### Module relationships

**`core/`** - Foundation layer shared by all modules:
- `utils/rigid_utils.py`: SE(3) rigid body transforms (Rigid class) used throughout the model
- `data/data_transforms.py`: Protein coordinate transformations
- `model/primitives.py`: AlphaFold-style components (Linear, LayerNorm)
- `np/`: NumPy utilities for data preprocessing

**`data/`** - Data loading and flow interpolation:
- `datasets*.py`: Dataset classes per task (binder, antibody, monomer, partial). Load PDB structures, create batched tensors with padding/masking.
- `interpolant*.py`: Flow matching interpolant classes. Define noise schedules and interpolation between noise and structure on SO(3) (rotations) and R^3 (translations).
- `so3_utils.py`: SO(3) rotation sampling and geodesic interpolation
- `residue_constants.py`: Amino acid chemistry constants

**`models/`** - Neural network architecture:
- `flow_module_*.py`: PyTorch Lightning modules (one per task variant). Handle training/inference loops, loss computation, checkpoint loading.
- `flow_model_*.py`: Flow model wrappers that compose the network components.
- `networks.py`: Backbone network assembling PairFormer + IPA stack.
- `pairformer.py`: Attention-based pair representation updates.
- `denoising.py` / `ipa_pytorch.py`: Invariant Point Attention (IPA) for 3D structure prediction.
- `node_feature_net.py` / `edge_feature_net.py`: Input feature extraction from residue types, coordinates, timestep.

**`experiments/`** - Inference runners:
- `inference_*.py`: Each contains an `Experiment` class that loads config + checkpoint and runs sampling via Lightning's `.test()` method.

**`analysis/`** - Metrics and validation:
- `metrics.py`: Structure quality metrics (RMSD, TM-score, clash detection)
- `antibody_metric.py`: Antibody-specific CDR metrics

**`pipeline_scripts/`** - Pipeline step implementations:
- `step01_cif_to_pdb.py`, `step02_interface_energy.py`, `step06_thread_and_repack.py`
- `pipeline_common.py`: Shared utilities (chain renaming, etc.)

**`ProteinMPNN/`** - Integrated sequence design tool (AbMPNN variant)

**`demo_scripts/`** - Supporting scripts for AF3Score, Rosetta analysis, CDR extraction, DockQ scoring

### Data flow through the system

```
Input PDB/CIF → preprocessing (chain renaming, interface analysis)
  → NodeFeatureNet + EdgeFeatureNet (embed residues, pairs, timestep)
  → PairformerStack (attention on pair representations)
  → IPAStack (3D structure generation via invariant point attention)
  → Flow Interpolant (SO(3) geodesic for rotations + Gaussian for translations)
  → Generated backbone PDBs
  → AbMPNN (constrained sequence design)
  → PyRosetta (sidechain packing, relax, interface scoring)
  → AF3Score / DockQ (quality assessment and ranking)
```

### Key design patterns

- **Config management**: YAML configs loaded via OmegaConf, overridable by CLI args. Each `sample_*.py` script has a `ConfigManager` that merges config file + command-line overrides.
- **Task variants via config**: The same model architecture handles different tasks (binder, antibody, monomer, scaffolding) by switching the config `task` field and using different `flow_module_*` / `datasets_*` classes.
- **Partial flow**: Affinity maturation uses `start_t` parameter (0-1) to control how much of the original structure is preserved. Higher `start_t` = more conservative refinement.
- **Chain convention in pipeline**: Nanobody/binder is always chain A, receptor/antigen is always chain C (internal convention after prepare_input step).
- **Protein representation**: Each residue is represented by N-CA-C-CB atom coordinates, an SE(3) rigid frame, residue type, and masks. The `Rigid` class in `core/utils/rigid_utils.py` handles all rigid body math.

## Demo Notebook

`demo_vhh.ipynb` contains the complete end-to-end VHH-CCL2 nanobody design workflow and is the best reference for understanding the full pipeline flow.

## Working Conventions

### Git and completion discipline
- Ask about git commit only when a feature has been tested and proven to work.
- Never mark something as done before proving it works.
- Challenge your own work before presenting it. Ask yourself: would a staff engineer approve this?

### Lessons and self-correction
- Maintain `tasks/lessons.md`. Update it after corrections from the user.
- Write rules for yourself that prevent the same mistake from recurring.
- Ruthlessly iterate on these rules until error rate drops.
- Review lessons at session start.

### Code quality
- Prioritize maintainability over quick fixes (unless directed otherwise).
- Find root causes; do not hack around problems unless necessary.
- Make as few changes as possible. Touch minimal code. Try not to introduce bugs.
- For non-trivial challenges, pause and ask if there is a more elegant way.
- If a solution feels hacky, ask: "If I knew everything I know now, would I implement the same solution?"
- Do not overengineer. Output is good, but excessive printing is usually not necessary.

### Workflow
- Use subagents to keep the main window clean.
- Enter plan mode for non-trivial tasks (3+ steps or architectural changes) without asking.
- Simple tasks: just do them directly.
