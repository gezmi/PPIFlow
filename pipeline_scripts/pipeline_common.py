"""Shared utilities and step functions for PPIFlow pipeline scripts.

Both run_pipeline.py (full pipeline with AF3Score) and
run_affinity_maturation.py (no AF3Score, CIF/PDB input) import from here.
"""

import argparse
import csv
import json
import os
import shutil
import subprocess

# Root of the PPIFlow repository (one level up from pipeline_scripts/)
PPIFLOW_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def run_cmd(cmd: str, description: str, dry_run: bool = False):
    """Run a shell command, printing it first."""
    print(f"\n{'='*60}")
    print(f"  {description}")
    print(f"{'='*60}")
    print(f"$ {cmd}\n")
    if dry_run:
        print("  [DRY RUN — skipped]")
        return 0
    ret = os.system(cmd)
    if ret != 0:
        print(f"  WARNING: command exited with code {ret}")
    return ret


def write_lsf_script(path: str, job_name: str, body: str,
                      gpu: bool = False, queue: str = None,
                      conda_env: str = None, n_cpus: int = 4,
                      log_dir: str = None):
    """Write a bsub-compatible LSF submission script."""
    if queue is None:
        queue = "gpu" if gpu else "short"
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

    lines = ["#!/bin/bash"]
    lines.append(f"#BSUB -q {queue}")
    lines.append(f"#BSUB -J {job_name}")
    lines.append(f"#BSUB -n {n_cpus}")
    if gpu:
        lines.append(f'#BSUB -gpu "num=1"')
    if log_dir:
        lines.append(f"#BSUB -o {log_dir}/{job_name}-%J.out")
        lines.append(f"#BSUB -e {log_dir}/{job_name}-%J.err")
    lines.append("")
    if conda_env:
        lines.append(f"source activate {conda_env} 2>/dev/null || conda activate {conda_env}")
        lines.append("")
    lines.append(f"cd {PPIFLOW_ROOT}")
    lines.append("")
    lines.append(body)
    lines.append("")

    with open(path, "w") as f:
        f.write("\n".join(lines))
    os.chmod(path, 0o755)
    return path


def submit_lsf(script_path: str, dry_run: bool = False):
    """Submit an LSF job via bsub."""
    cmd = f"bsub < {script_path}"
    print(f"Submitting: {cmd}")
    if dry_run:
        print("  [DRY RUN — skipped]")
        return
    ret = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(ret.stdout.strip())
    if ret.returncode != 0:
        print(f"  bsub stderr: {ret.stderr.strip()}")


# ---------------------------------------------------------------------------
# Common argparse arguments
# ---------------------------------------------------------------------------

def add_common_args(parser):
    """Add arguments shared by both pipeline scripts."""
    parser.add_argument("--output_dir", required=True,
                        help="Output base directory for all pipeline results")
    parser.add_argument("--name", required=True,
                        help="Target name (used for job names and file prefixes)")

    # Chain IDs
    parser.add_argument("--nanobody_chain", default="B",
                        help="Nanobody/VHH chain ID in the input file (default: B)")
    parser.add_argument("--receptor_chain", default="A",
                        help="Receptor/antigen chain ID in the input file (default: A)")

    # Energy analysis
    parser.add_argument("--energy_threshold", type=float, default=-5.0,
                        help="Energy threshold for key contact residues (default: -5.0 REU)")

    # Partial flow
    parser.add_argument("--start_t", type=float, default=0.6,
                        help="Partial flow start_t (default: 0.6)")
    parser.add_argument("--samples_per_target", type=int, default=8,
                        help="Number of backbone samples (default: 8)")

    # AbMPNN
    parser.add_argument("--num_seq_per_target", type=int, default=4,
                        help="AbMPNN sequences per backbone (default: 4)")
    parser.add_argument("--sampling_temp", type=float, default=0.1,
                        help="AbMPNN sampling temperature (default: 0.1)")

    # LSF / environment
    parser.add_argument("--gpu_queue", default="gpu",
                        help="LSF queue for GPU jobs (default: gpu)")
    parser.add_argument("--cpu_queue", default="short",
                        help="LSF queue for CPU jobs (default: short)")
    parser.add_argument("--conda_env", default="",
                        help="Conda environment name to activate in LSF scripts")

    # Execution control
    parser.add_argument("--dry_run", action="store_true",
                        help="Print commands without executing them")


# ---------------------------------------------------------------------------
# Step 1: Input preparation
# ---------------------------------------------------------------------------

def prepare_input_from_cif(args):
    """Convert CIF to PDB with chain renaming (nanobody->A, receptor->C)."""
    input_dir = os.path.join(args.output_dir, "00_input")
    os.makedirs(input_dir, exist_ok=True)
    output_pdb = os.path.join(input_dir, f"{args.name}.pdb")

    cmd = (
        f"python {PPIFLOW_ROOT}/pipeline_scripts/step01_cif_to_pdb.py "
        f"--input_cif {args.input_cif} "
        f"--output_pdb {output_pdb} "
        f"--nanobody_chain {args.nanobody_chain} "
        f"--receptor_chain {args.receptor_chain}"
    )
    run_cmd(cmd, "Step 1: CIF -> PDB + Chain Rename", args.dry_run)
    return output_pdb


def prepare_input_from_pdb(args):
    """Prepare PDB input, renaming chains if needed (nanobody->A, receptor->C)."""
    input_dir = os.path.join(args.output_dir, "00_input")
    os.makedirs(input_dir, exist_ok=True)
    output_pdb = os.path.join(input_dir, f"{args.name}.pdb")

    nb_chain = args.nanobody_chain
    rec_chain = args.receptor_chain

    if nb_chain == "A" and rec_chain == "C":
        print(f"\n{'='*60}")
        print(f"  Step 1: PDB Input (chains already A/C — copying)")
        print(f"{'='*60}")
        if args.dry_run:
            print("  [DRY RUN — skipped]")
        else:
            shutil.copy2(args.input_pdb, output_pdb)
            print(f"  Copied {args.input_pdb} -> {output_pdb}")
    else:
        print(f"\n{'='*60}")
        print(f"  Step 1: PDB Chain Rename ({nb_chain}->A, {rec_chain}->C)")
        print(f"{'='*60}")
        if args.dry_run:
            print("  [DRY RUN — skipped]")
        else:
            _rename_pdb_chains(args.input_pdb, output_pdb, nb_chain, rec_chain)
            print(f"  Wrote {output_pdb}")

    return output_pdb


def _rename_pdb_chains(input_pdb: str, output_pdb: str,
                       nanobody_chain: str, receptor_chain: str):
    """Rename and reorder chains in a PDB file so nanobody=A, receptor=C."""
    import gemmi

    st = gemmi.read_structure(input_pdb)
    model = st[0]

    nb_chain_obj = None
    rec_chain_obj = None
    for chain in model:
        if chain.name == nanobody_chain:
            nb_chain_obj = chain
        elif chain.name == receptor_chain:
            rec_chain_obj = chain

    if nb_chain_obj is None:
        raise ValueError(f"Nanobody chain '{nanobody_chain}' not found in PDB")
    if rec_chain_obj is None:
        raise ValueError(f"Receptor chain '{receptor_chain}' not found in PDB")

    nb_chain_obj.name = "_NB_"
    rec_chain_obj.name = "_REC_"
    nb_chain_obj.name = "A"
    rec_chain_obj.name = "C"

    for chain in model:
        for residue in chain:
            if chain.name == "A":
                residue.subchain = "A"
            elif chain.name == "C":
                residue.subchain = "C"

    new_model = gemmi.Model(model.name)
    for chain in model:
        if chain.name == "A":
            new_model.add_chain(chain)
    for chain in model:
        if chain.name == "C":
            new_model.add_chain(chain)
    for chain in model:
        if chain.name not in ("A", "C"):
            new_model.add_chain(chain)

    st[0] = new_model
    st.write_pdb(output_pdb)

    for chain in st[0]:
        polymer_res = [r for r in chain
                       if gemmi.find_tabulated_residue(r.name).is_amino_acid()]
        print(f"  Chain {chain.name}: {len(polymer_res)} amino acid residues")


# ---------------------------------------------------------------------------
# Step 2: Interface energy analysis
# ---------------------------------------------------------------------------

def step_interface_energy(args, input_pdb: str):
    output_dir = os.path.join(args.output_dir, "01_rosetta_interface")
    os.makedirs(output_dir, exist_ok=True)

    cmd = (
        f"python {PPIFLOW_ROOT}/pipeline_scripts/step02_interface_energy.py "
        f"--input_pdb {input_pdb} "
        f"--output_dir {output_dir} "
        f"--binder_chain A --target_chain C "
        f"--distance_cutoff 10.0"
    )
    run_cmd(cmd, "Step 2: Interface Energy Analysis (PyRosetta)", args.dry_run)
    return output_dir


# ---------------------------------------------------------------------------
# Step 3: CDR extraction + fixed positions
# ---------------------------------------------------------------------------

def step_cdr_and_fixed_positions(args, input_pdb: str, energy_dir: str):
    cdr_dir = os.path.join(args.output_dir, "02_cdr_analysis")
    os.makedirs(cdr_dir, exist_ok=True)
    input_dir = os.path.dirname(input_pdb)

    cdr_csv = os.path.join(cdr_dir, "cdr_idx.csv")
    cmd = (
        f"python {PPIFLOW_ROOT}/demo_scripts/get_cdr.py "
        f"{input_dir} {cdr_csv}"
    )
    run_cmd(cmd, "Step 3a: Extract CDR positions", args.dry_run)

    if args.dry_run:
        fixed_json = os.path.join(cdr_dir, "fixed_positions.json")
        placeholder = {
            "key_contact_residues": [],
            "partial_flow_fixed": "",
            "cdr_positions": "",
            "hotspots": "",
            "abmpnn_fw_and_key_fixed": "",
        }
        with open(fixed_json, "w") as f:
            json.dump(placeholder, f, indent=2)
        return cdr_dir

    # Parse energy CSV
    energy_csv = os.path.join(energy_dir, "residue_energy.csv")
    key_residues = []
    if os.path.exists(energy_csv):
        with open(energy_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                energy = float(row["binder_energy"])
                if energy < args.energy_threshold:
                    key_residues.append(row["residue"])
    print(f"Key contact residues (energy < {args.energy_threshold}): {key_residues}")

    # Parse CDR CSV
    fw_indices = []
    cdr_positions = ""
    if os.path.exists(cdr_csv):
        import pandas as pd
        df = pd.read_csv(cdr_csv)
        if not df.empty:
            row = df.iloc[0]
            fw_str = str(row.get("fw_index", ""))
            fw_indices = fw_str.split() if fw_str else []
            cdr_positions = str(row.get("r2_cdr_pos", ""))
    print(f"Framework indices: {len(fw_indices)} positions")
    print(f"CDR positions: {cdr_positions}")

    # Read hotspots
    hotspots_file = os.path.join(energy_dir, "hotspots.txt")
    hotspots = ""
    if os.path.exists(hotspots_file):
        with open(hotspots_file) as f:
            hotspots = f.read().strip()
    print(f"Hotspots: {hotspots}")

    # Compute fixed positions
    partial_flow_fixed = ",".join(key_residues)

    key_res_nums = []
    for r in key_residues:
        num = "".join(c for c in r if c.isdigit())
        if num:
            key_res_nums.append(int(num))

    fw_nums = []
    for idx in fw_indices:
        num = "".join(c for c in idx if c.isdigit())
        if num:
            fw_nums.append(int(num))

    abmpnn_fixed = sorted(set(fw_nums) | set(key_res_nums))
    abmpnn_fixed_str = " ".join(str(x) for x in abmpnn_fixed)

    fixed_data = {
        "key_contact_residues": key_residues,
        "partial_flow_fixed": partial_flow_fixed,
        "cdr_positions": cdr_positions,
        "hotspots": hotspots,
        "abmpnn_fw_and_key_fixed": abmpnn_fixed_str,
        "fw_indices": fw_indices,
        "key_res_nums": key_res_nums,
    }
    fixed_json = os.path.join(cdr_dir, "fixed_positions.json")
    with open(fixed_json, "w") as f:
        json.dump(fixed_data, f, indent=2)
    print(f"Wrote fixed positions to {fixed_json}")

    return cdr_dir


# ---------------------------------------------------------------------------
# Step 4: Partial flow (GPU, LSF)
# ---------------------------------------------------------------------------

def step_partial_flow(args, input_pdb: str, cdr_dir: str):
    output_dir = os.path.join(args.output_dir, "03_partial_flow")
    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(args.output_dir, "logs")

    fixed_json = os.path.join(cdr_dir, "fixed_positions.json")
    with open(fixed_json) as f:
        fixed_data = json.load(f)

    fixed_pos = fixed_data["partial_flow_fixed"]
    cdr_pos = fixed_data["cdr_positions"]
    hotspots = fixed_data["hotspots"]

    body = (
        f"python {PPIFLOW_ROOT}/sample_antibody_nanobody_partial.py \\\n"
        f"    --complex_pdb {input_pdb} \\\n"
        f"    --start_t {args.start_t} \\\n"
        f"    --fixed_positions \"{fixed_pos}\" \\\n"
        f"    --antigen_chain C \\\n"
        f"    --heavy_chain A \\\n"
        f"    --config {PPIFLOW_ROOT}/configs/inference_nanobody.yaml \\\n"
        f"    --samples_per_target {args.samples_per_target} \\\n"
        f"    --cdr_position \"{cdr_pos}\" \\\n"
        f"    --specified_hotspots \"{hotspots}\" \\\n"
        f"    --retry_Limit 10 \\\n"
        f"    --model_weights {PPIFLOW_ROOT}/checkpoints/nanobody.ckpt \\\n"
        f"    --output_dir {output_dir} \\\n"
        f"    --name {args.name}"
    )

    script = write_lsf_script(
        os.path.join(args.output_dir, "lsf_step04_partial_flow.sh"),
        job_name=f"pflow_{args.name}",
        body=body, gpu=True, queue=args.gpu_queue,
        conda_env=args.conda_env, log_dir=log_dir
    )
    print(f"Wrote LSF script: {script}")
    submit_lsf(script, args.dry_run)
    return output_dir


# ---------------------------------------------------------------------------
# Step 5: AbMPNN sequence redesign (GPU, LSF)
# ---------------------------------------------------------------------------

def step_abmpnn(args, backbone_dir: str, cdr_dir: str):
    output_dir = os.path.join(args.output_dir, "04_abmpnn")
    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(args.output_dir, "logs")

    fixed_json = os.path.join(cdr_dir, "fixed_positions.json")
    with open(fixed_json) as f:
        fixed_data = json.load(f)

    abmpnn_fixed = fixed_data["abmpnn_fw_and_key_fixed"]
    fixpos_csv = os.path.join(output_dir, "fixed_positions.csv")

    body = f"""# Generate fixed positions CSV from partial flow output
python -c "
import glob, os
backbone_dir = '{backbone_dir}'
fixpos_csv = '{fixpos_csv}'
abmpnn_fixed = '{abmpnn_fixed}'
pdbs = sorted(glob.glob(os.path.join(backbone_dir, '*.pdb')))
with open(fixpos_csv, 'w') as f:
    f.write('pdb_name,motif_index\\n')
    for pdb in pdbs:
        name = os.path.splitext(os.path.basename(pdb))[0]
        f.write(f'{{name}},{{abmpnn_fixed}}-\\n')
print(f'Wrote {{len(pdbs)}} entries to {{fixpos_csv}}')
"

# Run ProteinMPNN with abmpnn model
python {PPIFLOW_ROOT}/ProteinMPNN/protein_mpnn_run.py \\
    --path_to_model_weights "{PPIFLOW_ROOT}/ProteinMPNN/model_weights/" \\
    --model_name "abmpnn" \\
    --folder_with_pdbs_path "{backbone_dir}" \\
    --out_folder "{output_dir}" \\
    --chain_list "A" \\
    --position_list "{fixpos_csv}" \\
    --num_seq_per_target {args.num_seq_per_target} \\
    --sampling_temp {args.sampling_temp} \\
    --seed 37 \\
    --batch_size {args.num_seq_per_target} \\
    --omit_AAs C"""

    script = write_lsf_script(
        os.path.join(args.output_dir, "lsf_step05_abmpnn.sh"),
        job_name=f"abmpnn_{args.name}",
        body=body, gpu=True, queue=args.gpu_queue,
        conda_env=args.conda_env, log_dir=log_dir
    )
    print(f"Wrote LSF script: {script}")
    submit_lsf(script, args.dry_run)
    return output_dir


# ---------------------------------------------------------------------------
# Step 6: Sidechain packing (PyRosetta, LSF CPU)
# ---------------------------------------------------------------------------

def step_sidechain_packing(args, backbone_dir: str, abmpnn_dir: str):
    output_dir = os.path.join(args.output_dir, "05_sidechain_packing")
    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(args.output_dir, "logs")

    fasta_dir = os.path.join(abmpnn_dir, "seqs")

    body = (
        f"python {PPIFLOW_ROOT}/pipeline_scripts/step06_thread_and_repack.py \\\n"
        f"    --backbone_dir {backbone_dir} \\\n"
        f"    --fasta_dir {fasta_dir} \\\n"
        f"    --output_dir {output_dir} \\\n"
        f"    --design_chain A \\\n"
        f"    --max_iter 200"
    )

    script = write_lsf_script(
        os.path.join(args.output_dir, "lsf_step06_repack.sh"),
        job_name=f"repack_{args.name}",
        body=body, gpu=False, queue=args.cpu_queue,
        conda_env=args.conda_env, n_cpus=4, log_dir=log_dir
    )
    print(f"Wrote LSF script: {script}")
    submit_lsf(script, args.dry_run)
    return output_dir


# ---------------------------------------------------------------------------
# Rosetta relax + interface scoring (CPU, LSF)
# ---------------------------------------------------------------------------

def step_rosetta_relax(args, input_pdb_dir: str, output_dir: str,
                       lsf_script_name: str = "lsf_step_relax.sh"):
    """Run Rosetta relax on PDBs in input_pdb_dir.

    Args:
        input_pdb_dir: Directory containing PDBs to relax.
        output_dir: Where to write relaxed structures and scores.
        lsf_script_name: Name for the LSF submission script.
    """
    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(args.output_dir, "logs")

    body = (
        f"python {PPIFLOW_ROOT}/demo_scripts/relax_complex.py \\\n"
        f"    --pdb_dir {input_pdb_dir} \\\n"
        f"    --output_dir {output_dir} \\\n"
        f"    --batch_idx 0 \\\n"
        f"    --dump_pdb True \\\n"
        f"    --relax True \\\n"
        f"    --fixbb False \\\n"
        f"    --max_iter 170"
    )

    script = write_lsf_script(
        os.path.join(args.output_dir, lsf_script_name),
        job_name=f"relax_{args.name}",
        body=body, gpu=False, queue=args.cpu_queue,
        conda_env=args.conda_env, n_cpus=4, log_dir=log_dir
    )
    print(f"Wrote LSF script: {script}")
    submit_lsf(script, args.dry_run)
    return output_dir


# ---------------------------------------------------------------------------
# AF3 refold + DockQ (semi-manual)
# ---------------------------------------------------------------------------

def step_af3_refold_dockq(args, reference_dir: str, relax_dir: str,
                          refold_dir: str, dockq_dir: str,
                          lsf_script_name: str = "lsf_step_dockq.sh",
                          af3score_csv: str = None):
    """Set up AF3 refold + DockQ scoring and generate the final ranking script.

    Args:
        reference_dir: Directory with designed PDBs (used as AF3 input and
                       DockQ reference).
        relax_dir: Directory with Rosetta relax output (interface scores).
        refold_dir: Where AF3 refolding outputs should be placed.
        dockq_dir: Where DockQ results will be written.
        lsf_script_name: Name for the DockQ LSF script.
        af3score_csv: Path to AF3Score results CSV. If provided, final ranking
                      includes AF3Score metrics; otherwise ranks by DockQ +
                      interface score only.
    """
    os.makedirs(refold_dir, exist_ok=True)
    os.makedirs(dockq_dir, exist_ok=True)
    log_dir = os.path.join(args.output_dir, "logs")

    step_label = "AF3 Refold + DockQ (Semi-Manual)"
    print(f"\n{'='*60}")
    print(f"  {step_label}")
    print(f"{'='*60}")
    print(f"""
This step requires running AlphaFold3 refolding externally.

1. Run AF3 prediction (without template) on designs from:
   {reference_dir}
   Use 20 seeds x 5 models. Place outputs in:
   {refold_dir}

2. Then run DockQ scoring:
   python {PPIFLOW_ROOT}/demo_scripts/run_DockQv2_eachfolder_v2.py \\
       --input_dir {refold_dir} \\
       --reference_dir {reference_dir} \\
       --output_dir {dockq_dir}

3. Parse DockQ scores:
   python {PPIFLOW_ROOT}/demo_scripts/parse_dockq_scores.py {dockq_dir}

4. Final ranking:
   python {args.output_dir}/final_ranking.py
""")

    # Write DockQ LSF script
    body_dockq = (
        f"python {PPIFLOW_ROOT}/demo_scripts/run_DockQv2_eachfolder_v2.py \\\n"
        f"    --input_dir {refold_dir} \\\n"
        f"    --reference_dir {reference_dir} \\\n"
        f"    --output_dir {dockq_dir}\n\n"
        f"python {PPIFLOW_ROOT}/demo_scripts/parse_dockq_scores.py {dockq_dir}"
    )

    script = write_lsf_script(
        os.path.join(args.output_dir, lsf_script_name),
        job_name=f"dockq_{args.name}",
        body=body_dockq, gpu=False, queue=args.cpu_queue,
        conda_env=args.conda_env, log_dir=log_dir
    )
    print(f"Wrote DockQ LSF script: {script}")

    # Write final ranking script
    ranking_script = os.path.join(args.output_dir, "final_ranking.py")
    if af3score_csv:
        _write_ranking_with_af3score(ranking_script, dockq_dir, relax_dir,
                                     af3score_csv, args.output_dir)
    else:
        _write_ranking_dockq_only(ranking_script, dockq_dir, relax_dir,
                                  args.output_dir)
    os.chmod(ranking_script, 0o755)
    print(f"Wrote final ranking script: {ranking_script}")

    return dockq_dir


def _write_ranking_with_af3score(path, dockq_dir, relax_dir,
                                 af3score_csv, output_dir):
    with open(path, "w") as f:
        f.write(f"""#!/usr/bin/env python3
\"\"\"Final ranking of designs.

Filters:  DockQ > 0.49, AF3Score pTM > 0.8, AF3Score ipTM > 0.7
Ranking:  Score = (AF3Score ipTM * 100) - Interface Score
\"\"\"
import pandas as pd
import os

dockq_csv = "{dockq_dir}/summary_dockq_scores.csv"
relax_dir = "{relax_dir}"
af3score_csv = "{af3score_csv}"

df_dockq = pd.read_csv(dockq_csv)

relax_csvs = [f for f in os.listdir(relax_dir) if f.startswith("rosetta_complex_")]
df_relax = pd.concat([pd.read_csv(os.path.join(relax_dir, f)) for f in relax_csvs])

df_af3 = pd.read_csv(af3score_csv)

df = df_dockq.merge(df_relax[["pdb_name", "interface_score"]],
                     left_on="FolderName", right_on="pdb_name", how="left")
df = df.merge(df_af3[["description", "iptm", "ptm_A"]],
              left_on="FolderName", right_on="description", how="left")

filtered = df[
    (df["Overall_Avg_DockQ"] > 0.49) &
    (df["ptm_A"] > 0.8) &
    (df["iptm"] > 0.7)
].copy()

filtered["score"] = (filtered["iptm"] * 100) - filtered["interface_score"]
filtered = filtered.sort_values("score", ascending=False)

output_csv = "{output_dir}/final_ranked_designs.csv"
filtered.to_csv(output_csv, index=False)
print(f"{{len(filtered)}} designs passed all filters")
print(f"Results saved to {{output_csv}}")
if not filtered.empty:
    print("\\nTop designs:")
    print(filtered[["FolderName", "Overall_Avg_DockQ", "iptm", "ptm_A",
                     "interface_score", "score"]].head(10).to_string(index=False))
""")


def _write_ranking_dockq_only(path, dockq_dir, relax_dir, output_dir):
    with open(path, "w") as f:
        f.write(f"""#!/usr/bin/env python3
\"\"\"Final ranking of designs.

Filters:  DockQ > 0.49
Ranking:  Score = (DockQ * 100) - Interface Score
\"\"\"
import pandas as pd
import os

dockq_csv = "{dockq_dir}/summary_dockq_scores.csv"
relax_dir = "{relax_dir}"

df_dockq = pd.read_csv(dockq_csv)

relax_csvs = [f for f in os.listdir(relax_dir) if f.startswith("rosetta_complex_")]
df_relax = pd.concat([pd.read_csv(os.path.join(relax_dir, f)) for f in relax_csvs])

df = df_dockq.merge(df_relax[["pdb_name", "interface_score"]],
                     left_on="FolderName", right_on="pdb_name", how="left")

filtered = df[df["Overall_Avg_DockQ"] > 0.49].copy()

filtered["score"] = (filtered["Overall_Avg_DockQ"] * 100) - filtered["interface_score"]
filtered = filtered.sort_values("score", ascending=False)

output_csv = "{output_dir}/final_ranked_designs.csv"
filtered.to_csv(output_csv, index=False)
print(f"{{len(filtered)}} designs passed DockQ > 0.49 filter")
print(f"Results saved to {{output_csv}}")
if not filtered.empty:
    print("\\nTop designs:")
    print(filtered[["FolderName", "Overall_Avg_DockQ",
                     "interface_score", "score"]].head(10).to_string(index=False))
""")
