#!/usr/bin/env python3
"""Step 6: Sequence threading and sidechain packing with PyRosetta.

Threads each AbMPNN-designed sequence onto its backbone PDB using
pyrosetta.toolbox.mutate_residue(), then runs FastRelax with backbone
fixed and sidechains free (chi optimization only).

Usage:
    python step06_thread_and_repack.py \
        --backbone_dir path/to/03_partial_flow/ \
        --fasta_dir path/to/04_abmpnn/seqs/ \
        --output_dir path/to/05_sidechain_packing/ \
        --design_chain A
"""

import argparse
import glob
import os
import re

import pyrosetta
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.toolbox import mutate_residue


AA_3TO1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}
AA_1TO3 = {v: k for k, v in AA_3TO1.items()}


def parse_fasta(fasta_path: str) -> list:
    """Parse a FASTA file and return list of (header, sequence) tuples."""
    sequences = []
    header = None
    seq_lines = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    sequences.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        sequences.append((header, "".join(seq_lines)))
    return sequences


def extract_chain_sequence(full_sequence: str, pose, chain_id: str) -> str:
    """Extract the subsequence corresponding to a specific chain from
    a concatenated multi-chain FASTA sequence."""
    pdb_info = pose.pdb_info()
    chain_start = None
    chain_end = None
    for i in range(1, pose.total_residue() + 1):
        if pdb_info.chain(i) == chain_id:
            if chain_start is None:
                chain_start = i
            chain_end = i

    if chain_start is None:
        raise ValueError(f"Chain {chain_id} not found in pose")

    # The FASTA from ProteinMPNN concatenates all chains in PDB order
    offset = chain_start - 1
    length = chain_end - chain_start + 1
    return full_sequence[offset:offset + length], chain_start, chain_end


def thread_sequence(pose, new_seq: str, chain_start: int, chain_end: int):
    """Thread a new sequence onto the pose for residues chain_start..chain_end."""
    length = chain_end - chain_start + 1
    if len(new_seq) != length:
        raise ValueError(
            f"Sequence length {len(new_seq)} != chain residue count {length}")

    for idx, resnum in enumerate(range(chain_start, chain_end + 1)):
        current_aa = AA_3TO1.get(pose.residue(resnum).name3(), "X")
        target_aa = new_seq[idx]
        if current_aa != target_aa:
            mutate_residue(pose, resnum, target_aa)


def repack_sidechains(pose, scorefxn, chain_start: int, chain_end: int,
                      max_iter: int = 200):
    """Run FastRelax with backbone fixed, sidechains free."""
    fr = FastRelax()
    fr.set_scorefxn(scorefxn)
    fr.max_iter(max_iter)

    movemap = pyrosetta.MoveMap()
    # Fix all backbone, allow chi angles
    movemap.set_bb(False)
    movemap.set_chi(True)
    fr.set_movemap(movemap)

    fr.apply(pose)


def find_backbone_for_fasta(fasta_name: str, backbone_dir: str) -> str:
    """Match a FASTA file to its source backbone PDB.

    ProteinMPNN names FASTA files after the input PDB. The FASTA file for
    backbone 'sample_0.pdb' would be 'sample_0.fa'.
    """
    # Try exact match first
    base = os.path.splitext(fasta_name)[0]
    pdb_path = os.path.join(backbone_dir, base + ".pdb")
    if os.path.exists(pdb_path):
        return pdb_path

    # Try removing common suffixes added by ProteinMPNN
    for suffix in ["_0001", "_redesigned"]:
        if base.endswith(suffix):
            pdb_path = os.path.join(backbone_dir, base[:-len(suffix)] + ".pdb")
            if os.path.exists(pdb_path):
                return pdb_path

    return None


def main():
    parser = argparse.ArgumentParser(
        description="Thread AbMPNN sequences onto backbones and repack sidechains")
    parser.add_argument("--backbone_dir", required=True,
                        help="Directory with backbone PDBs from partial flow")
    parser.add_argument("--fasta_dir", required=True,
                        help="Directory with FASTA files from AbMPNN")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for full-atom PDBs")
    parser.add_argument("--design_chain", default="A",
                        help="Chain to thread new sequences onto (default: A)")
    parser.add_argument("--max_iter", type=int, default=200,
                        help="Max iterations for FastRelax (default: 200)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    pyrosetta.init(
        "-use_input_sc -ignore_unrecognized_res "
        "-ignore_zero_occupancy false -load_PDB_components false "
        "-relax:default_repeats 2 -no_fconfig",
        silent=True
    )

    scorefxn = ScoreFunctionFactory.create_score_function("ref2015")

    fasta_files = sorted(glob.glob(os.path.join(args.fasta_dir, "*.fa")))
    if not fasta_files:
        fasta_files = sorted(glob.glob(os.path.join(args.fasta_dir, "*.fasta")))

    if not fasta_files:
        print(f"No FASTA files found in {args.fasta_dir}")
        return

    total_designed = 0

    for fasta_path in fasta_files:
        fasta_name = os.path.basename(fasta_path)
        sequences = parse_fasta(fasta_path)

        if not sequences:
            print(f"  Skipping empty FASTA: {fasta_name}")
            continue

        # The first entry in ProteinMPNN output is the native/input sequence
        # Designed sequences start from index 1
        backbone_pdb = find_backbone_for_fasta(fasta_name, args.backbone_dir)
        if backbone_pdb is None:
            print(f"  WARNING: No backbone PDB found for {fasta_name}, skipping")
            continue

        print(f"Processing {fasta_name} -> backbone {os.path.basename(backbone_pdb)}")

        for seq_idx, (header, full_seq) in enumerate(sequences):
            # Skip native sequence (first entry)
            if seq_idx == 0:
                continue

            # Extract the designed chain sequence
            pose = pose_from_pdb(backbone_pdb)
            try:
                chain_seq, chain_start, chain_end = extract_chain_sequence(
                    full_seq, pose, args.design_chain)
            except ValueError as e:
                print(f"  WARNING: {e}, skipping seq {seq_idx}")
                continue

            # Thread the new sequence
            thread_sequence(pose, chain_seq, chain_start, chain_end)

            # Repack sidechains
            repack_sidechains(pose, scorefxn, chain_start, chain_end,
                              args.max_iter)

            # Save output
            base = os.path.splitext(os.path.basename(backbone_pdb))[0]
            out_name = f"{base}_seq{seq_idx}.pdb"
            out_path = os.path.join(args.output_dir, out_name)
            pose.dump_pdb(out_path)
            total_designed += 1
            print(f"  Wrote {out_name}")

    print(f"\nDone. {total_designed} full-atom PDBs written to {args.output_dir}")


if __name__ == "__main__":
    main()
