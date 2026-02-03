#!/usr/bin/env python
"""Prepare AlphaFold 3 JSON input files for refolding designed complexes.

Reads PDB files from a sidechain-packing directory, extracts per-chain
sequences using BioPython, and writes one AF3-format JSON per design.

Usage:
    python pipeline_scripts/prepare_af3_refold_jsons.py \
        --pdb_dir 05_sidechain_packing \
        --output_dir 07_af3_refold/json \
        --num_seeds 20
"""
import argparse
import json
import os
import sys
from glob import glob

from Bio.PDB import PDBParser

PROTEIN_LETTERS_3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E",
    "PHE": "F", "GLY": "G", "HIS": "H", "ILE": "I",
    "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
    "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S",
    "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M",
}


def extract_chain_sequences(pdb_path):
    """Return {chain_id: sequence} for all chains in a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    chain_seqs = {}
    for chain in structure[0]:
        seq = []
        for residue in chain:
            if residue.id[0] == " ":  # skip HETATMs / water
                aa = PROTEIN_LETTERS_3TO1.get(residue.get_resname().upper(), "X")
                seq.append(aa)
        chain_seqs[chain.id] = "".join(seq)
    return chain_seqs


def build_af3_json(name, chain_seqs, model_seeds):
    """Build an AF3-format JSON dict for a single design."""
    sequences = []
    for chain_id, seq in sorted(chain_seqs.items()):
        sequences.append({
            "protein": {
                "id": [chain_id],
                "sequence": seq,
            }
        })
    return {
        "name": name,
        "modelSeeds": model_seeds,
        "sequences": sequences,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Prepare AF3 JSON inputs from designed PDB complexes."
    )
    parser.add_argument(
        "--pdb_dir", required=True,
        help="Directory containing designed PDB files (e.g. 05_sidechain_packing)."
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Directory to write AF3 JSON files."
    )
    parser.add_argument(
        "--num_seeds", type=int, default=20,
        help="Number of AF3 model seeds (default: 20)."
    )
    parser.add_argument(
        "--chains", default="A,C",
        help="Comma-separated chain IDs to include (default: A,C)."
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    wanted_chains = set(args.chains.split(","))
    model_seeds = list(range(1, args.num_seeds + 1))

    pdb_files = sorted(glob(os.path.join(args.pdb_dir, "*.pdb")))
    if not pdb_files:
        print(f"No PDB files found in {args.pdb_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(pdb_files)} PDB files from {args.pdb_dir}")

    for pdb_path in pdb_files:
        stem = os.path.splitext(os.path.basename(pdb_path))[0]
        chain_seqs = extract_chain_sequences(pdb_path)

        # Keep only wanted chains
        chain_seqs = {c: s for c, s in chain_seqs.items() if c in wanted_chains}
        if not chain_seqs:
            print(f"  WARNING: {stem} has no chains matching {wanted_chains}, skipping")
            continue

        data = build_af3_json(stem, chain_seqs, model_seeds)
        out_path = os.path.join(args.output_dir, f"{stem}.json")
        with open(out_path, "w") as f:
            json.dump(data, f, indent=2)

    written = len(glob(os.path.join(args.output_dir, "*.json")))
    print(f"Wrote {written} JSON files to {args.output_dir}")


if __name__ == "__main__":
    main()
