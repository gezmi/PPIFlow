#!/usr/bin/env python
"""Convert AF3 CIF outputs to PDB for a single design.

AF3 creates timestamped output directories:
  af3_raw/{design}/{design}_{timestamp}/

This script finds those directories, converts all *_model.cif files
to PDB format using BioPython, and writes them to a flat output directory.

Usage:
    python pipeline_scripts/convert_af3_output.py \
        --af3_output_dir 07_af3_refold/af3_raw \
        --design sample0_seq1 \
        --output_dir 07_af3_refold/pdb/sample0_seq1
"""
import argparse
import os
import sys
from glob import glob

from Bio.PDB import MMCIFParser, PDBIO


def convert_cif_to_pdb(cif_path, pdb_path):
    """Convert a single CIF file to PDB format."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("s", cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)


def main():
    parser = argparse.ArgumentParser(
        description="Convert AF3 CIF outputs to PDB for a single design."
    )
    parser.add_argument(
        "--af3_output_dir", required=True,
        help="Parent dir containing AF3 output (e.g. af3_raw/)."
    )
    parser.add_argument(
        "--design", required=True,
        help="Design name (e.g. sample0_seq1)."
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Where to write PDB files."
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # AF3 output structure: {af3_output_dir}/{design}/{design}_{timestamp}/
    design_dir = os.path.join(args.af3_output_dir, args.design)
    if not os.path.isdir(design_dir):
        print(f"ERROR: Design directory not found: {design_dir}", file=sys.stderr)
        sys.exit(1)

    # Find all CIF model files across timestamped subdirs
    cif_files = sorted(glob(os.path.join(design_dir, "**", "*_model.cif"),
                            recursive=True))
    if not cif_files:
        # Fallback: try any .cif file
        cif_files = sorted(glob(os.path.join(design_dir, "**", "*.cif"),
                                recursive=True))

    if not cif_files:
        print(f"ERROR: No CIF files found in {design_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Converting {len(cif_files)} CIF files for design '{args.design}'")

    for cif_path in cif_files:
        # Extract seed and sample info from the parent directory name
        # AF3 structure: {design}_{timestamp}/seed-{N}_sample-{M}/{name}_model.cif
        parts = cif_path.split(os.sep)
        # Find the seed-N_sample-M directory
        seed_sample = None
        for part in parts:
            if part.startswith("seed-") and "_sample-" in part:
                seed_sample = part
                break

        if seed_sample:
            pdb_name = f"{args.design}_{seed_sample}.pdb"
        else:
            # Fallback: use CIF filename
            cif_stem = os.path.splitext(os.path.basename(cif_path))[0]
            pdb_name = f"{cif_stem}.pdb"

        pdb_path = os.path.join(args.output_dir, pdb_name)
        try:
            convert_cif_to_pdb(cif_path, pdb_path)
        except Exception as e:
            print(f"  WARNING: Failed to convert {cif_path}: {e}", file=sys.stderr)
            continue

    written = len(glob(os.path.join(args.output_dir, "*.pdb")))
    print(f"Wrote {written} PDB files to {args.output_dir}")


if __name__ == "__main__":
    main()
