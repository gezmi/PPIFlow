#!/usr/bin/env python3
"""Step 2: Interface energy analysis using PyRosetta.

Replaces the Rosetta XML pipeline (interface_analysis/0-run_rosetta.sh) with
a pure PyRosetta approach.  Computes per-residue pair energies from the
EnergyGraph for inter-chain pairs (binder A <-> target C) within a distance
cutoff, identifies hotspot residues on the antigen, and outputs:

  - residue_energy.csv  (per-binder-residue summed interaction energies)
  - hotspots.txt         (comma-separated antigen hotspot residue IDs)

Usage:
    python step02_interface_energy.py \
        --input_pdb path/to/complex.pdb \
        --output_dir path/to/01_rosetta_interface/ \
        --binder_chain A \
        --target_chain C \
        --distance_cutoff 10.0
"""

import argparse
import csv
import os
from collections import defaultdict

import pyrosetta
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory


def get_interface_energies(pose, scorefxn, binder_chain: str, target_chain: str,
                           distance_cutoff: float = 10.0):
    """Extract per-residue pair energies across binder-target interface.

    Returns:
        binder_energies: dict mapping binder residue labels (e.g. 'A42') to
            summed two-body interaction energy with target residues
        target_energies: dict mapping target residue labels to summed
            interaction energy with binder residues
    """
    # Score the pose to populate the EnergyGraph
    scorefxn(pose)
    energies = pose.energies()
    energy_graph = energies.energy_graph()
    pdb_info = pose.pdb_info()

    binder_energies = defaultdict(float)
    target_energies = defaultdict(float)

    # Identify binder and target residue indices (1-based Rosetta numbering)
    binder_residues = []
    target_residues = []
    for i in range(1, pose.total_residue() + 1):
        chain = pdb_info.chain(i)
        if chain == binder_chain:
            binder_residues.append(i)
        elif chain == target_chain:
            target_residues.append(i)

    binder_set = set(binder_residues)
    target_set = set(target_residues)

    # Iterate over all binder–target pairs and look up edges in the EnergyGraph
    for i in binder_residues:
        for j in target_residues:
            edge = energy_graph.find_energy_edge(i, j)
            if edge is None:
                continue

            # Check distance between CA atoms
            if pose.residue(i).has("CA") and pose.residue(j).has("CA"):
                ca_i = pose.residue(i).xyz("CA")
                ca_j = pose.residue(j).xyz("CA")
                dist = ca_i.distance(ca_j)
                if dist > distance_cutoff:
                    continue

            # Sum the two-body energy for this pair
            emap = edge.fill_energy_map()
            pair_energy = emap.dot(scorefxn.weights())

            binder_label = f"{pdb_info.chain(i)}{pdb_info.number(i)}"
            target_label = f"{pdb_info.chain(j)}{pdb_info.number(j)}"

            binder_energies[binder_label] += pair_energy
            target_energies[target_label] += pair_energy

    return dict(binder_energies), dict(target_energies)


def detect_hotspots(target_energies: dict, energy_threshold: float = -1.0) -> list:
    """Identify target residues at the interface with strong interactions."""
    hotspots = []
    for res_label, energy in sorted(target_energies.items(),
                                     key=lambda x: x[1]):
        if energy < energy_threshold:
            hotspots.append(res_label)
    return hotspots


def main():
    parser = argparse.ArgumentParser(
        description="PyRosetta per-residue interface energy analysis")
    parser.add_argument("--input_pdb", required=True,
                        help="Input PDB file (binder + target complex)")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for energy CSV and hotspots")
    parser.add_argument("--binder_chain", default="A",
                        help="Chain ID of the binder/nanobody (default: A)")
    parser.add_argument("--target_chain", default="C",
                        help="Chain ID of the target/receptor (default: C)")
    parser.add_argument("--distance_cutoff", type=float, default=10.0,
                        help="CA-CA distance cutoff in Angstroms (default: 10.0)")
    parser.add_argument("--hotspot_energy_threshold", type=float, default=-1.0,
                        help="Energy threshold for hotspot detection (default: -1.0)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Initialize PyRosetta
    pyrosetta.init(
        "-use_input_sc -ignore_unrecognized_res "
        "-ignore_zero_occupancy false -load_PDB_components false -no_fconfig",
        silent=True
    )

    # Load structure and score
    pose = pose_from_pdb(args.input_pdb)
    scorefxn = ScoreFunctionFactory.create_score_function("ref2015")

    print(f"Loaded {args.input_pdb}: {pose.total_residue()} residues")

    binder_energies, target_energies = get_interface_energies(
        pose, scorefxn, args.binder_chain, args.target_chain,
        args.distance_cutoff
    )

    # Write per-binder-residue energies
    energy_csv = os.path.join(args.output_dir, "residue_energy.csv")
    with open(energy_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["residue", "binder_energy"])
        for res_label in sorted(binder_energies.keys(),
                                key=lambda x: int(x[1:]) if x[1:].isdigit() else 0):
            writer.writerow([res_label, f"{binder_energies[res_label]:.4f}"])
    print(f"Wrote binder residue energies to {energy_csv}")
    print(f"  {len(binder_energies)} binder residues at interface")

    # Detect and write hotspots
    hotspots = detect_hotspots(target_energies, args.hotspot_energy_threshold)
    hotspot_file = os.path.join(args.output_dir, "hotspots.txt")
    with open(hotspot_file, "w") as f:
        f.write(",".join(hotspots))
    print(f"Wrote {len(hotspots)} hotspot residues to {hotspot_file}")

    # Print summary
    if binder_energies:
        sorted_binder = sorted(binder_energies.items(), key=lambda x: x[1])
        print("\nTop 10 binder residues by interaction energy:")
        for res, energy in sorted_binder[:10]:
            print(f"  {res}: {energy:.3f} REU")
    if hotspots:
        print(f"\nHotspot residues on target: {','.join(hotspots)}")


if __name__ == "__main__":
    main()
