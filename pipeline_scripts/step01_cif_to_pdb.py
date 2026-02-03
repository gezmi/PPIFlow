#!/usr/bin/env python3
"""Step 1: Convert CIF to PDB with chain renaming.

Reads a CIF file, renames chains so the nanobody becomes chain A and the
receptor becomes chain C, reorders so chain A appears first in the PDB
(required by the partial flow script which computes ab_length from the
first chain), and writes the result as a PDB file.

Usage:
    python step01_cif_to_pdb.py \
        --input_cif path/to/complex.cif \
        --output_pdb path/to/output.pdb \
        --nanobody_chain B \
        --receptor_chain A
"""

import argparse
import gemmi


def cif_to_pdb(input_cif: str, output_pdb: str,
               nanobody_chain: str, receptor_chain: str) -> None:
    doc = gemmi.cif.read(input_cif)
    st = gemmi.make_structure_from_block(doc[0])
    st.setup_entities()

    model = st[0]

    nb_chain = None
    rec_chain = None
    for chain in model:
        if chain.name == nanobody_chain:
            nb_chain = chain
        elif chain.name == receptor_chain:
            rec_chain = chain

    if nb_chain is None:
        raise ValueError(f"Nanobody chain '{nanobody_chain}' not found in CIF")
    if rec_chain is None:
        raise ValueError(f"Receptor chain '{receptor_chain}' not found in CIF")

    # Rename: nanobody -> A, receptor -> C
    # Use a temporary name to avoid collision if one is already 'A' or 'C'
    nb_chain.name = "_NB_"
    rec_chain.name = "_REC_"
    nb_chain.name = "A"
    rec_chain.name = "C"

    # Also update subchain (entity) names for consistency
    for chain in model:
        for residue in chain:
            if chain.name == "A":
                residue.subchain = "A"
            elif chain.name == "C":
                residue.subchain = "C"

    # Reorder chains: A (nanobody) first, then C (receptor)
    new_model = gemmi.Model(model.name)
    for chain in model:
        if chain.name == "A":
            new_model.add_chain(chain)
    for chain in model:
        if chain.name == "C":
            new_model.add_chain(chain)
    # Include any other chains unchanged (unlikely but safe)
    for chain in model:
        if chain.name not in ("A", "C"):
            new_model.add_chain(chain)

    st[0] = new_model

    st.write_pdb(output_pdb)

    # Report chain lengths
    for chain in st[0]:
        n_res = sum(1 for r in chain if r.entity_type == gemmi.EntityType.Polymer
                    or r.name in gemmi.expand_protein_one_letter_groups(""))
        polymer_res = [r for r in chain
                       if gemmi.find_tabulated_residue(r.name).is_amino_acid()]
        print(f"  Chain {chain.name}: {len(polymer_res)} amino acid residues")


def main():
    parser = argparse.ArgumentParser(
        description="Convert CIF to PDB with chain renaming for PPIFlow pipeline")
    parser.add_argument("--input_cif", required=True, help="Input CIF file path")
    parser.add_argument("--output_pdb", required=True, help="Output PDB file path")
    parser.add_argument("--nanobody_chain", default="B",
                        help="Chain ID of the nanobody in the CIF (default: B)")
    parser.add_argument("--receptor_chain", default="A",
                        help="Chain ID of the receptor in the CIF (default: A)")
    args = parser.parse_args()

    cif_to_pdb(args.input_cif, args.output_pdb,
               args.nanobody_chain, args.receptor_chain)
    print(f"Written PDB to {args.output_pdb}")


if __name__ == "__main__":
    main()
