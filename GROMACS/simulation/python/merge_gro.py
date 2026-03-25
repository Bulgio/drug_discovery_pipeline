#!/usr/bin/env python3
"""
merge_gro.py
Merges protein_processed.gro and ligand.gro preserving all original coordinates.
Replaces gmx insert-molecules which randomizes ligand placement.

Usage: python merge_gro.py <protein.gro> <ligand.gro> <output.gro>
"""

import sys

def read_gro(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    title   = lines[0]
    n_atoms = int(lines[1].strip())
    atoms   = lines[2:2 + n_atoms]
    box     = lines[2 + n_atoms]
    return title, n_atoms, atoms, box

def write_gro(path, title, atoms, box):
    with open(path, 'w') as f:
        f.write(title.strip() + '\n')
        f.write(f'{len(atoms)}\n')
        for atom in atoms:
            f.write(atom)
        f.write(box)

if len(sys.argv) != 4:
    print("Usage: merge_gro.py <protein.gro> <ligand.gro> <output.gro>")
    sys.exit(1)

prot_path = sys.argv[1]
lig_path  = sys.argv[2]
out_path  = sys.argv[3]

prot_title, prot_n, prot_atoms, prot_box = read_gro(prot_path)
lig_title,  lig_n,  lig_atoms,  lig_box  = read_gro(lig_path)

# Renumber atoms sequentially — residue numbers preserved from each file
merged_atoms = []
atom_counter = 1
for line in prot_atoms + lig_atoms:
    # GRO format: cols 0-4 resnum, 5-9 resname, 10-14 atomname, 15-19 atomnum, then coords
    new_line = line[:15] + f'{atom_counter:5d}' + line[20:]
    merged_atoms.append(new_line)
    atom_counter += 1

# Use protein box dimensions (the ligand box from editconf is meaningless here)
title = "Protein-Ligand Complex"
write_gro(out_path, title, merged_atoms, prot_box)

print(f"✓ Merged {prot_n} protein atoms + {lig_n} ligand atoms = {len(merged_atoms)} total")
print(f"✓ Box dimensions from protein: {prot_box.strip()}")
print(f"✓ Written to {out_path}")
