#!/usr/bin/env python3

from Script_finale.GROMACS.fix_rdkit import Chem
import sys
import os

if len(sys.argv) < 3:
    print("ERROR: No input file provided")
    print("Usage: pdb2sdf.py <pdb_file> <out_sdf>")
    sys.exit(1)

pdb_file = sys.argv[1]
out_path = sys.argv[2]

def fix_pdb_columns(path):
    """
    Fix non-standard 4-char residue names (e.g. LIG1) that shift all
    subsequent PDB columns by one, causing RDKit to fail reading resnum.
    """
    lines = open(path).readlines()
    fixed = []
    changed = False
    for line in lines:
        if line.startswith(('HETATM', 'ATOM')):
            if line[20:21] not in (' ', '\n', '\r'):
                resname = line[17:20]
                chain   = line[22:23]
                resnum  = line[23:27]
                line = (line[0:16] + ' ' + resname + ' ' +
                        chain + resnum + ' ' + '   ' + line[30:])
                changed = True
        fixed.append(line)

    if not changed:
        return path

    fixed_path = path.replace('.pdb', '_fixed.pdb')
    with open(fixed_path, 'w') as f:
        f.writelines(fixed)
    print(f"PDB column fix applied (4-char resname → 3-char)")
    return fixed_path

def read_mol(path):
    mol = Chem.MolFromPDBFile(path, removeHs=True, proximityBonding=False)
    if mol is not None:
        return mol, "CONECT-based"

    mol = Chem.MolFromPDBFile(path, removeHs=True, proximityBonding=True)
    if mol is not None:
        return mol, "proximity-based"

    with open(path, 'r') as f:
        pdb_block = f.read()
    mol = Chem.MolFromPDBBlock(pdb_block, removeHs=True, proximityBonding=True)
    if mol is not None:
        return mol, "block+proximity-based"

    return None, None

try:
    pdb_path = fix_pdb_columns(pdb_file)
    mol, strategy = read_mol(pdb_path)

    if mol is not None:
        print(f"Molecule read successfully ({strategy})")
        writer = Chem.SDWriter(out_path)
        writer.write(mol)
        writer.close()
        print(f"Successfully converted {pdb_file} to {out_path}")
        print(f"Molecule has {mol.GetNumAtoms()} atoms")
        sys.exit(0)
    else:
        raise ValueError("All RDKit reading strategies failed")

except Exception as e:
    print(f"ERROR in sdf creation: {e}")
    sys.exit(1)
