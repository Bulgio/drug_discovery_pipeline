"""
fix_connectivity.py
Called by simulation.sh as a fallback when antechamber reports bad connectivity.
Reads the ligand PDB, sanitizes it, computes Gasteiger charges, and writes a
cleaned SDF that antechamber can re-read.

Usage: python fix_connectivity.py <ligand_pdb> <out_sdf>
"""
import sys
from Script_finale.GROMACS.fix_rdkit import Chem
from rdkit.Chem import AllChem, SDWriter

pdb_path = sys.argv[1]
sdf_path = sys.argv[2]

mol = Chem.MolFromPDBFile(pdb_path, removeHs=False, proximityBonding=True)
if mol is not None:
    Chem.SanitizeMol(mol)
    AllChem.ComputeGasteigerCharges(mol)
    w = SDWriter(sdf_path)
    w.write(mol)
    w.close()
    print("Created cleaned SDF file")
else:
    print("ERROR: Could not read molecule for connectivity fix")
    sys.exit(1)
