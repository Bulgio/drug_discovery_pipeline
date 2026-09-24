"""
Function to protonate ligands using RDKit package functions

INPUTS:
argv[1] = full path to the raw ligand PDB  (e.g. $WORKDIR/ligand_raw.pdb)
argv[2] = work_dir                         (e.g. $WORKDIR)

Called from simulation.sh as:
    python protonation.py "$WORKDIR/ligand_raw.pdb" "$WORKDIR"

Output:
    $WORKDIR/ligand_protonated.pdb
"""
#!/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem import AllChem
import os
import sys

if len(sys.argv) < 3:
    print("Usage: protonation.py <ligand_raw_pdb> <work_dir>")
    sys.exit(1)

# argv[1] is the FULL PATH to the raw ligand PDB
# argv[2] is the work directory where output is written
# This matches exactly how simulation.sh calls this script:
#   python protonation.py "$WORKDIR/ligand_raw.pdb" "$WORKDIR"
raw_ligand  = sys.argv[1]
work_dir    = sys.argv[2]
prot_ligand = os.path.join(work_dir, 'ligand_protonated.pdb')

def read_mol_heavy_only(path):
    """
    Read a PDB file returning only the heavy-atom skeleton with correct
    bond orders (aromaticity inferred by RDKit).

    Why removeHs=True matters
    -------------------------
    MolFromPDBFile reads CONECT records as single bonds only — the PDB
    format does not encode bond order.  If the molecule is read with
    removeHs=False, the explicit H atoms in the CONECT graph prevent
    RDKit from running its aromaticity/Kekulization algorithm, so every
    ring carbon appears sp3 with 2 implicit H instead of sp2 with 1H.
    When AddHs is then called, it fills the inflated sp3 valences and
    produces ~3.5 H per heavy atom instead of the expected ~1.5.

    Reading with removeHs=True lets RDKit ignore the H atoms in the
    connectivity graph, correctly perceive aromaticity, and assign the
    right number of implicit H before AddHs is called.
    """
    # Primary: proximity bonding (for PDB files without CONECT records)
    mol = Chem.MolFromPDBFile(path, removeHs=True,
                              sanitize=True, proximityBonding=True)
    if mol is not None:
        return mol, "proximity-based"
    
    # Fallback: CONECT-based, strip H so aromaticity is inferred correctly
    mol = Chem.MolFromPDBFile(path, removeHs=True,
                              sanitize=True, proximityBonding=False)
    if mol is not None:
        return mol, "CONECT-based"
    # Last resort: PDB block reader
    with open(path, 'r') as f:
        pdb_block = f.read()
    mol = Chem.MolFromPDBBlock(pdb_block, removeHs=True, sanitize=True)
    if mol is not None:
        return mol, "block-based"

    return None, None

try:
    mol, strategy = read_mol_heavy_only(raw_ligand)

    if mol is None:
        raise ValueError("All RDKit reading strategies failed")

    print(f"Molecule read successfully ({strategy}), "
          f"{mol.GetNumAtoms()} heavy atoms")

    # Ensure at least one conformer exists so AddHs can place H coordinates
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

    # Add hydrogens ONCE at pH 7.4.
    # FIX 1: previous code called AddHs twice, producing ~3.5 H/heavy
    #        instead of the expected ~1.5 and corrupting valences.
    # FIX 2: reading with removeHs=True (above) ensures RDKit perceives
    #        aromaticity correctly before AddHs runs.
    mol = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)

    Chem.MolToPDBFile(mol, prot_ligand)

    n_heavy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() != 1)
    n_h     = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
    print(f"RDKit protonation successful: {mol.GetNumAtoms()} atoms total "
          f"(heavy={n_heavy}, H={n_h})")

except Exception as e:
    print(f"ERROR in protonation: {e}")
    import shutil
    shutil.copyfile(raw_ligand, prot_ligand)
    print("Using raw ligand as a fallback")
