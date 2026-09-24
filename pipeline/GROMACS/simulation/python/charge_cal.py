"""
charge_calc_complete
"""
#!/usr/bin/env python3

from Script_finale.GROMACS.fix_rdkit import Chem
import sys
import os


def elements(lig_path):
    """
    Fallback charge estimator based on element recognition.
    Only handles isolated metal cations reliably.
    For organic molecules RDKit (rdkit_cal) should always be preferred.
    """
    positive = 0

    with open(lig_path, 'r') as f:
        for line in f:
            if line.startswith(('HETATM', 'ATOM')):
                element = line[76:78].strip()

                # Common divalent cations
                if element in ['MG', 'CA', 'ZN', 'FE', 'MN', 'CU']:
                    positive += 2
                # Common monovalent cations
                elif element in ['NA', 'K']:
                    positive += 1

    # Anion detection removed: atom-name heuristics (e.g. 'O' in name) are
    # unreliable for organic ligands and produce large negative errors.
    # If no metals are found, 0 is the safest fallback for organic molecules.
    if positive == 0:
        print("WARNING: No metal ions detected — assuming organic molecule, ELEMENT_CHARGE=0")

    net_charge = positive
    # FIX: removed trailing \n — bash `cut -d'=' -f2` would capture an empty
    # second line, causing FORMAL_CHARGE / ELEMENT_CHARGE to be read as "".
    print(f"ELEMENT_CHARGE={net_charge}")
    return net_charge

def rdkit_cal(protonated_path):
    if not os.path.exists(protonated_path):
        print(f"ERROR: File {protonated_path} not found")
        sys.exit(1)

    try:
        mol = Chem.MolFromPDBFile(protonated_path, removeHs=False)
        if mol is None:
            # Try alternative reading
            with open(protonated_path, 'r') as f:
                pdb_block = f.read()
            mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False)

        if mol is None:
            print("ERROR: Failed to read molecule from PDB")
            sys.exit(1)

        formal_charge = Chem.GetFormalCharge(mol)
        # FIX: removed trailing \n (same reason as ELEMENT_CHARGE above)
        print(f"FORMAL_CHARGE={formal_charge}")

        # Additional diagnostics
        smi = Chem.MolToSmiles(mol)
        if '[NH3+]' in smi or '[NH4+]' in smi or '[nH+]' in smi:
            print("POSITIVE_GROUP_DETECTED")
        elif '[O-]' in smi or '[C-]' in smi or '[N-]' in smi or '[S-]' in smi:
            print("NEGATIVE_GROUP_DETECTED")
        return formal_charge
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)

"""
Main code
"""

if len(sys.argv) < 2:
    print("ERROR: No input file provided")
    print("Usage: charge_calc.py <pdb_file>")
    sys.exit(1)

lig_path = sys.argv[1]

el_charge = elements(lig_path)
formal_charge = rdkit_cal(lig_path)
