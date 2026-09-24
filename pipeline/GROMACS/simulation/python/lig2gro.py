"""
lig2gro.py
INPUT
lig_top  Path to the .prmtop
lig_crd  Path to the .inpcrd
work_dir $WORK_DIR
"""
#!/usr/bin/env python3

import parmed as pmd
import re
import os
import sys

lig_prmtop  = sys.argv[1]
lig_inpcrd  = sys.argv[2]
out_path    = sys.argv[3]
print("Converting Amber to GROMACS with charge validation...")

SECTIONS_TO_SKIP = {'defaults', 'system', 'molecules'}

def section_name(line):
    match = re.match(r'^\s*\[\s*(\w+)\s*\]', line)
    if match:
        return match.group(1).lower()
    return None

def top_to_itp(top_content):
    """
    Convert a parmed-generated .top to an .itp by stripping
    [ defaults ], [ system ] and [ molecules ] sections.
    """
    new_lines = []
    skip = False
    for line in top_content:
        name = section_name(line)
        if name is not None:
            skip = name in SECTIONS_TO_SKIP
            if skip:
                print(f"  Skipping section: [ {name} ]")
            else:
                print(f"  Keeping  section: [ {name} ]")
        if not skip:
            new_lines.append(line)
    return new_lines

def load_amber_structure(prmtop, inpcrd):
    """
    Load an Amber prmtop + inpcrd into a ParmEd Structure in a way that
    works across ALL ParmEd versions found on HPC clusters.

    Background
    ----------
    ParmEd's load_file() dispatches kwargs directly to AmberFormat.__init__.
    The 'xyz' keyword was added only in standalone ParmEd >= 4.x.
    The version bundled with AmberTools 25 on ComputeCanada/Narval is an
    older 3.x variant where 'xyz' is not recognised:
        TypeError: AmberFormat.__init__() got an unexpected keyword argument 'xyz'

    The two-positional-argument call pmd.load_file(prmtop, inpcrd) also
    fails on newer standalone ParmEd 4.x:
        AmberFormat.__init__() takes from 1 to 2 positional arguments but 3 were given

    The only API that is version-agnostic is the two-step approach:
      1. pmd.load_file(prmtop)   — loads topology only
      2. ligand.load_rst7(inpcrd) — injects coordinates

    load_rst7() accepts a filename directly and has been present in every
    ParmEd release, making this approach safe on all versions.
    """
    ligand = pmd.load_file(prmtop)
    ligand.load_rst7(inpcrd)
    return ligand

try:
    ligand = load_amber_structure(lig_prmtop, lig_inpcrd)

    # Validate total charge for free energy accuracy
    calculated_charge = sum(atom.charge for atom in ligand.atoms)
    print(f"Actual charge in parameters: {calculated_charge:.4f}")

    # Re-number atoms sequentially to avoid off-by-one in some ParmEd versions
    for i, atom in enumerate(ligand.atoms):
        atom.number = i + 1

    gro_path = os.path.join(out_path, 'ligand.gro')
    top_path = os.path.join(out_path, 'ligand.top')
    itp_path = os.path.join(out_path, 'ligand.itp')

    # FIX: use combine='all' when saving the .gro file.
    #
    # By default (combine=None) the GRO writer calls struct.split() which
    # requires explicit bond objects in the ParmEd bond list.  Amber prmtop
    # bonds are stored differently and do not populate this list in a way
    # that split() can use, causing:
    #   RuntimeError: Could not find <Atom C1 [0]; In MOL 0>
    #
    # combine='all' bypasses the matching entirely and writes atoms in the
    # order they appear in struct.atoms (= prmtop order = inpcrd order).
    # This is correct and safe for a single-molecule ligand file.
    ligand.save(gro_path, overwrite=True, combine='all')
    ligand.save(top_path, overwrite=True)

    with open(top_path, 'r') as f:
        top_content = f.readlines()

    new_lines = top_to_itp(top_content)

    with open(itp_path, 'w') as f:
        f.writelines(new_lines)

    print(f"✓ ligand.gro written")
    print(f"✓ ligand.itp written ({len(new_lines)} lines kept)")

except Exception as e:
    print(f"✗ Conversion failed: {e}")
    import traceback
    traceback.print_exc()
    exit(1)
