"""
Docstring for GROMACS.job_scripts.top_writer

Inputs: 
ff_path     :       str     "$WORKDIR/forcefield.itp"
ligand_path :       str     "$WORKDIR/ligand.itp"
water_path  :       str     "$WORKDIR/spce.itp"
protein_path:       str     "$WORKDIR/protein.itp"

ion_path    :       str     (may not be used if is already present)
water_custom:       str     (optional)
out_path    :       str

Objective: 
#include "amber99sb-ildn.ff/forcefield.itp"
#include "ligand.itp"      ; Must come before protein
#include "protein.itp"     ; Remove duplicate forcefield/water includes!
#include "amber99sb-ildn.ff/spce.itp"
#include "amber99sb-ildn.ff/ions.itp"

[ system ]
Protein-Ligand Complex

[ molecules ]
Protein_chain_A   1   ; Coordinates first in gro
Ligand            1   ; Coordinates second in gro
SOL               X   ; Coordinates third in gro
NA/CL             Y   ; Coordinates last in gro
"""
import sys
import os

#Functions

def file_reader(path :str ):

    with open(path,'r') as f:
        lines = f.readlines()
    return lines

def line_include(path):
    # Use only the filename — grompp resolves includes relative to its working
    # directory ($WORKDIR), where all .itp files have already been copied.
    line = f'#include "{os.path.basename(path)}"\n'

    return line

def top2itp(prot_path,dir_path):

    top_lines=file_reader(prot_path)
    new_lines=[]
    capture=False
    name='protein.itp'
    out_path=os.path.join(dir_path,name)
    
    for line in top_lines:
        if line.strip() == '[ moleculetype ]': #Eliminates everything until it encounters this
            capture=True
        elif line.strip() == '; Include water topology':
            capture = False
            break
        
        if capture:
            new_lines.append(line)
        else:
            continue
    
    if out_path:
        with open(out_path, 'w') as f:
            f.writelines(new_lines)
        print(f"✓ protein.itp")
        return out_path

def get_molecule_line(itp_path):
    """Extract molecule name from [ moleculetype ] line"""
    with open(itp_path, 'r') as f:
        for line in f:
            if line.strip().startswith('[ moleculetype ]'):
                # Next non-empty, non-comment line has the name
                for next_line in f:
                    next_line = next_line.strip()
                    if next_line and not next_line.startswith(';'):
                        return next_line.split()[0]  # First word is the name
    return None  # Fallback

#Main
ion_check = False
water_check = False

if len(sys.argv) < 6:
    print("ERROR: Need at least: ff_path ligand_path water_path protein_path out_path")
    sys.exit(1)
elif len(sys.argv) == 6:
    ion_check = False      # No ion file
    water_check = False    # Default water
elif len(sys.argv) == 7:
    ion_check = True       # Ion file provided
    water_check = False    # Use default water
elif len(sys.argv) == 8:
    ion_check = True       # Ion file provided
    water_check = True     # Use water customized
else:
    print("ERROR: Too many arguments (max 8)")
    sys.exit(1)

ff_path     =sys.argv[1]#:       str     "$WORKDIR/forcefield.itp"
ligand_path =sys.argv[2]#:       str     "$WORKDIR/ligand.itp"
water_path  =sys.argv[3]#:       str     "$WORKDIR/spce.itp"
protein_path=sys.argv[4]#:       str     "$WORKDIR/protein.top"
out_path    =sys.argv[5]#:       str     "$WORKDIR"
if ion_check:
    ion_path    =sys.argv[6]#:       str     (may not be used if is already present)
if water_check:
    water_custom=sys.argv[7]#:       str     (optional)

dir_path=os.path.dirname(ff_path) #$WORKDIR
prot_itp_path=top2itp(protein_path,dir_path) #"$WORKDIR/protein.itp"
top_path=os.path.join(out_path,'complex.top')

with open(top_path,'w') as file:
    file.write(';HEADER\n\n;Forcefield : amber99sb-ildn\n')
    file.write(line_include(ff_path))
    file.write(line_include(ligand_path))
    file.write(line_include(prot_itp_path))
    if not water_check:
        file.write(line_include(water_path))
    else:
        file.write(line_include(water_custom))
    if ion_check:
        file.write(line_include(ion_path))
    
    # ===== SYSTEM DEFINITION - CRITICAL! =====
    file.write('\n[ system ]\n')
    file.write('; Name\n')
    file.write('Protein-Ligand Complex in Water\n\n')
    
    # ===== MOLECULES - CRITICAL! =====
    file.write('[ molecules ]\n')
    file.write('; Compound      #mols\n')
    file.write(f'{get_molecule_line(prot_itp_path)}         1\n')
    file.write(f'{get_molecule_line(ligand_path)}           1\n')
    #sol modified by solvate
    #ions added by genion

print(f"✓ Master topology written to {top_path}")
print(f"\n{'='*60}")