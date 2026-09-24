#!/bin/bash
#SBATCH --job-name=simulation
#SBATCH --account=def-jtus
#SBATCH --time=36:00:00
#SBATCH --gpus=a100_3g.20gb:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/home/toscana/scratch/new_GROMACS/logs/sim/output_%j.txt
#SBATCH --error=/home/toscana/scratch/new_GROMACS/logs/sim/error_%j.txt
#SBATCH --mail-user=matteo.bulgini@studenti.polito.it
#SBATCH --mail-type=ALL

# Module loading
module purge
module load StdEnv/2023
module load gcc/12.3
module load python/3.12.4
module load cuda/12.6
module load openmpi/4.1.5
module load gromacs/2025.4
module load rdkit/2024.09.6
module load ambertools/25.0
module load openbabel/3.1.1

check_file() {
    if [ ! -s "$1" ]; then
        echo "❌ ERROR: $1 is empty or missing!"
        return 1
    fi
    echo "✓ $1 looks good"
    return 0
}

PDB_FILE="$1"
OUTPUT_DIR="$2"
FILENAME="$3"
WORKDIR="$OUTPUT_DIR/temp"
JOB_DIR="$4"
INP_DIR="$5"
RESULTS_DIR="$OUTPUT_DIR/results"

echo "=== Starting GROMACS simulation ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "PDB File: $PDB_FILE"
echo "Output Directory: $OUTPUT_DIR"

cd "$OUTPUT_DIR"

mkdir -p $WORKDIR
mkdir -p $RESULTS_DIR

echo "=== EXTRACTING PROTEIN AND LIGAND ==="

if [ ! -f "$PDB_FILE" ]; then
    echo "ERROR: PDB file $PDB_FILE not found!"
    exit 1
fi

# Strip CRLF from input PDB (files from docking tools often have Windows line endings)
sed 's/\r//' "$PDB_FILE" > "$WORKDIR/input_clean.pdb"

# Extract protein (ATOM records + crystal waters)
grep -E '(^ATOM  |^HETATM.*HOH)' "$WORKDIR/input_clean.pdb" > "$WORKDIR/protein.pdb" 2>/dev/null
# Extract ligand (non-water HETATM only)
grep '^HETATM' "$WORKDIR/input_clean.pdb" | grep -v 'HOH' > "$WORKDIR/ligand_raw.pdb" 2>/dev/null
echo 'END' >> "$WORKDIR/ligand_raw.pdb"

if [ ! -s "$WORKDIR/ligand_raw.pdb" ]; then
    echo "❌ FATAL ERROR: No ligand found in PDB file. This pipeline requires a ligand."
    exit 1
fi

echo "=== PROTONATION OF THE LIGAND ==="

echo "=== PROCESSING ==="
echo "Ligand atoms: $(grep -c "HETATM" "$WORKDIR/ligand_raw.pdb")"
echo "Ligand hydrogens: $(grep " H " "$WORKDIR/ligand_raw.pdb" | wc -l)"

cp "$WORKDIR/ligand_raw.pdb" "$WORKDIR/ligand_original.pdb"
LIGAND_H_COUNT=$(grep " H " "$WORKDIR/ligand_raw.pdb" | wc -l)
if [ "$LIGAND_H_COUNT" -eq 0 ]; then
    echo "Ligand has no hydrogens. Adding hydrogens at pH 7.4..."
    python "$JOB_DIR/fix_rdkit.py" "$WORKDIR/ligand_raw.pdb" "$WORKDIR/ligand_fixed.pdb"
    python "$JOB_DIR/protonation.py" "$WORKDIR/ligand_fixed.pdb" "$WORKDIR"
else
    echo "Ligand already has hydrogens. Copying as-is."
    cp "$WORKDIR/ligand_raw.pdb" "$WORKDIR/ligand_protonated.pdb"
fi

cp "$WORKDIR/ligand_protonated.pdb" "$WORKDIR/ligand.pdb"

echo "=== VALIDATING PROTONATION ==="
FINAL_ATOMS=$(grep -c "^HETATM" "$WORKDIR/ligand.pdb")
HEAVY_ATOMS=$(grep "^HETATM" "$WORKDIR/ligand.pdb" | grep -v " H " | wc -l)
HYDROGENS=$(grep "^HETATM" "$WORKDIR/ligand.pdb" | grep " H " | wc -l)
echo "Final ligand: Total=$FINAL_ATOMS Heavy=$HEAVY_ATOMS Hydrogens=$HYDROGENS"

echo "=== STEP 2: CHARGE DETERMINATION ==="
python "$JOB_DIR/charge_cal.py" "$WORKDIR/ligand_protonated.pdb" > "$WORKDIR/charge.txt"

FORMAL_CHARGE=$(grep "FORMAL_CHARGE=" "$WORKDIR/charge.txt" | cut -d'=' -f2)
ELEMENT_CHARGE=$(grep "ELEMENT_CHARGE=" "$WORKDIR/charge.txt" | cut -d'=' -f2)

echo "Formal charge (RDKit): |$FORMAL_CHARGE|"
echo "Element heuristic:     |$ELEMENT_CHARGE|"

if [ -n "$FORMAL_CHARGE" ] && [[ "$FORMAL_CHARGE" =~ ^-?[0-9]+$ ]]; then
    CHARGE="$FORMAL_CHARGE"
    CHARGE_METHOD="formal"
elif [ -n "$ELEMENT_CHARGE" ] && [[ "$ELEMENT_CHARGE" =~ ^-?[0-9]+$ ]]; then
    CHARGE="$ELEMENT_CHARGE"
    CHARGE_METHOD="element"
else
    CHARGE="0"
    CHARGE_METHOD="assumed_neutral"
    echo "WARNING: Could not determine charge, assuming 0"
fi
echo "Final charge: $CHARGE (method: $CHARGE_METHOD)"

echo "=== STEP 3: BCC PARAMETERIZATION ==="

python "$JOB_DIR/pdb2sdf.py" "$WORKDIR/ligand_protonated.pdb" "$WORKDIR/ligand_py.sdf"

MAX_RETRIES=3
RETRY_COUNT=0
BCC_SUCCESS=false

while [ "$RETRY_COUNT" -le "$MAX_RETRIES" ] && [ "$BCC_SUCCESS" = false ]; do
    echo "BCC attempt $((RETRY_COUNT + 1)) of $((MAX_RETRIES + 1))"
    rm -f "$WORKDIR/ligand.mol2" "$WORKDIR/sqm.out" "$WORKDIR/sqm.in"

    if antechamber \
        -i "$WORKDIR/ligand_py.sdf" \
        -fi sdf \
        -o "$WORKDIR/ligand.mol2" \
        -fo mol2 \
        -c bcc \
        -nc "$CHARGE" \
        -at gaff2 \
        -j 4 \
        -s 2 \
        2> "$WORKDIR/antechamber_attempt_${RETRY_COUNT}.log"; then

        BCC_SUCCESS=true
        echo "✓ BCC successful!"
    else
        RETRY_COUNT=$((RETRY_COUNT + 1))
        if [ "$RETRY_COUNT" -le "$MAX_RETRIES" ]; then
            echo "AM1-BCC failed, retrying..."
            if grep -q "bad connectivity\|Cannot assign bond" "$WORKDIR/antechamber_attempt_$((RETRY_COUNT-1)).log"; then
                echo "Fixing connectivity issues..."
                python "$JOB_DIR/fix_connectivity.py" "$WORKDIR/ligand.pdb" "$WORKDIR/ligand_clean.sdf"
                if [ -f "$WORKDIR/ligand_clean.sdf" ]; then
                    obabel "$WORKDIR/ligand_clean.sdf" -O "$WORKDIR/ligand.pdb" 2>/dev/null
                fi
            fi
            obabel "$WORKDIR/ligand.pdb" -O "$WORKDIR/ligand.pdb" --minimize --steps 1000 2>/dev/null
        else
            echo "BCC failed after $MAX_RETRIES attempts. Falling back to Gasteiger charges..."
            antechamber \
                -i "$WORKDIR/ligand.pdb" \
                -fi pdb \
                -o "$WORKDIR/ligand.mol2" \
                -fo mol2 \
                -c gas \
                -nc "$CHARGE" \
                2>"$WORKDIR/antechamber_gas.log"
            if [ $? -eq 0 ]; then
                echo "Using Gasteiger charges (accuracy reduced)"
                BCC_SUCCESS=true
            else
                echo "FATAL: All parameterization methods failed"
                exit 1
            fi
        fi
    fi
done

    echo "=== STEP 4: GENERATING FF PARAMETERS ==="
    if ! parmchk2 -i "$WORKDIR/ligand.mol2" -f mol2 -o "$WORKDIR/ligand.frcmod" -s gaff2 2>"$WORKDIR/parmchk.log"; then
        echo "WARNING: parmchk2 had issues. Creating minimal frcmod..."
        echo "Remark line created by tleap" > "$WORKDIR/ligand.frcmod"
        echo "MASS" >> "$WORKDIR/ligand.frcmod"
    fi

    if grep -q "ATTN: need revision" "$WORKDIR/ligand.frcmod"; then
        MISSING_PARAMS=$(grep -c "ATTN" "$WORKDIR/ligand.frcmod")
        echo "⚠️  Warning: $MISSING_PARAMS missing parameters detected"
    fi

    cat > "$WORKDIR/leap.in" << EOF
source leaprc.gaff2
source leaprc.water.spce

LIG = loadmol2 $WORKDIR/ligand.mol2
loadamberparams $WORKDIR/ligand.frcmod

check LIG
saveamberparm LIG $WORKDIR/ligand.prmtop $WORKDIR/ligand.inpcrd

quit
EOF

    if ! tleap -f "$WORKDIR/leap.in" > "$WORKDIR/tleap.log" 2>&1; then
        echo "ERROR: tleap failed! Check $WORKDIR/tleap.log"
        exit 1
    fi

    echo "=== STEP 5: LIGAND GROMACS CONVERSION ==="
    # FIX: lig2gro.py converts from Amber to GROMACS preserving crystal coordinates.
    # We do NOT run gmx editconf on the ligand — that would destroy the binding pose.
    python "$JOB_DIR/lig2gro.py" "$WORKDIR/ligand.prmtop" "$WORKDIR/ligand.inpcrd" "$WORKDIR"

    if [ ! -f "$WORKDIR/ligand.gro" ]; then
        echo "❌ ERROR: ligand.gro not created!"
        exit 1
    fi
    echo "✓ ligand.gro created with crystal coordinates preserved"

    # FIX: Generate position restraints for the ligand heavy atoms.
    # These are activated during NVT and NPT equilibration via -DPOSRES_LIG.
    echo "=== GENERATING LIGAND POSITION RESTRAINTS ==="
    echo "0" | gmx genrestr \
        -f "$WORKDIR/ligand.gro" \
        -o "$WORKDIR/posre_ligand.itp" \
        -fc 1000 1000 1000 \
        -quiet 2>"$WORKDIR/genrestr.log"

    if [ ! -f "$WORKDIR/posre_ligand.itp" ]; then
        echo "❌ ERROR: posre_ligand.itp not created!"
        exit 1
    fi
    echo "✓ posre_ligand.itp created"

    # FIX: Append the posre_ligand.itp include block to ligand.itp so it activates
    # automatically when -DPOSRES_LIG is defined in nvt.mdp / npt.mdp.
    cat >> "$WORKDIR/ligand.itp" << 'EOF'

; Position restraints for ligand — activated by define = -DPOSRES_LIG
#ifdef POSRES_LIG
#include "posre_ligand.itp"
#endif
EOF
    echo "✓ posre_ligand.itp linked into ligand.itp"

# PROTEIN PROCESS
echo "=== PROCESSING PROTEIN ==="

if ! gmx pdb2gmx \
    -f "$WORKDIR/protein.pdb" \
    -o "$WORKDIR/protein_processed.gro" \
    -p "$WORKDIR/protein_processed.top" \
    -water spce \
    -ignh \
    -i "$WORKDIR/posre.itp" \
    -ff amber99sb-ildn \
    -quiet 2>"$WORKDIR/pdb2gmx.log"; then
    echo "ERROR: pdb2gmx failed! Check $WORKDIR/pdb2gmx.log"
    exit 1
fi

echo "=== PROTEIN PROCESSED ==="

PROTEIN_NAME=$(grep -m1 '\[ molecules \]' "$WORKDIR/protein_processed.top" -A 2 | tail -n 1 | awk '{print $1}')
echo "Protein molecule name: $PROTEIN_NAME"

echo "=== CREATING PROTEIN-LIGAND COMPLEX ==="

# FIX: Use merge_gro.py instead of gmx insert-molecules.
# This preserves the crystal binding pose for both protein and ligand.
python "$JOB_DIR/merge_gro.py" \
    "$WORKDIR/protein_processed.gro" \
    "$WORKDIR/ligand.gro" \
    "$WORKDIR/complex_raw.gro"

if [ ! -s "$WORKDIR/complex_raw.gro" ]; then
    echo "❌ ERROR: merge_gro.py failed to create complex_raw.gro"
    exit 1
fi
echo "✓ Complex assembled with crystal coordinates"

echo "=== CREATING COMBINED TOPOLOGY ==="
FF_PATH="$(dirname $(which gmx))/../share/gromacs/top/amber99sb-ildn.ff/forcefield.itp"
FF_DIR="$(dirname $(which gmx))/../share/gromacs/top/amber99sb-ildn.ff"
WATER_PATH="$(dirname $(which gmx))/../share/gromacs/top/amber99sb-ildn.ff/spce.itp"
ION_PATH="$(dirname $(which gmx))/../share/gromacs/top/amber99sb-ildn.ff/ions.itp"

cp "$FF_PATH"              "$WORKDIR/forcefield.itp"
cp "$FF_DIR/ffnonbonded.itp" "$WORKDIR/ffnonbonded.itp"
cp "$FF_DIR/ffbonded.itp"  "$WORKDIR/ffbonded.itp"
cp "$WATER_PATH"           "$WORKDIR/spce.itp"
cp "$ION_PATH"             "$WORKDIR/ions.itp"
FF_PATH="$WORKDIR/forcefield.itp"
WATER_PATH="$WORKDIR/spce.itp"
ION_PATH="$WORKDIR/ions.itp"

python "$JOB_DIR/top_writer.py" \
    "$FF_PATH" \
    "$WORKDIR/ligand.itp" \
    "$WATER_PATH" \
    "$WORKDIR/protein_processed.top" \
    "$WORKDIR" \
    "$ION_PATH"

if [ ! -f "$WORKDIR/complex.top" ]; then
    echo "ERROR: top_writer.py failed to create complex.top"
    exit 1
fi

CURRENT_GRO="$WORKDIR/complex_raw.gro"
CURRENT_TOP="$WORKDIR/complex.top"

echo "=== COMPLEX PREPARATION ENDED ==="

check_file "$CURRENT_GRO"
check_file "$CURRENT_TOP"
check_file "$INP_DIR/ion.mdp"

echo "=== BOX CREATION ==="
# Center the complex and define solvent box (1.2 nm padding, cubic)
gmx editconf \
    -f "$CURRENT_GRO" \
    -o "$WORKDIR/complex_centered.gro" \
    -c -d 1.2 -bt cubic \
    -quiet 2>"$WORKDIR/editconf.log"

BOX_INFO=$(tail -n 1 "$WORKDIR/complex_centered.gro")
echo "Box dimensions: $BOX_INFO"

echo "=== SOLVATE THE BOX ==="
gmx solvate \
    -cp "$WORKDIR/complex_centered.gro" \
    -cs spc216.gro \
    -o "$WORKDIR/complex_solvated.gro" \
    -p "$CURRENT_TOP" \
    -quiet 2>"$WORKDIR/solvate.log"

WATER_COUNT=$(grep -c "SOL" "$WORKDIR/complex_solvated.gro" || echo "0")
echo "Water molecules added: $WATER_COUNT"

echo "=== RUNNING GROMPP FOR IONS ==="
gmx grompp \
    -f "$INP_DIR/ion.mdp" \
    -c "$WORKDIR/complex_solvated.gro" \
    -p "$CURRENT_TOP" \
    -o "$WORKDIR/ions.tpr" \
    -pp \
    -maxwarn 5 \
    -quiet 2>"$WORKDIR/grompp_ions.log"

check_file "$WORKDIR/complex_solvated.gro"
if ! check_file "$WORKDIR/ions.tpr"; then
    cat "$WORKDIR/grompp_ions.log"
    exit 1
fi

echo "=== ADDING IONS ==="
echo "SOL" | gmx genion \
    -s "$WORKDIR/ions.tpr" \
    -o "$WORKDIR/complex_final.gro" \
    -p "$CURRENT_TOP" \
    -pname NA -nname CL \
    -neutral -conc 0.15 \
    -quiet 2>"$WORKDIR/genion.log"

if [ ! -s "$WORKDIR/complex_final.gro" ]; then
    echo "❌ Ion placement failed!"
    tail -20 "$WORKDIR/genion.log"
    exit 1
fi
echo "✓ Ions added successfully"

check_file "$CURRENT_TOP"
check_file "$WORKDIR/ligand.itp"

FINAL_ATOMS=$(tail -n 1 "$WORKDIR/complex_final.gro" | awk '{print $1}')
echo "Total atoms in final system: $FINAL_ATOMS"

TOP_MOLECULES=$(grep -A 100 "\[ molecules \]" "$CURRENT_TOP" | grep -v "\[" | grep -v "^;" | grep -v "^$")
echo "Topology molecules: $TOP_MOLECULES"

echo "=== START MD SIMULATION ==="

# STEP 1: ENERGY MINIMIZATION
echo "=== ENERGY MINIMIZATION ==="
gmx grompp \
    -f "$INP_DIR/min.mdp" \
    -c "$WORKDIR/complex_final.gro" \
    -p "$CURRENT_TOP" \
    -o "$WORKDIR/em.tpr" \
    -maxwarn 5 \
    -quiet 2>"$WORKDIR/grompp_em.log"

# FIX: Replace -deffnm (deprecated in 2025) with explicit file flags
gmx mdrun \
    -s "$WORKDIR/em.tpr" \
    -o "$WORKDIR/em.trr" \
    -c "$WORKDIR/em.gro" \
    -e "$WORKDIR/em.edr" \
    -g "$WORKDIR/em.log" \
    -ntmpi 1 -ntomp 8

if [ ! -s "$WORKDIR/em.gro" ]; then
    echo "❌ ERROR: Energy minimization failed!"
    tail -20 "$WORKDIR/em.log"
    exit 1
fi
echo "✓ Energy minimization complete"

# STEP 2: NVT EQUILIBRATION
# FIX: nvt.mdp must have: define = -DPOSRES -DPOSRES_LIG
# This restrains protein heavy atoms (posre.itp) and ligand heavy atoms (posre_ligand.itp)
echo "=== NVT EQUILIBRATION ==="
gmx grompp \
    -f "$INP_DIR/nvt.mdp" \
    -c "$WORKDIR/em.gro" \
    -p "$CURRENT_TOP" \
    -r "$WORKDIR/complex_final.gro" \
    -o "$WORKDIR/nvt.tpr" \
    -maxwarn 5 \
    -quiet 2>"$WORKDIR/grompp_nvt.log"

gmx mdrun \
    -s "$WORKDIR/nvt.tpr" \
    -o "$WORKDIR/nvt.trr" \
    -c "$WORKDIR/nvt.gro" \
    -e "$WORKDIR/nvt.edr" \
    -g "$WORKDIR/nvt.log" \
    -cpo "$WORKDIR/nvt.cpt" \
    -ntmpi 1 -ntomp 8 -nb gpu

if [ ! -s "$WORKDIR/nvt.gro" ]; then
    echo "❌ ERROR: NVT equilibration failed!"
    tail -20 "$WORKDIR/nvt.log"
    exit 1
fi
echo "✓ NVT equilibration complete"

# STEP 3: NPT EQUILIBRATION
# FIX: -r reference added; nvt.mdp must have: define = -DPOSRES -DPOSRES_LIG
echo "=== NPT EQUILIBRATION ==="
gmx grompp \
    -f "$INP_DIR/npt.mdp" \
    -c "$WORKDIR/nvt.gro" \
    -p "$CURRENT_TOP" \
    -r "$WORKDIR/complex_final.gro" \
    -t "$WORKDIR/nvt.cpt" \
    -o "$WORKDIR/npt.tpr" \
    -maxwarn 5 \
    -quiet 2>"$WORKDIR/grompp_npt.log"

gmx mdrun \
    -s "$WORKDIR/npt.tpr" \
    -o "$WORKDIR/npt.trr" \
    -c "$WORKDIR/npt.gro" \
    -e "$WORKDIR/npt.edr" \
    -g "$WORKDIR/npt.log" \
    -cpo "$WORKDIR/npt.cpt" \
    -ntmpi 1 -ntomp 8 -nb gpu

if [ ! -s "$WORKDIR/npt.gro" ]; then
    echo "❌ ERROR: NPT equilibration failed!"
    tail -20 "$WORKDIR/npt.log"
    exit 1
fi
echo "✓ NPT equilibration complete"

# STEP 4: PRODUCTION MD
echo "=== PRODUCTION MD ==="
gmx grompp \
    -f "$INP_DIR/md.mdp" \
    -c "$WORKDIR/npt.gro" \
    -p "$CURRENT_TOP" \
    -t "$WORKDIR/npt.cpt" \
    -o "$WORKDIR/md.tpr" \
    -maxwarn 5 \
    -quiet 2>"$WORKDIR/grompp_md.log"

gmx mdrun \
    -s "$WORKDIR/md.tpr" \
    -x "$WORKDIR/md.xtc" \
    -c "$WORKDIR/md.gro" \
    -e "$WORKDIR/md.edr" \
    -g "$WORKDIR/md.log" \
    -cpo "$WORKDIR/md.cpt" \
    -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu

if [ ! -s "$WORKDIR/md.xtc" ]; then
    echo "❌ ERROR: Production MD failed!"
    tail -20 "$WORKDIR/md.log"
    exit 1
fi
echo "✓ Production MD complete"

echo "=== SIMULATION COMPLETE ==="

R_LOG_DIR="$RESULTS_DIR/logs"
mkdir -p $R_LOG_DIR

echo "=== POST-SIMULATION ANALYSIS ==="

# FIX: Copy essential files to results BEFORE analysis
# so all analysis references consistent paths
cp "$WORKDIR/md.gro"          "$RESULTS_DIR/final_structure.gro"
cp "$WORKDIR/md.xtc"          "$RESULTS_DIR/trajectory.xtc"
cp "$WORKDIR/md.xtc"          "$RESULTS_DIR/md_prod.xtc"
cp "$WORKDIR/md.edr"          "$RESULTS_DIR/energies.edr"
cp "$WORKDIR/md.log"          "$RESULTS_DIR/simulation.log"
cp "$WORKDIR/md.tpr"          "$RESULTS_DIR/md.tpr"
cp "$WORKDIR/complex.top"     "$RESULTS_DIR/system.topology"
cp "$WORKDIR/complex_final.gro" "$RESULTS_DIR/initial_system.gro"

# Create master index file
echo "=== CREATING MASTER INDEX FILE ==="

# Step 1: create the MOL residue group without renaming yet,
# so we can read the group number that make_ndx assigns dynamically.
echo -e "keep 1\nr MOL\nq" | gmx make_ndx \
    -f "$WORKDIR/complex_final.gro" \
    -o "$RESULTS_DIR/index.ndx" \
    2>"$R_LOG_DIR/make_ndx_master.log"

# Step 2: find the group number assigned to MOL by parsing the log.
# make_ndx prints lines like: "  13 MOL          :   42 atoms"
# Note: ParmEd names the ligand residue 'MOL' (Amber convention), not 'LIG'
LIG_GROUP=$(grep -oP '(?<=^\s{0,4})\d+(?=\s+MOL)' "$R_LOG_DIR/make_ndx_master.log" | tail -1)
if [ -z "$LIG_GROUP" ]; then
    echo "❌ ERROR: Could not determine LIG group number from make_ndx output."
    cat "$R_LOG_DIR/make_ndx_master.log"
    exit 1
fi
echo "✓ Ligand index group number: $LIG_GROUP"

# Step 3: rename that group to LIG in the index file.
echo -e "name ${LIG_GROUP} LIG\nq" | gmx make_ndx \
    -f "$WORKDIR/complex_final.gro" \
    -n "$RESULTS_DIR/index.ndx" \
    -o "$RESULTS_DIR/index.ndx" \
    -quiet 2>>"$R_LOG_DIR/make_ndx_master.log"

# Copy index to workdir for gmx_MMPBSA compatibility
cp "$RESULTS_DIR/index.ndx" "$WORKDIR/index.ndx"

# 1. PROTEIN BACKBONE RMSD
# FIX: Use md.tpr as reference (the actual production start, not the minimized structure)
echo "=== PROTEIN BACKBONE RMSD ==="
echo -e "4\n4" | gmx rms \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -o "$RESULTS_DIR/rmsd_protein.xvg" \
    -tu ns \
    -quiet 2>"$R_LOG_DIR/rmsd_protein.log"

# 2. LIGAND RMSD
echo "=== LIGAND RMSD ==="
# Create ligand index — written to RESULTS_DIR
echo -e "r MOL\nq" | gmx make_ndx \
    -f "$WORKDIR/complex_final.gro" \
    -o "$RESULTS_DIR/ligand.ndx" \
    -quiet 2>"$R_LOG_DIR/make_ndx_ligand.log"

# FIX: Use md.tpr as reference; use $RESULTS_DIR/ligand.ndx (consistent path)
echo "MOL" | gmx rms \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -n "$RESULTS_DIR/ligand.ndx" \
    -o "$RESULTS_DIR/rmsd_ligand.xvg" \
    -tu ns \
    -quiet 2>"$R_LOG_DIR/rmsd_ligand.log"

# 3. RMSF
echo "=== RESIDUE RMSF ==="
echo "1" | gmx rmsf \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -o "$RESULTS_DIR/rmsf_residues.xvg" \
    -res \
    -quiet 2>"$R_LOG_DIR/rmsf.log"

# 4. HYDROGEN BONDS (protein-ligand)
# Use new gmx hbond (GROMACS 2024+) with text selections — no index groups needed.
echo "=== HYDROGEN BONDS ==="
gmx hbond \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -ref "Protein" \
    -sel "resname MOL" \
    -num "$RESULTS_DIR/hbonds_num.xvg" \
    -quiet 2>"$R_LOG_DIR/hbond_calc.log"

# 5. ENERGY COMPONENTS
echo "=== ENERGY COMPONENTS ==="
echo "11 12 13 14 15 16 17 18 19 20 21 22" | gmx energy \
    -f "$WORKDIR/md.edr" \
    -o "$RESULTS_DIR/energy_components.xvg" \
    -quiet 2>"$R_LOG_DIR/energy.log"

# 6. SASA
echo "=== SASA ==="
echo "1" | gmx sasa \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -o "$RESULTS_DIR/sasa_protein.xvg" \
    -or "$RESULTS_DIR/sasa_residue.xvg" \
    -quiet 2>"$R_LOG_DIR/sasa.log"

# 7. SNAPSHOTS
echo "=== EXTRACTING SNAPSHOTS ==="
mkdir -p "$RESULTS_DIR/snapshots"
echo "0" | gmx trjconv \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -o "$RESULTS_DIR/snapshots/snapshots_.pdb" \
    -dt 1000 \
    -sep \
    -quiet 2>"$R_LOG_DIR/trjconv_snapshots.log"

# 8. RADIUS OF GYRATION
echo "=== RADIUS OF GYRATION ==="
echo "1" | gmx gyrate \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -o "$RESULTS_DIR/gyrate.xvg" \
    -quiet 2>"$R_LOG_DIR/gyrate.log"

# 9. SECONDARY STRUCTURE
# FIX: gmx do_dssp is deprecated in GROMACS 2025; use gmx dssp instead
echo "=== SECONDARY STRUCTURE ==="
gmx dssp \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -o "$RESULTS_DIR/ss.dat" \
    -num "$RESULTS_DIR/ss_count.xvg" \
    -quiet 2>"$R_LOG_DIR/dssp.log"

# 10. HYDROGEN BOND MAP
# gmx hbond (new) does not support -hbm; use gmx hbond-legacy for the XPM map.
# LIG_GROUP was determined dynamically during make_ndx above.
echo "=== HYDROGEN BOND MAP ==="
echo -e "1\n${LIG_GROUP}" | gmx hbond-legacy \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -n "$WORKDIR/index.ndx" \
    -hbm "$RESULTS_DIR/hbmap.xpm" \
    -quiet 2>"$R_LOG_DIR/hbmap.log"

# 11. CONTACTS (minimum distance protein-ligand)
echo "=== CONTACT ANALYSIS ==="
echo -e "1\n${LIG_GROUP}" | gmx mindist \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -n "$WORKDIR/index.ndx" \
    -od "$RESULTS_DIR/contacts.xvg" \
    -quiet 2>"$R_LOG_DIR/contacts.log"

# 12. PCA
echo "=== PCA ==="
echo "1" | gmx covar \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -o "$RESULTS_DIR/eigenvalues.xvg" \
    -v "$RESULTS_DIR/eigenvec.trr" \
    -quiet 2>"$R_LOG_DIR/covar.log"

echo "1" | gmx anaeig \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -v "$RESULTS_DIR/eigenvec.trr" \
    -first 1 -last 2 \
    -proj "$RESULTS_DIR/projection.xvg" \
    -quiet 2>"$R_LOG_DIR/anaeig.log"

# 13. REDUCED TRAJECTORY
echo "=== REDUCED TRAJECTORY ==="
echo "1" | gmx trjconv \
    -s "$RESULTS_DIR/md.tpr" \
    -f "$WORKDIR/md.xtc" \
    -o "$RESULTS_DIR/md_prod_reduced.xtc" \
    -dt 100 \
    -quiet 2>"$R_LOG_DIR/trjconv_reduced.log"

echo "=== ADDITIONAL ANALYSIS COMPLETE ==="

echo "=== GENERATING SUMMARY REPORT ==="

SUMMARY_FILE="$RESULTS_DIR/simulation_summary.txt"

cat > "$SUMMARY_FILE" << EOF
===============================================================================
SIMULATION SUMMARY REPORT
===============================================================================
Job ID:           $SLURM_JOB_ID
PDB File:         $(basename "$PDB_FILE")
System:           Protein-Ligand Complex
Date:             $(date)

===============================================================================
SYSTEM COMPOSITION
===============================================================================
Total atoms:      $FINAL_ATOMS
Water molecules:  $(grep -c "SOL" "$WORKDIR/complex_final.gro")
Ions:             NA:$(grep -c " NA " "$WORKDIR/complex_final.gro") CL:$(grep -c " CL " "$WORKDIR/complex_final.gro")
Box dimensions:   $BOX_INFO

Ligand info:
  Charge:          $CHARGE (method: $CHARGE_METHOD)
  Parameterization: $([ "$BCC_SUCCESS" = true ] && echo "AM1-BCC/GAFF2" || echo "Gasteiger/GAFF2")

===============================================================================
SIMULATION PROTOCOL
===============================================================================
Minimization:     steep descent, emtol=1000 kJ/mol/nm
NVT:              200 ps, V-rescale, 300K, POSRES + POSRES_LIG
NPT:              200 ps, Berendsen, 300K/1bar, POSRES + POSRES_LIG
Production:       10 ns, Parrinello-Rahman, 300K/1bar

===============================================================================
PERFORMANCE
===============================================================================
Total walltime:   $(echo "$SECONDS / 3600" | bc -l | xargs printf "%.2f") hours

===============================================================================
FILES GENERATED
===============================================================================
Trajectory:           results/trajectory.xtc
Final structure:      results/final_structure.gro
Energy file:          results/energies.edr
Index file:           results/index.ndx

Analysis:
  rmsd_protein.xvg    : Protein backbone RMSD (ref: production start)
  rmsd_ligand.xvg     : Ligand RMSD (ref: production start)
  rmsf_residues.xvg   : Per-residue flexibility
  hbonds_num.xvg      : Protein-ligand H-bonds over time
  sasa_protein.xvg    : Solvent accessible surface area
  contacts.xvg        : Minimum protein-ligand distance
  ss_count.xvg        : Secondary structure over time
===============================================================================
EOF

echo "Summary saved to: $SUMMARY_FILE"
echo "=== ALL DONE ==="
