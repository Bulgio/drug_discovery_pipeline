#!/bin/bash

# ============================================================
# analysis.sh
# Post-simulation analysis for a completed protein-ligand MD.
#
# Usage:
#   sbatch (via analysis_launcher.sh) or directly:
#   bash analysis.sh <output_dir>
# ============================================================

module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
module load gromacs/2025.4

OUTPUT_DIR="${1%/}"

if [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: provide output directory as argument"
    exit 1
fi

RESULTS_DIR="$OUTPUT_DIR/results"
WORKDIR="$OUTPUT_DIR/temp"
R_LOG_DIR="$RESULTS_DIR/logs"

# ── sanity checks ────────────────────────────────────────────
for f in "$RESULTS_DIR/md.tpr" "$RESULTS_DIR/trajectory.xtc" \
          "$RESULTS_DIR/initial_system.gro" "$RESULTS_DIR/energies.edr"; do
    if [ ! -f "$f" ]; then
        echo "❌ Missing: $f — aborting"
        exit 1
    fi
done

echo "=== ANALYSIS: $OUTPUT_DIR ==="
echo "Job ID: $SLURM_JOB_ID"
mkdir -p "$R_LOG_DIR"

TPR="$RESULTS_DIR/md.tpr"
XTC="$RESULTS_DIR/trajectory.xtc"
GRO="$RESULTS_DIR/initial_system.gro"
EDR="$RESULTS_DIR/energies.edr"
NDX="$RESULTS_DIR/index.ndx"

# ── helper: skip if output already exists ────────────────────
# Usage: run_step "label" "output_file" <command...>
run_step() {
    local label="$1"
    local outfile="$2"
    shift 2
    if [ -f "$outfile" ]; then
        echo "⏭  SKIP: $label ($(basename "$outfile") already exists)"
        return 0
    fi
    echo "→  Running: $label"
    if "$@" >"$R_LOG_DIR/${label}.log" 2>&1; then
        echo "✓  Done:    $label"
    else
        echo "❌ FAILED:  $label — see $R_LOG_DIR/${label}.log"
        tail -5 "$R_LOG_DIR/${label}.log" | sed 's/^/   | /'
    fi
}

# ============================================================
# STEP 1 — MASTER INDEX FILE
# Always regenerate — cheap and avoids stale group numbers.
# Probe first to detect LIG/MOL, then write clean index.
# ============================================================
echo "=== CREATING MASTER INDEX FILE ==="

NDX_PROBE=$(echo "q" | gmx make_ndx -f "$GRO" -o /dev/null 2>&1 || true)

LIG_GROUP=$(echo "$NDX_PROBE" | awk '/^[ \t]*[0-9]+ LIG[ \t]*:/ \
    { match($0,/[0-9]+/); print substr($0,RSTART,RLENGTH); exit }')
MOL_GROUP=$(echo "$NDX_PROBE" | awk '/^[ \t]*[0-9]+ MOL[ \t]*:/ \
    { match($0,/[0-9]+/); print substr($0,RSTART,RLENGTH); exit }')

if [ -n "$LIG_GROUP" ]; then
    printf "q\n" | gmx make_ndx -f "$GRO" -o "$NDX" \
        -quiet 2>"$R_LOG_DIR/make_ndx.log"
    echo "✓ LIG group $LIG_GROUP found directly"
elif [ -n "$MOL_GROUP" ]; then
    printf "name %s LIG\nq\n" "$MOL_GROUP" | gmx make_ndx -f "$GRO" -o "$NDX" \
        -quiet 2>"$R_LOG_DIR/make_ndx.log"
    LIG_GROUP="$MOL_GROUP"
    echo "✓ MOL group $MOL_GROUP renamed to LIG"
else
    N_GROUPS=$(echo "$NDX_PROBE" | awk \
        '/^[ \t]*[0-9]+ [A-Za-z]/ { match($0,/[0-9]+/); n=substr($0,RSTART,RLENGTH) } \
        END { print n+0 }')
    LIG_GROUP=$((N_GROUPS + 1))
    printf "r LIG\nname %d LIG\nq\n" "$LIG_GROUP" | gmx make_ndx -f "$GRO" -o "$NDX" \
        -quiet 2>"$R_LOG_DIR/make_ndx.log"
    echo "✓ LIG added as new group $LIG_GROUP"
fi

if [ -z "$LIG_GROUP" ]; then
    echo "❌ ERROR: could not determine LIG group. Check $R_LOG_DIR/make_ndx.log"
    exit 1
fi

cp "$NDX" "$WORKDIR/index.ndx" 2>/dev/null || true

# ============================================================
# STEP 2 — PROTEIN BACKBONE RMSD
# ============================================================
echo "=== PROTEIN BACKBONE RMSD ==="
run_step "rmsd_protein" "$RESULTS_DIR/rmsd_protein.xvg" \
    bash -c "echo -e 'Backbone\nBackbone' | gmx rms \
        -s '$TPR' -f '$XTC' -n '$NDX' \
        -o '$RESULTS_DIR/rmsd_protein.xvg' -tu ns -quiet"

# ============================================================
# STEP 3 — LIGAND RMSD
# Fit on Backbone to remove protein drift, measure on LIG.
# Uses master index directly — avoids duplicate LIG group bug.
# ============================================================
echo "=== LIGAND RMSD ==="
run_step "rmsd_ligand" "$RESULTS_DIR/rmsd_ligand.xvg" \
    bash -c "echo -e 'Backbone\nLIG' | gmx rms \
        -s '$TPR' -f '$XTC' -n '$NDX' \
        -o '$RESULTS_DIR/rmsd_ligand.xvg' -tu ns -quiet"

# ============================================================
# STEP 4 — RMSF (per residue)
# ============================================================
echo "=== RESIDUE RMSF ==="
run_step "rmsf" "$RESULTS_DIR/rmsf_residues.xvg" \
    bash -c "echo 'Protein' | gmx rmsf \
        -s '$TPR' -f '$XTC' -n '$NDX' \
        -o '$RESULTS_DIR/rmsf_residues.xvg' -res -quiet"

# ============================================================
# STEP 5 — HYDROGEN BONDS + MAP
# Single hbond-legacy call produces both -num and -hbm.
# Replaces the broken gmx hbond (new) + separate hbond-legacy.
# ============================================================
echo "=== HYDROGEN BONDS ==="
run_step "hbonds" "$RESULTS_DIR/hbonds_num.xvg" \
    bash -c "echo -e '1\n${LIG_GROUP}' | gmx hbond-legacy \
        -s '$TPR' -f '$XTC' -n '$NDX' \
        -num '$RESULTS_DIR/hbonds_num.xvg' \
        -hbm '$RESULTS_DIR/hbmap.xpm' \
        -quiet"

# ============================================================
# STEP 6 — ENERGY COMPONENTS (autodiscovery)
# ============================================================
echo "=== ENERGY COMPONENTS ==="
if [ ! -f "$RESULTS_DIR/energy_components.xvg" ]; then
    EDR_TERMS=$(echo "0" | gmx energy -f "$EDR" 2>&1 || true)
    get_edr_index() {
        echo "$EDR_TERMS" | grep -iE "$1" | \
            awk 'NR==1 { match($0,/[0-9]+/); print substr($0,RSTART,RLENGTH) }'
    }
    ENERGY_IDX=""
    for term in "Potential" "Kinetic-En\." "Total-Energy" "Conserved-En\." \
                "Temperature" "^[[:space:]]+[0-9]+ Pressure" "Box-X" "Volume"; do
        idx=$(get_edr_index "$term")
        [ -n "$idx" ] && ENERGY_IDX="$ENERGY_IDX$idx "
    done
    ENERGY_IDX="${ENERGY_IDX% }"
    echo "  Thermodynamic indices: [$ENERGY_IDX]"
    if [ -n "$ENERGY_IDX" ]; then
        printf "%s\n0\n" "$(echo "$ENERGY_IDX" | tr ' ' '\n')" | gmx energy \
            -f "$EDR" -o "$RESULTS_DIR/energy_components.xvg" \
            -quiet 2>"$R_LOG_DIR/energy.log" \
        && echo "✓ energy_components.xvg" \
        || { echo "❌ energy failed — see $R_LOG_DIR/energy.log"
             tail -5 "$R_LOG_DIR/energy.log" | sed 's/^/   | /'; }
    else
        echo "⚠ No energy terms found — skipping"
    fi
else
    echo "⏭  SKIP: energy_components (already exists)"
fi

# ============================================================
# STEP 7 — SASA
# ============================================================
echo "=== SASA ==="
run_step "sasa" "$RESULTS_DIR/sasa_protein.xvg" \
    bash -c "echo 'Protein' | gmx sasa \
        -s '$TPR' -f '$XTC' -n '$NDX' \
        -o '$RESULTS_DIR/sasa_protein.xvg' \
        -or '$RESULTS_DIR/sasa_residue.xvg' -quiet"

# ============================================================
# STEP 8 — SNAPSHOTS (non-Water only — excludes solvent/ions)
# ============================================================
echo "=== EXTRACTING SNAPSHOTS ==="
mkdir -p "$RESULTS_DIR/snapshots"
run_step "snapshots" "$RESULTS_DIR/snapshots/snapshots_0.pdb" \
    bash -c "echo 'non-Water' | gmx trjconv \
        -s '$TPR' -f '$XTC' -n '$NDX' \
        -o '$RESULTS_DIR/snapshots/snapshots_.pdb' \
        -dt 1000 -sep -quiet"

# ============================================================
# STEP 9 — RADIUS OF GYRATION
# ============================================================
echo "=== RADIUS OF GYRATION ==="
run_step "gyrate" "$RESULTS_DIR/gyrate.xvg" \
    bash -c "echo 'Protein' | gmx gyrate \
        -s '$TPR' -f '$XTC' -n '$NDX' \
        -o '$RESULTS_DIR/gyrate.xvg' -quiet"

# ============================================================
# STEP 10 — SECONDARY STRUCTURE
# ============================================================
echo "=== SECONDARY STRUCTURE ==="
run_step "dssp" "$RESULTS_DIR/ss_count.xvg" \
    bash -c "echo 'Protein' | gmx dssp \
        -s '$TPR' -f '$XTC' -n '$NDX' \
        -o '$RESULTS_DIR/ss.dat' \
        -num '$RESULTS_DIR/ss_count.xvg' -quiet"

# ============================================================
# STEP 11 — CONTACTS (minimum protein-ligand distance)
# ============================================================
echo "=== CONTACT ANALYSIS ==="
run_step "contacts" "$RESULTS_DIR/contacts.xvg" \
    bash -c "echo -e '1\n${LIG_GROUP}' | gmx mindist \
        -s '$TPR' -f '$XTC' -n '$NDX' \
        -od '$RESULTS_DIR/contacts.xvg' -quiet"

# ============================================================
# STEP 12 — PCA (C-alpha, lighter than full Protein)
# ============================================================
echo "=== PCA ==="
run_step "pca_covar" "$RESULTS_DIR/eigenvalues.xvg" \
    bash -c "echo -e 'C-alpha\nC-alpha' | gmx covar \
        -s '$TPR' -f '$XTC' -n '$NDX' \
        -o '$RESULTS_DIR/eigenvalues.xvg' \
        -v '$RESULTS_DIR/eigenvec.trr' -quiet"

run_step "pca_proj" "$RESULTS_DIR/projection.xvg" \
    bash -c "echo 'C-alpha' | gmx anaeig \
        -s '$TPR' -f '$XTC' -n '$NDX' \
        -v '$RESULTS_DIR/eigenvec.trr' \
        -first 1 -last 2 \
        -proj '$RESULTS_DIR/projection.xvg' -quiet"

# ============================================================
# STEP 13 — REDUCED TRAJECTORY
# ============================================================
echo "=== REDUCED TRAJECTORY ==="
run_step "trjconv_reduced" "$RESULTS_DIR/md_prod_reduced.xtc" \
    bash -c "echo 'Protein' | gmx trjconv \
        -s '$TPR' -f '$XTC' -n '$NDX' \
        -o '$RESULTS_DIR/md_prod_reduced.xtc' \
        -dt 100 -quiet"

# ============================================================
# SUMMARY REPORT
# ============================================================
echo "=== GENERATING SUMMARY REPORT ==="
SUMMARY_FILE="$RESULTS_DIR/simulation_summary.txt"
FINAL_ATOMS=$(tail -n 1 "$GRO" | awk '{print $1}' 2>/dev/null || echo "N/A")

cat > "$SUMMARY_FILE" << EOF
===============================================================================
SIMULATION SUMMARY REPORT
===============================================================================
Job ID:           $SLURM_JOB_ID
Output directory: $OUTPUT_DIR
Date:             $(date)

===============================================================================
SYSTEM COMPOSITION
===============================================================================
Total atoms (complex_final): $FINAL_ATOMS
Ligand index group:          $LIG_GROUP (LIG)

===============================================================================
SIMULATION PROTOCOL
===============================================================================
Minimization: steep descent, emtol=1000 kJ/mol/nm
NVT:          200 ps, V-rescale, 300K, POSRES + POSRES_LIG
NPT:          200 ps, Berendsen, 300K/1bar, POSRES + POSRES_LIG
Production:   10 ns, Parrinello-Rahman, 300K/1bar

===============================================================================
FILES GENERATED (results/)
===============================================================================
index.ndx               : Master index (all groups + LIG)
rmsd_protein.xvg        : Protein backbone RMSD vs production start
rmsd_ligand.xvg         : Ligand RMSD (fit Backbone, measure LIG)
rmsf_residues.xvg       : Per-residue flexibility
hbonds_num.xvg          : Protein-ligand H-bond count over time
hbmap.xpm               : H-bond occupancy map
contacts.xvg            : Minimum protein-ligand distance
sasa_protein.xvg        : Solvent accessible surface area
sasa_residue.xvg        : Per-residue SASA
gyrate.xvg              : Radius of gyration
ss.dat / ss_count.xvg   : Secondary structure over time
eigenvalues.xvg         : PCA eigenvalues
projection.xvg          : PCA projection (PC1 vs PC2)
md_prod_reduced.xtc     : Reduced trajectory (every 100 ps)
snapshots/              : PDB snapshots every 1 ns (non-Water)
===============================================================================
EOF

# Quick stats
echo ""
if [ -f "$RESULTS_DIR/rmsd_protein.xvg" ]; then
    FINAL_RMSD=$(grep -v "^[@#]" "$RESULTS_DIR/rmsd_protein.xvg" | \
        tail -1 | awk '{printf "%.3f nm", $2}')
    echo "  Protein final RMSD : $FINAL_RMSD"
fi
if [ -f "$RESULTS_DIR/rmsd_ligand.xvg" ]; then
    FINAL_LIG=$(grep -v "^[@#]" "$RESULTS_DIR/rmsd_ligand.xvg" | \
        tail -1 | awk '{printf "%.3f nm", $2}')
    echo "  Ligand  final RMSD : $FINAL_LIG"
fi
if [ -f "$RESULTS_DIR/hbonds_num.xvg" ]; then
    AVG_HB=$(grep -v "^[@#]" "$RESULTS_DIR/hbonds_num.xvg" | \
        awk '{sum+=$2; n++} END {printf "%.2f", sum/n}')
    echo "  Avg H-bonds        : $AVG_HB"
fi

echo ""
echo "Summary saved to: $SUMMARY_FILE"
echo "=== ANALYSIS COMPLETE: $OUTPUT_DIR ==="
