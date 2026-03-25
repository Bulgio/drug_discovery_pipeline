#!/bin/bash 
#SBATCH --mail-user=matteo.bulgini@studenti.polito.it
#SBATCH --mail-type=ALL


# Module loading
module purge
module load StdEnv/2023
module load gcc/12.3
module load python/3.12.4
module load cuda/12.9
module load openmpi/4.1.5
module load gromacs/2025.4
module load rdkit/2024.09.6
module load ambertools/25.0
module load openbabel/3.1.1

DIR_NAME="$1"
INP_FILE_DIR="$2"
NEC_FILES="$INP_FILE_DIR/file_necessary.txt"
BASE_DIR=$(dirname "$DIR_NAME")
LOG_DIR="$BASE_DIR/logs/$(basename "$DIR_NAME")"
VENV_PATH="/home/toscana/scratch/venv_mmpbsa/bin/activate"
export AMBERHOME=/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/MPI/gcc12/openmpi4/ambertools/25.0
if [[ -n "$SLURM_NTASKS" ]]; then
    NCPUS=$SLURM_NTASKS
else
    NCPUS=${SLURM_CPUS_PER_TASK:-1}
fi
mkdir -p $LOG_DIR

# Print job info
echo "=========================================="
echo "Job started at: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $DIR_NAME"
echo "Necessary files file: $NEC_FILES"
echo "=========================================="

# Add this:
echo "--- CPU/Task info ---"
echo "SLURM_NTASKS: ${SLURM_NTASKS:-not set}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-not set}"
echo "SLURM_JOB_CPUS_PER_NODE: ${SLURM_JOB_CPUS_PER_NODE:-not set}"
echo "NCPUS resolved to: $NCPUS"
echo "---------------------"

dos2unix $NEC_FILES

# Change to the working directory
cd "$DIR_NAME" || { echo "Cannot cd to $DIR_NAME"; exit 1; }
echo "Current working directory: $(pwd)"
# Read paths for MMPBSA
mapfile -t PATHS < "$NEC_FILES"

# Clean up paths (remove comments and empty lines)
CLEAN_PATHS=()
while IFS= read -r line; do
    [[ -z "$line" ]] && continue
    [[ "$line" =~ ^#.*$ ]] && continue
    # Remove carriage returns and trim whitespace
    line=$(echo "$line" | tr -d '\r' | xargs)
    CLEAN_PATHS+=("$line")
done < "$NEC_FILES"

# Make sure we have at least 4 paths
if [[ ${#CLEAN_PATHS[@]} -lt 4 ]]; then
    echo "ERROR: Not enough paths in $NEC_FILES. Need at least 4."
    echo "Found: ${CLEAN_PATHS[*]}"
    exit 1
fi

MMPBSA_INPUT="${CLEAN_PATHS[0]}"
TPR_FILE="${CLEAN_PATHS[1]}"
XTC_FILE="${CLEAN_PATHS[2]}"
NDX_FILE="${CLEAN_PATHS[3]}"

# Construct full paths (files should be in the working directory)
TPR_FILE="$DIR_NAME/$TPR_FILE"
XTC_FILE="$DIR_NAME/$XTC_FILE"
NDX_FILE="$DIR_NAME/$NDX_FILE"
MOL_FILE="$DIR_NAME/temp/ligand.mol2"

echo "=========================================="
echo "MM/PBSA input file: $MMPBSA_INPUT"
echo "TPR: $TPR_FILE"
echo "XTC: $XTC_FILE"
echo "NDX: $NDX_FILE"
echo "=========================================="

# Check if files exist
missing_files=0
for file in "$MMPBSA_INPUT" "$TPR_FILE" "$XTC_FILE" "$NDX_FILE"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file not found: $file"
        ((missing_files++))
    fi
done

if [[ $missing_files -gt 0 ]]; then
    echo "Missing $missing_files required files. Exiting."
    exit 1
fi

mapfile -t GROUP_NAMES < <(grep "^\[" "$NDX_FILE" | tr -d '[]' | sed 's/^ *//;s/ *$//')

PROT_GROUP=""
LIG_GROUP=""
COM_GROUP=""
for i in "${!GROUP_NAMES[@]}"; do
    name="${GROUP_NAMES[$i]}"
    [[ "$name" =~ Protein ]] && [[ -z "$PROT_GROUP" ]] && PROT_GROUP=$i
    [[ "$name" =~ MOL|LIG|UNK ]] && [[ -z "$LIG_GROUP" ]] && LIG_GROUP=$i
    [[ "$name" =~ non-Water ]] && [[ -z "$COM_GROUP" ]] && COM_GROUP=$i
done

if [[ -z "$PROT_GROUP" ]] || [[ -z "$LIG_GROUP" ]]; then
    echo "ERROR: Could not find protein or ligand group in $NDX_FILE"
    echo "Protein group: $PROT_GROUP"
    echo "Ligand group: $LIG_GROUP"
    echo "Available groups in index file:"
    grep -n "^\[" "$NDX_FILE"
    exit 1
fi

echo "Found groups - Protein: $PROT_GROUP, Ligand: $LIG_GROUP"

OUTPUT_FILE="$DIR_NAME/MMPBSA.dat"

# Ensure xtc output directory exists before trjconv
mkdir -p "$DIR_NAME/xtc"

# COM_GROUP fallback: use non-Water if available, otherwise System (0)
if [[ -z "$COM_GROUP" ]]; then
    COM_GROUP=$(echo "$NDX_PROBE" 2>/dev/null | awk '/^[ \t]*[0-9]+ non-Water[ \t]*:/ { match($0,/[0-9]+/); print substr($0,RSTART,RLENGTH); exit }')
fi
if [[ -z "$COM_GROUP" ]]; then
    COM_GROUP=0
    echo "⚠ COM_GROUP not found — using System (0) for PBC centering"
fi

echo -e "${PROT_GROUP}\n${COM_GROUP}" | gmx trjconv \
    -s "$TPR_FILE" \
    -f "$XTC_FILE" \
    -o "${DIR_NAME}/xtc/trajectory_pbc.xtc" \
    -n "$NDX_FILE" \
    -pbc mol \
    -center \
    -quiet 2>/dev/null

XTC_FILE="${DIR_NAME}/xtc/trajectory_pbc.xtc"

if [[ ! -f "${DIR_NAME}/xtc/trajectory_pbc.xtc" ]]; then
    echo "ERROR: trjconv PBC correction failed for $(basename $DIR_NAME)"
    exit 1
fi
echo "PBC-corrected trajectory generated"

# Activate virtual environment and run MMPBSA
if [[ -f "$VENV_PATH" ]]; then
    echo "Activating virtual environment: $VENV_PATH"
    source "$VENV_PATH" || { echo "Failed to activate virtual environment"; exit 1; }
else
    echo "WARNING: Virtual environment not found at $VENV_PATH"
fi

unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU SLURM_MEM_PER_NODE

echo "Running gmx_MMPBSA in parallel mode..."
srun -n $NCPUS gmx_MMPBSA MPI -O \
    -i "$MMPBSA_INPUT" \
    -cs "$TPR_FILE" \
    -ct "$XTC_FILE" \
    -ci "$NDX_FILE" \
    -cg $PROT_GROUP $LIG_GROUP \
    -lm "$MOL_FILE" \
    -o "$OUTPUT_FILE" \
    -eo "$DIR_NAME/MMPBSA.csv" \
    -nogui \
    2>&1 | tee "$LOG_DIR/mmpbsa_debug.log"

MMPBSA_EXIT_CODE=${PIPESTATUS[0]}

# Check if MMPBSA ran successfully
if [[ $MMPBSA_EXIT_CODE -eq 0 ]]; then
    echo "gmx_MMPBSA completed successfully"
else
    echo "gmx_MMPBSA failed with exit code $MMPBSA_EXIT_CODE"
fi

deactivate 2>/dev/null || true

echo "=========================================="
echo "Job completed at: $(date)"
echo "Exit code: $MMPBSA_EXIT_CODE"
echo "=========================================="

exit $MMPBSA_EXIT_CODE