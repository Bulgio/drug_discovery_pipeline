#!/bin/bash
#SBATCH --job-name=boltz2_BACE1_1
#SBATCH --account=def-jtus
#SBATCH --time=48:00:00
#SBATCH --gpus=a100_3g.20gb:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/output_%j.txt
#SBATCH --error=logs/error_%j.txt
#SBATCH --mail-user=matteo.bulgini@studenti.polito.it
#SBATCH --mail-type=ALL

# Caricamento dei moduli 
module purge 
module load python/3.12.4
module load cuda/12.9
module load rdkit/2024.09.6

# Attiva l'ambiente virtuale

source $HOME/py312/bin/activate 

# Posizionamento nella cartella di attività

cd $HOME/scratch/boltz2/affinity

# cartelle di attività

INPUT_DIR=$HOME/scratch/boltz2/affinity/yaml
OUTPUT_DIR=$HOME/scratch/boltz2/affinity/outputs
LOG_DIR=$HOME/scratch/boltz2/affinity/logs

# Crea directory se non esistono
mkdir -p $OUTPUT_DIR
mkdir -p $LOG_DIR

# Bash checks

echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "GPU assigned: $CUDA_VISIBLE_DEVICES"
echo ""
echo "Checking model weights..."
ls -lh $HOME/.boltz/ || echo "ERROR: Model weights not found!"
echo ""

python job_scripts/predict_all.py

echo ""
echo "Job finished at: $(date)"

echo ""
echo "Starting compression of results..."
echo "========================================"

# Create timestamp for unique archive name
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
JOB_NAME="boltz2_BACE1"

# Compress outputs
if [ -d "$OUTPUT_DIR" ] && [ "$(ls -A $OUTPUT_DIR 2>/dev/null)" ]; then
    echo "Compressing outputs..."
    tar -czf "${JOB_NAME}_outputs_${TIMESTAMP}.tar.gz" -C "$(dirname $OUTPUT_DIR)" "$(basename $OUTPUT_DIR)"
    echo "✓ Outputs compressed: ${JOB_NAME}_outputs_${TIMESTAMP}.tar.gz"
else
    echo "⚠ No outputs to compress"
fi

# Compress logs
if [ -d "$LOG_DIR" ] && [ "$(ls -A $LOG_DIR 2>/dev/null)" ]; then
    echo "Compressing logs..."
    tar -czf "${JOB_NAME}_logs_${TIMESTAMP}.tar.gz" -C "$(dirname $LOG_DIR)" "$(basename $LOG_DIR)"
    echo "✓ Logs compressed: ${JOB_NAME}_logs_${TIMESTAMP}.tar.gz"
else
    echo "⚠ No logs to compress"
fi

echo "========================================"
echo "Compression complete!"