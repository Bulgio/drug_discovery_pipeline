#!/bin/bash
#SBATCH --job-name=GROMACS_launcher
#SBATCH --account=def-jtus
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/home/toscana/scratch/new_GROMACS/logs/launch/output_%j.txt
#SBATCH --error=/home/toscana/scratch/new_GROMACS/logs/launch/error_%j.txt
#SBATCH --mail-user=matteo.bulgini@studenti.polito.it
#SBATCH --mail-type=ALL

# Module loading
module purge 
module load python/3.12.4
module load cuda/12.9
module load rdkit/2024.09.6
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
module load gromacs/2025.4

# Preparatory steps

GEN_INPUT_DIR="/home/toscana/scratch/new_GROMACS/input"
GEN_OUTPUT_DIR="/home/toscana/scratch/new_GROMACS/output"
JOB_PATH="/home/toscana/scratch/new_GROMACS/scripts/simulation.sh"
JOB_DIR="/home/toscana/scratch/new_GROMACS/scripts"
MDP_DIR="/home/toscana/scratch/new_GROMACS/input/mdp"

mkdir -p $GEN_OUTPUT_DIR
# FIND THE PDB FILES
find "$GEN_INPUT_DIR" -type f -name "*.pdb" > "$GEN_OUTPUT_DIR/pdb_files.txt" # Find all the pdbs files

while read pdb_file; do

    echo "PROCESSING: $pdb_file"
    filename=$(basename "$pdb_file" .pdb)
    code=$(echo "$filename" | cut -d'_' -f1)
    
    echo "PROCESSING: $pdb_file"

    # Create a corresponding output directory structure
    # Remove the INPUT_DIR part from the path
    output_subdir="$GEN_OUTPUT_DIR/$code"
    # Create the output directory

    # Check if the output directory already exists OR if a tar.gz file exists
    if [ -d "$output_subdir" ] || [ -f "${output_subdir}.tar.gz" ]; then
        echo "SKIPPING: Output directory (${output_subdir}) or tar.gz archive (${output_subdir}.tar.gz) already exists"
        echo "This job was likely already processed."
        echo "---"
        continue  # Skip to the next file
    fi

    mkdir -p "$output_subdir"
    job_name="gmx_${filename}"


    # Submit job with arguments
    sbatch \
        --job-name="$job_name" \
        --account=def-jtus \
        --time=36:00:00 \
        --gpus=a100_3g.20gb:1 \
        --cpus-per-task=8 \
        --mem=32G \
        "$JOB_PATH" "$pdb_file" "$output_subdir" "$filename" "$JOB_DIR" "$MDP_DIR"
    
    echo "Submitted job: $job_name for file: $filename"
    echo "Output directory: $output_subdir"
    echo "---"
    
done < "$GEN_OUTPUT_DIR/pdb_files.txt"

echo "All jobs submitted!"
