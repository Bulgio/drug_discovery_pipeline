#!/bin/bash 
#SBATCH --mail-user=matteo.bulgini@studenti.polito.it
#SBATCH --mail-type=ALL


# Function to compress a single directory
compress_dir() {
    local dir_path="$1"
    local main_dir="$2"
    
    [[ -e "$dir_path" ]] || return 0
    [[ -d "$dir_path" ]] || return 0
    
    local dir_name=$(basename "$dir_path")
    local output_path="$main_dir/${dir_name}.tar.gz"
    
    # Skip if already compressed
    [[ -f "$output_path" ]] && return 0
    
    echo "[$$] Starting: $dir_name ($(du -sh "$dir_path" | cut -f1))"
    local start_time=$(date +%s)
    
    # Change to main_dir for compression to avoid path issues
    local current_dir=$(pwd)
    cd "$main_dir" || return 1
    
    if tar -czf "$output_path" \
        "$dir_name" \
        --record-size=1M \
        --blocking-factor=256 \
        --no-acls \
        --no-selinux \
        --no-xattrs 2>/dev/null; then
        
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        local archive_size=$(du -h "$output_path" | cut -f1)
        
        echo "[$$] ✓ $dir_name compressed in ${duration}s (${archive_size})"
        
        # Verify archive
        echo "[$$] Verifying $dir_name..."
        if tar -tzf "$output_path" &>/dev/null; then
            echo "[$$] ✓ Verification successful for $dir_name"
            # Uncomment to remove original after verification
            # rm -rf "$dir_path"
            # echo "[$$] Removed original $dir_name"
            cd "$current_dir" || return 1
            return 0
        else
            echo "[$$] ✗ Verification failed for $dir_name"
            rm -f "$output_path"
            cd "$current_dir" || return 1
            return 1
        fi
    else
        echo "[$$] ✗ Failed to compress $dir_name"
        rm -f "$output_path"
        cd "$current_dir" || return 1
        return 1
    fi
}

# Function to organize files based on extension
organize_dir() {
    local dir_path="$1"
    local ext="$2"
    
    # Validation
    if [[ ! -d "$dir_path" ]]; then
        echo "Error: Directory '$dir_path' does not exist" >&2
        return 1
    fi
    
    if [[ -z "$ext" ]]; then
        echo "Error: No extension provided" >&2
        return 1
    fi
    
    # Normalize extension (remove leading dot if present)
    ext="${ext#.}"
    local ext_pattern=".$ext"
    
    # Create extension directory name
    local target_dir="$dir_path/$ext"
    
    # Find all files with the given extension in the directory 
    local files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(find "$dir_path" -type f -name "*$ext_pattern" \
        -not -path "$target_dir/*" \
        -print0 2>/dev/null)
    
    # Check if any files found
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "No *$ext_pattern files found in $dir_path"
        return 0
    fi
    
    # Count files
    local file_count=${#files[@]}
    echo "Found $file_count *$ext_pattern files"
    
    # Create target directory if it doesn't exist
    mkdir -p "$target_dir"
    
    # Move each file
    local moved_count=0
    for file in "${files[@]}"; do
        if [[ -n "$file" ]]; then
            local filename=$(basename "$file")
            mv "$file" "$target_dir/"
            echo "  Moved: $filename"
            ((moved_count++))
        fi
    done
    
    echo "✓ Moved $moved_count *$ext_pattern files to $ext/"
    
    if [[ $moved_count -gt 0 ]]; then
        local dir_size=$(du -sh "$target_dir" | cut -f1)
        echo "  Directory size: $dir_size"
    fi
}

check_file() {
    if [ ! -s "$1" ]; then
        echo "❌ ERROR: $1 is empty or missing!"
        return 1
    fi
    echo "✓ $1 looks good"
    return 0
}

# Function to process a directory with both organization and compression
process_directory() {
    local work_dir="$1"
    local extension="$2"
    
    echo "=== Processing directory: $work_dir for extension: $extension ==="
    
    organize_dir "$work_dir" "$extension"
    
    echo "=== Processing complete for $extension ==="
}


#Module loading 

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

#Initialization
DIR_NAME="$1"
INP_FILE_DIR="$2"
EXT_FILE="$INP_FILE_DIR/extension.txt"

echo "=========================================="
echo "Job started at: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $DIR_NAME"
echo "Extensions file: $EXT_FILE"
echo "=========================================="

dos2unix "$EXT_FILE"

# Change to the working directory
cd "$DIR_NAME" || { echo "Cannot cd to $DIR_NAME"; exit 1; }
echo "Current working directory: $(pwd)"

# Process each extension from the file
if [[ -f "$EXT_FILE" ]]; then
    while IFS= read -r extension || [[ -n "$extension" ]]; do
        # Skip empty lines and comments
        [[ -z "$extension" ]] && continue
        [[ "$extension" =~ ^#.*$ ]] && continue
        
        # Trim whitespace
        extension=$(echo "$extension" | xargs)
        
        echo "PROCESSING $extension files"
        process_directory "$DIR_NAME/results" "$extension"

    done < "$EXT_FILE"
else
    echo "WARNING: Extensions file not found: $EXT_FILE"
fi

compress_dir "$DIR_NAME/results" "$DIR_NAME"

