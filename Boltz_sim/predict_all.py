import os
import subprocess
from pathlib import Path

# Define absolute paths
BASE_DIR = Path("/home/toscana/scratch/boltz2/affinity")  
INPUT_DIR = BASE_DIR / "yaml"
OUTPUT_DIR = BASE_DIR / "outputs"
LOG_DIR = BASE_DIR / "logs"
cache_dir = "/home/toscana/.boltz"

# Check if directories exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

print(f"Base directory: {BASE_DIR}")
success_count = 0
for entry in os.scandir(INPUT_DIR):
    
    if entry.is_dir():
        log_file_e = LOG_DIR / f"{entry.name}.log"  
        print(f"\n Process Directory: {entry}")
        yaml_files = [f.name for f in os.scandir(entry.path) 
                     if f.is_file() and f.name.endswith(('.yaml', '.yml'))]
        if len(yaml_files) != 0:
            print(f"  Found {len(yaml_files)} YAML files")
            #Process .yaml files in the directory
            entry_counter = 0

            for file in yaml_files:
                
                print(f"Found file YAML {file}")
                # BEGIN PROCESSING FILES
                print(f"Processing: {file}")
                file_path = os.path.join(entry.path,file)
                # Build the command
                cmd = [
                    "boltz", 
                    "predict", 
                    file_path,  # Use full path or relative path
                    "--cache", cache_dir,
                    "--no_kernels",
                    "--out_dir", str(OUTPUT_DIR)
                ]

                # Create a unique log file for THIS specific YAML file
                log_file = LOG_DIR / f"{entry.name}_{Path(file).stem}.log"

                # Still print to console (optional, for monitoring)
                print(f"Processing {file} - see {log_file} for details")

                print(f"Running command: {' '.join(cmd)}")

                #Run the command
                try:
                    result = subprocess.run(
                    cmd,
                    capture_output=True,  # Capture stdout and stderr
                    text=True,            # Return strings, not bytes
                    check=False,          # Raise exception if command fails (False, because if not, a wrong yaml crashes everything)
                    timeout=7200          # Timeout maximum (2 hs per .yaml, in reality it does have a lot less because the whole job has a timeout in the .sh script)  
                    )
                except subprocess.TimeoutExpired:
                    print(f"TIMEOUT processing {file} after 7200 seconds")
                    with open(log_file, 'w') as log:
                        log.write(f"Command: {' '.join(cmd)}\n")
                        log.write("Return code: TIMEOUT\n")
                        log.write("="*60 + "\n")
                        log.write("TIMEOUT EXPIRED after 7200 seconds\n")
                    continue  # Skip to next file 

                if result.returncode == 0:

                    print("Command successful!")
                    success_count += 1
                    entry_counter += 1
                    print(f"Output:\n{result.stdout}")
                    print(f"Code result: {result.returncode}")
                    
                else:
                    print(f"Command failed with code {result.returncode}")
                    print(f"Error output:\n{result.stderr}")
                
                # Open the log file for writing
                with open(log_file, 'w') as log:
                    # Write all details to the log file
                    log.write(f"Processing: {file}\n")

                    # Save ALL details to log
                    log.write(f"Command: {' '.join(cmd)}\n")
                    log.write(f"Return code: {result.returncode}\n")
                    log.write("="*60 + "\n")
                    log.write("STDOUT from boltz:\n" + result.stdout + "\n")
                    log.write("="*60 + "\n")
                    log.write("STDERR from boltz:\n" + result.stderr + "\n")
                    log.write("="*60 + "\n")
        else:
            print(f"  No YAML files in {entry.name}")
            with open(log_file_e, 'w') as log:
                log.write(f"No YAML files found in {entry.name}\n")
    
        
        # Open the log file for writing
        with open(log_file_e, 'w') as log:
            # Write all details to the log file
            log.write(f"Processed: {entry.name}\n")

            # Save ALL details to log
            log.write(f"YAML file success: {entry_counter}\n")
            log.write("="*60 + "\n")
            log.write(f"YAML file failures: {len(yaml_files) - entry_counter} ""\n")
            log.write("="*60 + "\n")

print(f"\nTOTAL: {success_count} files processed successfully")
