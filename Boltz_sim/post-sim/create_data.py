"""
Script created to create the pdbs from the cif produced by Boltz2 simulations, from the top 20

Input files:
updated_selection.json
Example of the data inside:
  {
    "id": "16114585",
    "combined_score": 0.9963472885642035,
    "confidence_placing": 18,
    "affinity_placing": 13,
    "smiles": "CN1C(=O)C(N=C1N)(C2=CC=C(C=C2)OC(F)F)C3=CC=CC(=C3)C=CCCCO",
    "reference": "DB07874"
  }
keys: ["id","reference"] others are not important to this script
Directory to target: (Inside the cluster)
path: 
> /home/toscana/scratch/boltz2/affinity/outputs/similar_DBXXXXX/boltz_results_yyyyyyy_BACE1
The similar_DBXXXXX is a problem because it slows the identification of the compounds inside the cluster 
(it is not good in doing iteration in a lot of files)
> Resolved using the existing code

Directory structure:
out_dir/
├── lightning_logs/                                            # Logs generated during training or evaluation
├── predictions/                                               # Contains the model's predictions
    ├── [input_file1]/
        ├── [input_file1]_model_0.cif                          # The predicted structure in CIF format, with the inclusion of per token pLDDT scores
        ├── confidence_[input_file1]_model_0.json              # The confidence scores (confidence_score, ptm, iptm, ligand_iptm, protein_iptm, complex_plddt, complex_iplddt, chains_ptm, pair_chains_iptm)
        ├── affinity_[input_file1].json                        # The affinity scores (affinity_pred_value, affinity_probability_binary, affinity_pred_value1, affinity_probability_binary1, affinity_pred_value2, affinity_probability_binary2)

        ├── pae_[input_file1]_model_0.npz                      # The predicted PAE score for every pair of tokens
        ├── pde_[input_file1]_model_0.npz                      # The predicted PDE score for every pair of tokens
        ├── plddt_[input_file1]_model_0.npz                    # The predicted pLDDT score for every token
        ...
        └── [input_file1]_model_[diffusion_samples-1].cif      # The predicted structure in CIF format
        ...
    └── [input_file2]/
        ...
└── processed/                                                 # Processed data used during execution 

For the execution of the script Biopython is necessary, so pip install biopython and/or module load biopython is necessary to run it. 
Validation for the use of Biopython to convert mmCif intop PDB https://link.springer.com/article/10.1186/s12859-023-05388-9
"""
import os
import json
import tarfile
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from pathlib import Path

def json_reader(file_path: str):
    """Reads JSON file into a dictionary."""
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data

def tar_extractor(tar_path,input_dir):
        # Check if tar file exists
    if not os.path.exists(tar_path):
        print(f"Error: Tar file not found at {tar_path}")
        return
    
    # Extract the tar file
    print(f"Extracting {tar_path}...")
    try:
        with tarfile.open(tar_path, 'r:gz') as tar:
            tar.extractall(input_dir)
        print("Extraction complete.")
    except Exception as e:
        print(f"Error extracting tar file: {e}")
        return

def cif_to_pdb_converter(cif_path: str, pdb_id: str, output_dir: str):
    """Convert mmCIF file to PDB format."""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create unique structure name
    name = f"{pdb_id}_structure"
    parser = MMCIFParser()
    io = PDBIO()
    
    try:
        structure = parser.get_structure(name, cif_path)

        rename_chains_to_pdb_compliant(structure)

        io.set_structure(structure)
        
        pdb_name = f"{name}.pdb"
        pdb_path = os.path.join(output_dir, pdb_name)
        io.save(pdb_path)
        
        print(f"Converted {cif_path} to {pdb_path}")
        return structure
    except Exception as e:
        print(f"Error converting {cif_path}: {e}")
        return None

def rename_chains_to_pdb_compliant(structure):
    """
    Rename chain IDs to be PDB-compliant (single character).
    Handles multi-character chain IDs by mapping them to A-Z, a-z, 0-9.
    """
    # PDB-compliant chain IDs (62 possible characters)
    valid_chain_ids = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
    
    for model in structure:
        # Get all chain IDs in this model
        chain_ids = [chain.id for chain in model]
        
        # Create a mapping from old chain IDs to new ones
        chain_mapping = {}
        used_ids = set()
        
        for i, old_id in enumerate(chain_ids):
            if len(old_id) == 1 and old_id in valid_chain_ids:
                # Already PDB-compliant, keep it if not used
                if old_id not in used_ids:
                    chain_mapping[old_id] = old_id
                    used_ids.add(old_id)
                else:
                    # Need to find a new ID
                    for new_id in valid_chain_ids:
                        if new_id not in used_ids:
                            chain_mapping[old_id] = new_id
                            used_ids.add(new_id)
                            break
            else:
                # Multi-character or invalid chain ID, map to available ID
                for new_id in valid_chain_ids:
                    if new_id not in used_ids:
                        chain_mapping[old_id] = new_id
                        used_ids.add(new_id)
                        break
        
        # Apply the mapping
        for old_id, new_id in chain_mapping.items():
            if old_id != new_id:
                chain = model[old_id]
                chain.id = new_id
                # Also update the chain's full_id
                chain.full_id = (chain.full_id[0], chain.full_id[1], new_id)
        
        # If we couldn't map all chains (unlikely with 62 slots), warn
        if len(chain_mapping) < len(chain_ids):
            print(f"Warning: Could not map all {len(chain_ids)} chains to PDB-compliant IDs")


def main():
    # Use forward slashes or raw strings for Windows paths
    input_dir = 'Boltz_files/outputs'
    output_dir = 'Boltz_files/pdb_outputs'  # Separate output directory
    tar_file = 'selection.tar.gz'
    
    # Create directories if they don't exist
    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    tar_path = os.path.join(input_dir, tar_file)
    
    # tar_extractor(tar_path)
    
    # Process extracted content
    for entry in os.scandir(input_dir):
        # Skip the original tar file
        if entry.name.endswith('.tar.gz'):
            continue
        
        # Process JSON files if needed
        if entry.name.endswith('.json'):
            # data = json_reader(entry.path)
            # Process JSON data here if needed
            continue
        
        # Process directories (likely containing CIF files)
        if entry.is_dir():
            dir_name = entry.name
            parts = dir_name.split('_')
            base_id = parts[2] if parts else dir_name
            placeholder = parts[2] + '_' + parts[3]
            e_path = os.path.join(entry.path,placeholder)
            # Find all CIF files in the directory
            cif_files = [f for f in os.listdir(e_path) 
                        if f.endswith('.cif')]
            
            print(f"Found {len(cif_files)} CIF files in {dir_name}")
            
            # Convert each CIF file
            for idx, cif_file in enumerate(cif_files):
                cif_path = os.path.join(e_path, cif_file)
                unique_id = f"{base_id}_{idx}"
                cif_to_pdb_converter(cif_path, unique_id, output_dir)
    
    print("Processing complete.")

if __name__ == "__main__":
    main()