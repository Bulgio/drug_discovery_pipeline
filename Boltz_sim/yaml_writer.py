"""
Script that writes the .yaml required to run Boltz predictions. 
In order to run this script msa is necessarly generated and cleaned before, while the sequence is extracted before from Uniprot and analyzed according to the relevant data. 
"""


import pandas as pd
import re
import os
from pathlib import Path
from Script_finale.GROMACS.fix_rdkit import Chem
import json
import yaml
import tarfile
import shutil
from textwrap import dedent

def verify_sequence(msa_path: str,prot_seq: str):
    """
    Function to verify the sequence in the MSA is the same as in the YAML files
    """
    with open(msa_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        msa_query = lines[1]

    seq=prot_seq
    error_code = 999
    if msa_query == seq:
        print(f"The sequence is mantained in the msa")
        return
    else:
        print(f'ERROR! The sequences do not match')
        return error_code
    
def json_reader(json_path:str,chain_id='A1'):
    """
    Function to read the json summary of uniprot
    Expected input:
    json_path : str Relative path to the json summary
    Expected output:
    seq: mature sequence of the protein
    p_residues: list of protein residues that are identified as the pockets
    """
    with open(json_path,'r') as f:
        data = json.load(f) #Load the data in a dictionary from a .json file
    #Debug correct: print the dictionary
    print(data)
    print(type(data))
    #Extract the necessary data
    seq = data.get('mature_sequence') #Data relative to the mature sequence
    sites = data.get('active_sites') #Data relative to the active sites, with correct numbering
    p_residues = []
    r = 2
    # Obtain a list of dictionaries with the binding sites enumerated correctly
    for site in sites:
        for i in range(-r, r + 1):
            p_residues.append([chain_id, site + i])

    return seq,p_residues
        
def pubchem_smiles_canonical(smiles_code):
    """
    Function to ensure that the smiles code is in canonical shape
    input: smiles code (string)
    """
    


    mol=Chem.MolFromSmiles(smiles_code)
    if mol is None:
        return None
    Chem.Kekulize(mol, clearAromaticFlags=True)
    return Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=False
    )

def csv_reader(csv_path,max_lines=1000):

    """
    Function to read the ligand cvs and return the data needed to name and build the yamls
    Input:
    csv_path: path to the csv called
    max_lines: number of lines saved
    """

    parameters = pd.read_csv(csv_path, usecols=["smiles", "CIDs"],nrows=max_lines) #Read max_lines number of each csv
    ids = parameters['CIDs']
    smiles = parameters['smiles']
    return ids,smiles

def protein_compiler(seq,msa_path,id='A1'):
    """
    Function that compiles the YAML file section for the protein, in the correct formatting
    Inputs:
    seq: FASTA sequence of the protein
    msa_path: relative path to the .a3m file inside the cluster
    id: letter identifier, A is standard as there's more cases of 1 protein simulation
    Outputs:
    protein_dict: Dictionary containing all the informations
    """

    keys = ['id','sequence','msa']

    #Normalize the path
    norm_path = os.path.normpath(msa_path)
    norm_path = norm_path.replace('\\','/')
    # now the path should be in posix rules
    values = [str(id),str(seq),str(msa_path)]
    protein_dict = {k: v for k, v in zip(keys, values)}
    print(protein_dict) #Debug check
    return protein_dict

def ligand_compiler(id,smiles):
    """
    Function to compile the ligand part of the YAML file correctly formatted
    Input:
    id : CIDs from Pubchem as identifiers of molecules
    smiles: smiles codes from Pubchem
    Output:
    ligand_dict: Dictionary containing all the informations
    """
    keys = ['id','smiles']

    s_code = pubchem_smiles_canonical(smiles)
    values = [str(id),str(s_code)]
    ligand_dict = {k: v for k, v in zip(keys, values)}
    print(ligand_dict) #Debug check
    return ligand_dict

def extract_msa_info(msa_path):
    """
    Extract reference sequence and numbering from .a3m file.
    Returns the sequence and a list of residue positions.
    Useful to verify pocket position is aligned
    """
    

    with open(msa_path, 'r') as f:
        lines = f.readlines()
    
    # Find the reference sequence (usually first sequence after header)
    ref_sequence = ""
    for i, line in enumerate(lines):
        if line.startswith('>') and i+1 < len(lines):
            # Next line should be the reference sequence
            ref_sequence = lines[i+1].strip()
            break
    
    # Generate residue numbers (starting from 1)
    residue_numbers = list(range(1, len(ref_sequence) + 1))
    
    return ref_sequence, residue_numbers

def write_yaml(protein, ligand, p_residues, out_directory):
    """
    Function that writes the whole YAML file, combining the precedently elaborated information. 
    It works with the Boltz 2 Standard Formatting, combining one ligand with one protein and predicting their affinity
    INPUT: 
    >Protein dict structure from protein_compiler
    >Ligand dict structure from protein compiler
    >p_residues list containing the pocket residues of the protein, from json_reader
    >out_directory path (string) for the directory in which to write the yaml 
    """

    out_path = os.path.join(out_directory, f"{ligand['id']}_BACE1.yaml")
    
    # Format contacts as flow style
    contacts_str = str(p_residues).replace("'", "").replace('"', '')
    
    ligand_id = "B1"

    # Build YAML with dedent - This allows nice code formatting
    yaml_content = dedent(f"""\
        version: 1
        sequences:
          - protein:
              id: {protein['id']}
              sequence: {protein['sequence']}
              msa: {protein['msa']}
          - ligand:
              id: {ligand_id}
              smiles: {ligand['smiles']}
        constraints:
          - pocket:
              binder: {ligand_id}
              contacts: {contacts_str}
        properties:
          - affinity:
              binder: {ligand_id}
    """)
    
    with open(out_path, 'w') as yaml_file:
        yaml_file.write(yaml_content)
    
    print(f"YAML file created: {out_path}")
    return out_path
"""
Section 0: Paths of the file realized until now
"""
#path to the csv directory
csv_dir = '..\csv'
list_csv = os.listdir(csv_dir)
#path to the json file
json_path = '..\\P56817_summary.json'
#path to read the msa file
msa_file = '..\\msa\\BACE1_cleaned.a3m'
#path to write in the .yaml files in order for them to work
msa_path = msa_file #In the case of automated msa compiling this is not necessary, but manually compiling the MSA is more reliable.
out_directory='..\\yaml'
"""
Section 1: Compiling of the parts
"""
seq,p_res=json_reader(json_path)
b=verify_sequence(msa_file,seq)
total_successful_count = 0
for csv in list_csv:
    successful_count = 0
    csv_path=os.path.join(csv_dir,csv)
    ids,smiles=csv_reader(csv_path,500)
    csv_basename=csv.removesuffix('.csv')
    if b is None:
        p_dict=protein_compiler(seq,msa_path)
        num=len(ids)
        # Create output directory if it doesn't exist
        os.makedirs(out_directory, exist_ok=True)
        # Create parent directory if it doesn't exist
        dir_cache=os.path.join(out_directory,csv_basename)
        os.makedirs(dir_cache, exist_ok=True)
        for id,smile in zip(ids,smiles):
            l_dict=ligand_compiler(id,smile)
            """
            Section 2: Writing of the yaml files
            """
            write_yaml(p_dict,l_dict,p_res,dir_cache)
            successful_count += 1
            total_successful_count += 1
        filename,ex=os.path.splitext(os.path.basename(csv_path))
        out_filename = os.path.join(out_directory, f"{filename}.tar.gz")
        if num >= 50:
            with tarfile.open(out_filename,"w:gz") as tar:
                tar.add(dir_cache,arcname=csv_basename)
                # Remove entire cache directory and all contents
            shutil.rmtree(dir_cache)
            print(f"Removed cache directory: {dir_cache}")
    print(f"The csv {csv} resulted in {successful_count} yaml files")
print(f"Among all csvs, a total of {total_successful_count} yaml files has been written")
            