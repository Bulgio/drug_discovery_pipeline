"""
PubChem Similarity Query Script
Queries Pubchem database to search for similarity in respect of the csv of smiles it was given. Uses DrugBANK ids as identifiers (as names are bulky)

Input:
ranking_compounds_smiles.csv
structure: name,primary_id,all_ids,pdb_entries,SMILES
file containing infos on the compounds that were top ranked using an algorithm that considers GBVI, Internal energy and catalithic energies

Output:

Several csvs/other db formats of molecules with a similarity value to the compounds tested from DRUGBANK. 
WARNING: KEEP THE SIMILARITY VALUE HIGH AS THER MAY BE THOUSANDS OF MOLECULES VALID

Credits to https://colab.research.google.com/github/Ash100/DaS/blob/main/Similar_Compound_Search_PubChem.ipynb and Dr. Ashfaq Ahmad for the code snippets at the base of the code. 

"""

#Packages
import os
import time
from pathlib import Path
from urllib.parse import quote


import requests
import pandas as pd
from Script_finale.GROMACS.fix_rdkit import Chem
from rdkit.Chem import PandasTools

def pubchem_smiles_canonical(smiles_code):
    mol=Chem.MolFromSmiles(smiles_code)
    if mol is None:
        return None
    Chem.Kekulize(mol, clearAromaticFlags=True)
    return Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=False
    )

def query_pubchem_for_similar_compounds(smiles, threshold=90, n_records=0):
    """
    Query PubChem for similar compounds and return the job key.

    Parameters
    ----------
    smiles : str
        The canonical SMILES string for the given compound.
    threshold : int
        The threshold of similarity. In PubChem, the default threshold is 90%.
    n_records : int
        The maximum number of feedback records, here is disabled. !!

    Returns
    -------
    str
        The job key from the PubChem web service.
    """
    escaped_smiles = quote(smiles).replace("/", ".")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/{escaped_smiles}/JSON?Threshold={threshold}"
    r = requests.get(url)
    r.raise_for_status()
    key = r.json()["Waiting"]["ListKey"]
    return key

def check_and_download(key, attempts=30, delay=2, between=10):
    """
    Check job status and download PubChem CIDs when the job finished

    Parameters
    ----------
    key : str
        The job key of the PubChem service.
    attempts : int
        Number of attempts to check job status.
    delay : int
        Seconds before first check.
    between : int
        Seconds between subsequent checks.

    Returns
    -------
    list
        The PubChem CIDs of similar compounds.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{key}/cids/JSON"
    print(f"Querying for job {key} at URL {url}...", end="")
    # First Delay
    print(f"Waiting {delay} seconds before first check...")
    time.sleep(delay)
    for attempt in range(attempts):
        try:
            print(f"Attempt {attempt + 1}/{attempts}...", end="")
            r = requests.get(url,timeout=30)
            # Check for 400/404 errors
            if r.status_code == 400 or r.status_code == 404:
                    print(f"Job not ready or expired (HTTP {r.status_code})")
                    if attempt < attempts - 1:
                        time.sleep(between)
                        continue
                    else:
                        raise ValueError(f"Job key {key} expired after {attempts} attempts")
            r.raise_for_status() #Raise for other HTTP errors
            response = r.json()
            if "IdentifierList" in response:
                cids = response["IdentifierList"]["CID"]
                print(f"Success! Found {len(cids)} CIDs")
                return cids
            elif "Waiting" in response:
                #Job is processing...
                print(f"Job is still processing...")
            else:
                print(f"Unexpected response: {response.key()}")
        except requests.exceptions.RequestException as e:
            print(f"Network error: {e}")
            if attempt < attempts - 1:
                time.sleep(between)
                continue
        
        if attempt < attempts - 1:
            print(f"Waiting {between} seconds...")
            time.sleep(between)

    raise ValueError(f"Could not find matches for job key: {key} after {attempts} attempts")

def smiles_from_pubchem_cids(cids):
    """
    Get the canonical SMILES string from the PubChem CIDs.

    Parameters
    ----------
    cids : list
        A list of PubChem CIDs.

    Returns
    -------
    list
        The canonical SMILES strings of the PubChem CIDs.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{','.join(map(str, cids))}/property/CanonicalSMILES,IsomericSMILES,ConnectivitySMILES/JSON"
    r = requests.get(url)
    r.raise_for_status()
    smiles_list = []
    for item in r.json()["PropertyTable"]["Properties"]:
        # Use .get() to avoid KeyError if field is missing
        smiles = (item.get("CanonicalSMILES") or 
                  item.get("IsomericSMILES") or 
                  item.get("ConnectivitySMILES"))
        if smiles:  # Only add if not None/empty
            smiles_list.append(smiles)
        else:
            # Optionally log or handle missing SMILES
            smiles_list.append(None)  # or skip with continue
    
    return smiles_list

def batch_work_cids(cids,n_batch=10):

    """
    Docstring for batch_work_cids
    
    :param cids: List of cids from Pubchem
    :param n_batch: Number of cids in a singular selection
    """
    cid_batches = [cids[i:i+n_batch] for i in range(0, len(cids), n_batch)] #cid batches

    return cid_batches
    
## Section 0: Reading of csv
# Smiles parameters
csv_path_i = input("Insert the relative path of the .csv file containing the compounds I need to search for inside PubChem:")
csv_path=os.path.normpath(csv_path_i)
# reading csv file 
ranking_compounds = pd.read_csv(csv_path, usecols=["Primary_ID", "SMILES"])


for index, row in ranking_compounds.iterrows():
    target = row["SMILES"]  # Access SMILES column
    primary_id = row["Primary_ID"]  # Access primary_id column
    ## Section 1: Retrieval of queries 
    # Example:
    query =  pubchem_smiles_canonical(target)
    if query is None:
        print(f"Invalid SMILES skipped: {target}")
        continue
    # Request the job
    job_key = query_pubchem_for_similar_compounds(query,85)
    # Downloads similar CIDS from PubCHEM
    similar_cids = check_and_download(job_key,15)

    # Divide the cids into batches
    cids_batches = batch_work_cids(similar_cids,10)
    # After getting similar_cids
    print(f"Number of similar CIDs found: {len(similar_cids)}")
    #print(f"First 10 CIDs: {similar_cids[:10]}")

    smiles_list = []
    cids_list = []

    for num,batch in enumerate(cids_batches,start=1):
        # Check
        print(f"\n{'='*50}")
        progress = (num/len(cids_batches))*100
        print(f"Progress: {progress:.1f}%")
        print(f"\n{'='*50}")
        cids_list.extend(batch)
        batch_smiles = smiles_from_pubchem_cids(batch)
        #print(f"SMILES returned for batch: {batch_smiles}")  # Debug
        print(f"None values in batch: {sum(1 for s in batch_smiles if s is None)}")
        smiles_list.extend(batch_smiles)
        print(f"  ✓ Batch number {num} processed")
        time.sleep(1)  # Add delay to avoid rate limiting
    query_results_df = pd.DataFrame({"smiles":smiles_list,"CIDs":cids_list})
    query_results_df = query_results_df.dropna(subset=['smiles'])

    PandasTools.AddMoleculeColumnToFrame(query_results_df, smilesCol="smiles")
    query_results_df.head(6) #Shows 6 results
    #Save your data into CSV file
    query_results_df.to_csv(f'Pubchem_screening_results\\similar_{primary_id}.csv', index=True)