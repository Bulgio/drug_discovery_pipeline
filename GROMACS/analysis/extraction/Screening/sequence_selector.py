"""
Script to ensure that the sequence is the one assigned by UNIPROT as BACE1 human.
Code: P56817 · BACE1_HUMAN
Operations performed: 
>Download the sequence using the url API request. 
>Reads the JSON file downloaded, printing the sequences and various important informations
>Select the mature chain from the database
>Finds the pocket important residues
"""
import os
from os.path import exists
import shutil
import pycurl
import zipfile
import io
import xml.etree.ElementTree as ET
import pandas as pd
import tempfile
from Script_finale.GROMACS.fix_rdkit import Chem
import csv
import re
from urllib.request import urlretrieve
from Bio import SeqIO
import requests
import json

def rest_request(url):
    response = requests.get(url)
    if response.status_code == 200:
        uniprot_data = response.json()
        print("Success! File downloaded")
    
        with open(os.path.join('output\sequence','P56817.json'),'w') as f:
            json.dump(uniprot_data,f,indent=2)
        print(f'Data saved!')
    else:
        print(f"Download Failed! HTTP Status Code: {response.status_code}")
    return uniprot_data

def get_active_site_positions(data):
    """Extract just the active site positions"""
    positions = []
    
    for feature in data.get('features', []):
        if feature.get('type') == 'Active site':
            location = feature.get('location')
            if location:
                if 'position' in location:
                    pos = location['position'].get('value')
                    positions.append(pos)
                elif 'start' in location and 'end' in location:
                    start = location['start'].get('value')
                    end = location['end'].get('value')
                    if start == end:
                        positions.append(start)
                    else:
                        positions.append(start-end)
    
    return positions

def summary_entry(data):
    """
    Extracts key informations
    """
    chain_start = None
    for feature in data.get('features', []):
        if feature.get('type') == 'Chain':
            chain_start = feature['location']['start']['value']
            chain_end = feature['location']['end']['value']
            break

    sequence = str(data['sequence']['value'])
    mature_sequence = sequence[chain_start-1:chain_end]
    active_sites = [pos - (chain_start)+1 for pos in get_active_site_positions(data)]
    summary = {
        'id' : data.get('uniProtkbId'),
        'organism' : data.get('organism',{}).get('scientificName'),
        'name' : data['proteinDescription']['recommendedName']['fullName']['value'],
        'full_sequence' : sequence,
        'l_fullseq' : len(sequence),
        'chain_start' : chain_start,
        'chain_end' : chain_end,
        'mature_sequence' : mature_sequence,
        'l_matureseq' : len(mature_sequence),
        'active_sites' : active_sites
    }

    print(f'Summary: {summary}')
    with open(os.path.join('output\sequence','P56817_summary.json'),'w') as f:
            json.dump(summary,f,indent=2)
    print(f'Data saved!')
    return summary

def verify_sequence(msa_path: str,uniprot_summary: dict):
     
    with open(msa_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        msa_query = lines[1]

    seq=uniprot_summary['mature_sequence']
    chain_start=uniprot_summary['chain_start']
    error_code = 999
    if msa_query == seq:
        print(f"The sequence of {uniprot_summary['name']} is mantained in the msa, the first residue is listed as number {chain_start+1} in the Uniprot database")
        return chain_start
    else:
        print(f'ERROR! The sequences of the msa and of the Uniprot do not match')
        return error_code
    

url = 'https://rest.uniprot.org/uniprotkb/P56817.json'
data = rest_request(url)
#Debug showing keys 
print("Top-level keys in data:")
for key in data.keys():
    print(f"-{key}")
    if key == 'features':
        features = data.get('features',[])
        for feature in features:
            print(f"    -{feature}")
#

n=verify_sequence('BACE1_human_mature_sequence_d8829_cleaned_25_filtered.a3m',summary_entry(data))
"""
End of the script, the references are enclosed in the JSON. 
"""