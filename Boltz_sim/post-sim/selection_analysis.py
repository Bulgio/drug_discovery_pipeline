"""

"""
import os
import json
import pandas as pd


def json_reader(file_path : str):
# Reads the json into a dict structure, can be called in a for loop
    name=os.path.basename(file_path)

    with open(file_path,'r') as f:
        data = json.load(f)
    return data

def find_compound_by_id(compounds_data, identifier):
    """Find compound by identifier in list of compounds"""
    # Check if compounds_data is a list of dictionaries
    if isinstance(compounds_data, list):
        for compound in compounds_data:
            # Handle both dictionary compound format and CID string format
            if isinstance(compound, dict):
                if compound.get('id') == identifier:
                    return True
            elif str(compound) == str(identifier):
                return True
    # Also check if compounds_data is a DataFrame or Series
    elif hasattr(compounds_data, 'isin'):
        return identifier in compounds_data.values
    return False

def csv_reader(csv_path, max_lines=500):
    """
    Function to call to read the ligand csv and return the data needed to name and build the yamls
    Input:
    csv_path: path to the csv called
    max_lines: number of lines saved
    """
    try:
        # Read the CSV file
        parameters = pd.read_csv(csv_path, usecols=["smiles", "CIDs"], nrows=max_lines)
        return parameters
    except Exception as e:
        print(f"Error reading {csv_path}: {e}")
        return pd.DataFrame()

def data_origin(csv_path, sel_list, out_path,print_flag=True):
    """
    Function to look for the csv of origin of the compounds that made it in the top 20 of the ranking.
    """
    counter = []

    # Get all CSV files in the directory
    csv_files = [f for f in os.listdir(csv_path) if f.endswith('.csv')]
    
    for csv_file in csv_files:
        csv_full_path = os.path.join(csv_path, csv_file)
        
        db_counter = 0
        try:
            # Read compounds from CSV
            db_compounds_df = csv_reader(csv_full_path)
            
            # Extract CIDs from the DataFrame
            if not db_compounds_df.empty and 'CIDs' in db_compounds_df.columns:
                # Convert to list of CIDs (as strings)
                cid_list = db_compounds_df['CIDs'].dropna().astype(str).tolist()
                
                # Check each selected compound against the CID list
                for sel_compound in sel_list:
                    compound_id = str(sel_compound.get('id', ''))
                    if compound_id in cid_list:
                        db_counter += 1
            
            csv_name = csv_file.replace('similar_', '')
            
            counter_info = {
                'DrugBank_ID': csv_name,
                'Top_20_counter': db_counter,
                'is_present': 'yes' if db_counter > 0 else 'no'
            }
            counter.append(counter_info)
            
        except Exception as e:
            print(f"Error processing {csv_file}: {e}")
            continue
    if print_flag:
    # Save results to JSON file
        output_file = os.path.join(out_path, 'compound_origins.json')
        with open(output_file, 'w') as file:
            json.dump(counter, file, indent=2)

        print(f"\nFile with information about the top 20 saved to {output_file}")
    
    # Print summary
    print("\nSummary of compound origins:")
    for item in counter:
        if item['is_present'] == 'yes':
            print(f"  {item['DrugBank_ID']}: {item['Top_20_counter']} compounds")

    return counter

def inclusion_selection(list_origin,sort_list,sel_list,csv_path,out_path,integ_num=1,):
    """
    Docstring for inclusion_selection
    
    :param list_origin: Origin list of selection list
    :param sort_list: List of complete results
    :param sel_list: List of selections
    :param csv_path: Path to the directory
    :param integ_num: Number of compounds to be integrated for DB not included in the ranking std = 1
    """
    # Get all CSV files in the directory
    csv_files = [f for f in os.listdir(csv_path) if f.endswith('.csv')]
    db_with_selections = [item['DrugBank_ID'] for item in list_origin if item['is_present'] == 'yes']

    new_inclusions = []

    for csv_file in csv_files:
        csv_full_path = os.path.join(csv_path, csv_file)
        csv_name = csv_file.replace('similar_', '')

        if csv_name not in db_with_selections:
            included_count = 0
            print(f"\nProcessing database without representation: {csv_name}")
            # Read compounds from CSV
            db_compounds_df = csv_reader(csv_full_path)    
                # Extract CIDs from the DataFrame
            if not db_compounds_df.empty and 'CIDs' in db_compounds_df.columns:
                # Convert to list of CIDs (as strings)
                cid_list = db_compounds_df['CIDs'].dropna().astype(str).tolist()
                for cid in cid_list:
                    # Find this compound in the sorted list
                    for sorted_compound in sort_list:
                        # Assuming sort_list contains dictionaries with 'id' field
                        if str(sorted_compound.get('id', '')) == str(cid):
                            # Check if not already in selection
                            already_selected = False
                            for sel_compound in sel_list:
                                if str(sel_compound.get('id', '')) == str(cid):
                                    already_selected = True
                                    break
                            
                            if not already_selected:
                                new_inclusions.append(sorted_compound)
                                included_count += 1
                                print(f"  Adding compound {cid} from {csv_name}")
                                
                                if included_count >= integ_num:
                                    break
                    
                    if included_count >= integ_num:
                        break
        # Add the new inclusions to the selection list
    updated_selection = sel_list + new_inclusions
    
    # Create a summary
    summary = {
        'original_selection_count': len(sel_list),
        'new_inclusions_count': len(new_inclusions),
        'updated_selection_count': len(updated_selection),
        'new_inclusions': new_inclusions,
        'databases_added': list(set([csv_file.replace('similar_', '') for csv_file in csv_files 
                                     if csv_file.replace('similar_', '') not in db_with_selections]))
    }
    
        # Save results to JSON file
    output_file = os.path.join(out_path, 'summary_integration.json')
    with open(output_file, 'w') as file:
        json.dump(summary, file, indent=2)

    output_file = os.path.join(out_path, 'updated_selection.json')
    with open(output_file, 'w') as file:
        json.dump(updated_selection, file, indent=2)

    return updated_selection, summary                

def add_smiles_to_compounds(compounds_list, csv_path,out_path,print_flag=True):
    """
    Add SMILES codes to each compound in the list by looking them up in CSV files.
    
    :param compounds_list: List of compound dictionaries
    :param csv_path: Path to directory containing CSV files with SMILES data
    :return: Updated list of compounds with SMILES codes added
    """
    # Get all CSV files in the directory
    csv_files = [f for f in os.listdir(csv_path) if f.endswith('.csv')]
    
    # Create a mapping of CIDs to SMILES from all CSV files
    cid_to_smiles = {}
    
    for csv_file in csv_files:
        csv_full_path = os.path.join(csv_path, csv_file)
        try:
            df = csv_reader(csv_full_path)
            if not df.empty and 'CIDs' in df.columns and 'smiles' in df.columns:
                for _, row in df.iterrows():
                    cid = str(row['CIDs'])
                    smiles = row['smiles']
                    if pd.notna(cid) and pd.notna(smiles):
                        cid_to_smiles[cid] = smiles
        except Exception as e:
            print(f"Error processing {csv_file} for SMILES lookup: {e}")
            continue
    
    print(f"\nLoaded SMILES for {len(cid_to_smiles)} unique compounds from CSV files")
    
    # Add SMILES to each compound in the list
    updated_compounds = []
    compounds_without_smiles = []
    
    for compound in compounds_list:
        compound_id = str(compound.get('id', ''))
        
        # Create a copy of the compound dictionary
        updated_compound = compound.copy()
        
        # Look up SMILES
        if compound_id in cid_to_smiles:
            updated_compound['smiles'] = cid_to_smiles[compound_id]
        else:
            updated_compound['smiles'] = None
            compounds_without_smiles.append(compound_id)
        
        updated_compounds.append(updated_compound)
    
    # Print summary
    if compounds_without_smiles:
        print(f"\nWarning: Could not find SMILES for {len(compounds_without_smiles)} compounds:")
        for cid in compounds_without_smiles[:10]:  # Show first 10
            print(f"  CID: {cid}")
        if len(compounds_without_smiles) > 10:
            print(f"  ... and {len(compounds_without_smiles) - 10} more")
    else:
        print(f"\nSuccessfully added SMILES to all {len(updated_compounds)} compounds")
    
    if print_flag:
    # Save results to JSON file
        output_file = os.path.join(out_path, 'updated_compounds.json')
        with open(output_file, 'w') as file:
            json.dump(updated_compounds, file, indent=2)

    return updated_compounds

def add_reference_to_compounds(compounds_list, csv_path, out_path, print_flag=True):
    """
    Add reference (DBXXXXXX code) to each compound in the list by looking up their origin in CSV files.
    
    :param compounds_list: List of compound dictionaries
    :param csv_path: Path to directory containing CSV files with compound data
    :param out_path: Output path for saving JSON file
    :param print_flag: Whether to print output summary and save file
    :return: Updated list of compounds with reference keys added
    """
    # Get all CSV files in the directory
    csv_files = [f for f in os.listdir(csv_path) if f.endswith('.csv')]
    
    # Create a mapping of CIDs to their DrugBank references
    cid_to_reference = {}
    
    for csv_file in csv_files:
        csv_full_path = os.path.join(csv_path, csv_file)
        try:
            # Extract DrugBank ID from filename (remove 'similar_' prefix and '.csv' extension)
            drugbank_id = csv_file.replace('similar_', '').replace('.csv', '')
            
            df = csv_reader(csv_full_path)
            if not df.empty and 'CIDs' in df.columns:
                for _, row in df.iterrows():
                    cid = str(row['CIDs'])
                    if pd.notna(cid):
                        # If CID already has a reference, append this one as well
                        if cid in cid_to_reference:
                            if drugbank_id not in cid_to_reference[cid]:
                                cid_to_reference[cid].append(drugbank_id)
                        else:
                            cid_to_reference[cid] = [drugbank_id]
        except Exception as e:
            print(f"Error processing {csv_file} for reference lookup: {e}")
            continue
    
    print(f"\nLoaded references for {len(cid_to_reference)} unique compounds from CSV files")
    
    # Add reference to each compound in the list
    updated_compounds = []
    compounds_without_reference = []
    
    for compound in compounds_list:
        compound_id = str(compound.get('id', ''))
        
        # Create a copy of the compound dictionary
        updated_compound = compound.copy()
        
        # Look up reference
        if compound_id in cid_to_reference:
            # Join multiple references with semicolon if compound appears in multiple databases
            updated_compound['reference'] = ';'.join(cid_to_reference[compound_id])
        else:
            updated_compound['reference'] = 'unknown'
            compounds_without_reference.append(compound_id)
        
        updated_compounds.append(updated_compound)
    
    # Print summary
    if compounds_without_reference:
        print(f"\nWarning: Could not find reference for {len(compounds_without_reference)} compounds")
        for cid in compounds_without_reference[:10]:  # Show first 10
            print(f"  CID: {cid}")
        if len(compounds_without_reference) > 10:
            print(f"  ... and {len(compounds_without_reference) - 10} more")
    else:
        print(f"\nSuccessfully added references to all {len(updated_compounds)} compounds")
    
    if print_flag:
        # Save results to JSON file
        output_file = os.path.join(out_path, 'updated_selection.json')
        with open(output_file, 'w') as file:
            json.dump(updated_compounds, file, indent=2)
        
        print(f"\nUpdated selection with references saved to: {output_file}")
    
    return updated_compounds

#Selection path
sel_path='Boltz_files/outputs/selection/top_20_compounds.json'
#All the sorted compounds analyzed
sort_path='Boltz_files/outputs/selection/sorted_data.json'
#Raw compound data
csv_path='output/Scripts/csv'
#Output path
out_path='Boltz_files/outputs'

#Reading of the results
selection=json_reader(sel_path)
sel_list=selection['top_20_ranking']
#num=selection['metadata']['total_compounds'] #Number of the total compounds screened
data=json_reader(sort_path)
sort_list=data['ranking']

list_origin = data_origin(csv_path,sel_list,out_path) #Look for the origin of the compounds

up_sel,summ = inclusion_selection(list_origin,sort_list,sel_list,csv_path,out_path)

data_origin(csv_path,up_sel,out_path,False)

final_compounds=add_smiles_to_compounds(up_sel,csv_path,out_path)

final_compounds_with_references = add_reference_to_compounds(final_compounds, csv_path, out_path)