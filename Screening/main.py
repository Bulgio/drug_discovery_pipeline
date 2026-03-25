"""
MAIN SCRIPT FOR THE SCREENING PHASE OF THE THESIS

INPUT:
> p_data            | Path to the file containing the id and password for DrugBank 
> db_path           | Path to the file containing the title of the databases required 
> ranking_path      | Path to the file containing the DBids of the manually selected compounds (subset for thesis work)
"""
"""
MAIN SCRIPT FOR THE SCREENING PHASE OF THE THESIS
Complete version with target polypeptide CSV parsing and XML drug extraction
"""

import os
import shutil
import pycurl
import zipfile
import io
import xml.etree.ElementTree as ET
import pandas as pd
import tempfile
import csv
import re
import json
from os.path import exists

# ============================================================================
# DATABASE DOWNLOAD FUNCTIONS
# ============================================================================

def download_databases(db_path, p_data, force_download=False):
    """
    Download both full database and target polypeptide IDs CSV
    """
    
    # Create database directory if it doesn't exist
    db_directory = "database"
    os.makedirs(db_directory, exist_ok=True)
    
    # Read credentials
    with open(p_data, "r", encoding='latin-1', errors='replace') as pdatafile:
        data = pdatafile.readlines()
        username = data[0].strip('\n')
        password = data[1].strip('\n')
    
    # Read database types from db_path
    with open(db_path, 'r') as file:
        datalines = file.readlines()
    
    downloaded_files = []  # Will store (db_type, file_path)
    
    for line in datalines:
        url_base = line.strip()
        
        # Map URL base to filename and type
        if url_base == "all-full-database":
            output_filename = "drugbank_full_database.zip"
            db_type = "full_database"
            description = "Complete database"
        elif url_base == "target-all-polypeptide-ids":
            output_filename = "drugbank_target_polypeptides.zip"
            db_type = "target_polypeptides"
            description = "Target polypeptide IDs (CSV)"
        else:
            output_filename = f"{url_base}.zip"
            db_type = url_base
            description = url_base
        
        db_full_path = os.path.join(db_directory, output_filename)
        
        # Check if already exists
        if os.path.exists(db_full_path) and not force_download:
            print(f"✓ {description} already exists: {output_filename}")
            downloaded_files.append((db_type, db_full_path))
            continue
        
        print(f"\nDownloading: {description}")
        print(f"Output file: {output_filename}")
        
        url = f"https://go.drugbank.com/releases/5-1-13/downloads/{url_base}"
        
        try:
            # Initialize cURL object
            c = pycurl.Curl()
            c.setopt(c.URL, url)
            c.setopt(c.USERPWD, f"{username}:{password}")
            c.setopt(c.FOLLOWLOCATION, True)
            c.setopt(c.FAILONERROR, True)
            c.setopt(c.NOPROGRESS, False)
            
            # Download file
            with open(output_filename, 'wb') as file:
                c.setopt(c.WRITEDATA, file)
                c.perform()
            
            # Check HTTP status code
            http_code = c.getinfo(c.RESPONSE_CODE)
            c.close()
            
            if http_code == 200:
                # Move file to database directory
                shutil.move(output_filename, db_full_path)
                print(f"✓ Download completed: {output_filename}")
                downloaded_files.append((db_type, db_full_path))
            else:
                print(f"✗ Failed to download. HTTP Code: {http_code}")
                
        except pycurl.error as e:
            print(f"✗ Error downloading: {e}")
            if os.path.exists(output_filename):
                os.remove(output_filename)
        except Exception as e:
            print(f"✗ Unexpected error: {e}")
            if os.path.exists(output_filename):
                os.remove(output_filename)
    
    return downloaded_files

# ============================================================================
# TARGET DATABASE PARSING (CSV)
# ============================================================================

def extract_drugbank_ids_from_target_db(zip_path, prot_code):
    """
    Search for protein in target polypeptides database CSV
    Returns list of DrugBank IDs that target this protein
    """
    found = False
    drugbank_ids = []
    
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"\nExtracting target database...")
        
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(temp_dir)
        
        # Find CSV file
        csv_files = [f for f in os.listdir(temp_dir) if f.endswith('.csv')]
        if not csv_files:
            print("No CSV file found in the archive")
            return False, []
        
        csv_path = os.path.join(temp_dir, csv_files[0])
        print(f"Reading CSV file: {csv_files[0]}")
        print(f"Searching for Uniprot ID: {prot_code}")

        # Read the CSV file
        df = pd.read_csv(csv_path)
        print(f"CSV columns: {list(df.columns)}")

        # Column names (based on your output)
        uniprot_col = 'UniProt ID'
        db_col = 'Drug IDs'
        
        # Check if columns exist
        if uniprot_col not in df.columns:
            print(f"Error: Column '{uniprot_col}' not found in CSV")
            print(f"Available columns: {list(df.columns)}")
            return False, []
        
        if db_col not in df.columns:
            print(f"Error: Column '{db_col}' not found in CSV")
            return False, []
        
        # Find rows where Uniprot ID matches our protein
        matches = df[df[uniprot_col].astype(str).str.contains(prot_code, na=False, case=False)]
        
        if not matches.empty:
            found = True
            print(f"\n✓ Found {len(matches)} entries for protein {prot_code}")
            
            # Extract DrugBank IDs from the matching rows
            for idx, row in matches.iterrows():
                print(f"\nEntry {idx + 1}:")
                print(f"  Name: {row.get('Name', 'N/A')}")
                print(f"  Gene Name: {row.get('Gene Name', 'N/A')}")
                print(f"  UniProt ID: {row.get(uniprot_col, 'N/A')}")
                print(f"  Species: {row.get('Species', 'N/A')}")
                
                # Get Drug IDs - they are semicolon-separated
                drug_ids = row.get(db_col, '')
                if pd.notna(drug_ids) and drug_ids:
                    # Split by semicolon (;) to get individual DrugBank IDs
                    ids = str(drug_ids).split(';')
                    # Clean up each ID (remove whitespace)
                    clean_ids = [id.strip() for id in ids if id.strip()]
                    drugbank_ids.extend(clean_ids)
                    print(f"  DrugBank IDs found: {len(clean_ids)}")
                    # Show first 5 IDs as sample
                    for drug_id in clean_ids[:5]:
                        print(f"    - {drug_id}")
                    if len(clean_ids) > 5:
                        print(f"    ... and {len(clean_ids) - 5} more")
            
            # Remove duplicate DrugBank IDs
            drugbank_ids = list(set(drugbank_ids))
            print(f"\n✓ Total unique DrugBank IDs found: {len(drugbank_ids)}")
            
            # Show all DrugBank IDs
            if drugbank_ids:
                print("\nAll DrugBank IDs targeting this protein:")
                # Print in groups of 5 for better readability
                for i in range(0, len(drugbank_ids), 5):
                    group = drugbank_ids[i:i+5]
                    print(f"  {', '.join(group)}")
        else:
            print(f"\n✗ Protein {prot_code} not found in the database")
            
            # Show a few examples of what UniProt IDs look like in the database
            print("\nFirst few UniProt IDs in database (for reference):")
            for idx, val in df[uniprot_col].head(5).items():
                if pd.notna(val):
                    print(f"  {val}")
    
    return found, drugbank_ids

# ============================================================================
# FULL DATABASE PARSING (XML)
# ============================================================================

def extract_drug_info_from_xml(zip_data, drugbank_ids):
    """
    Parse DrugBank XML and extract SMILES and other info for specific DrugBank IDs
    
    Parameters:
    -----------
    zip_data : bytes
        The in-memory ZIP file data containing the DrugBank XML
    drugbank_ids : list
        List of DrugBank IDs to search for
    
    Returns:
    --------
    pd.DataFrame : DataFrame containing drug information
    """
    
    # Create a BytesIO stream from the in-memory ZIP
    with zipfile.ZipFile(io.BytesIO(zip_data)) as z:
        # Lists all the files in the ZIP archive
        for file_info in z.infolist():
            if file_info.filename.endswith('.xml'):
                # Reading the XML file
                with z.open(file_info.filename) as xml_file:
                    # Parsing
                    tree = ET.parse(xml_file)
                    root = tree.getroot()
                    
                    # Define namespaces
                    if root.tag.startswith('{'):
                        namespace_uri = root.tag.split('}')[0].strip('{')
                        ns = {'db': namespace_uri}
                    else:
                        ns = {}  # No namespace

                    results = []
                    id_set = set(drugbank_ids)  # For faster lookup

                    # Find all the drug elements
                    for drug in root.findall('.//db:drug', ns):
                        # Get all DrugBank IDs for this drug
                        drugbank_id_elements = drug.findall('db:drugbank-id', ns)
                        
                        all_drug_ids = []
                        primary_id = None
                        
                        for drug_id in drugbank_id_elements:
                            id_value = drug_id.text
                            if id_value:
                                all_drug_ids.append(id_value)
                                if drug_id.get('primary') == 'true':
                                    primary_id = id_value
                        
                        # If no primary ID found, use the first ID
                        if not primary_id and all_drug_ids:
                            primary_id = all_drug_ids[0]
                        
                        # Check if any of this drug's IDs are in our target list
                        matching_ids = [did for did in all_drug_ids if did in id_set]
                        
                        if matching_ids:  # This drug matches our search
                            drug_data = {}
                            
                            # Get drug name
                            name_element = drug.find('db:name', ns)
                            drug_data['Name'] = name_element.text if name_element is not None else 'Unknown'
                            
                            # Add ID information
                            drug_data['Primary_ID'] = primary_id
                            drug_data['All_DrugBank_IDs'] = '; '.join(all_drug_ids)
                            
                            # Extract SMILES (from calculated properties)
                            smiles = ""
                            calc_props = drug.find('db:calculated-properties', ns)
                            if calc_props is not None:
                                for prop in calc_props.findall('db:property', ns):
                                    kind_elem = prop.find('db:kind', ns)
                                    if kind_elem is not None and kind_elem.text == 'SMILES':
                                        value_elem = prop.find('db:value', ns)
                                        if value_elem is not None:
                                            smiles = value_elem.text
                                            break
                            
                            # If SMILES not found in calculated properties, try external identifiers
                            if not smiles:
                                ext_ids = drug.find('db:external-identifiers', ns)
                                if ext_ids is not None:
                                    for ext_id in ext_ids.findall('db:external-identifier', ns):
                                        resource = ext_id.findtext('db:resource', '', ns)
                                        if 'SMILES' in resource:
                                            smiles = ext_id.findtext('db:identifier', '', ns)
                                            break
                            
                            drug_data['SMILES'] = smiles
                            
                            # Extract drug type
                            drug_data['Drug_Type'] = drug.get('type', 'unknown')
                            
                            # Extract description (if available)
                            description = drug.findtext('db:description', '', ns)
                            drug_data['Description'] = description[:200] + '...' if len(description) > 200 else description
                            
                            # Extract CAS number (if available)
                            cas_number = drug.findtext('db:cas-number', '', ns)
                            drug_data['CAS_Number'] = cas_number
                            
                            # Extract groups (approved, experimental, etc.)
                            groups = []
                            for group in drug.findall('db:groups/db:group', ns):
                                if group.text:
                                    groups.append(group.text)
                            drug_data['Groups'] = '; '.join(groups)
                            
                            results.append(drug_data)
                    
                    # Create DataFrame
                    if results:
                        return pd.DataFrame(results)
                    else:
                        return pd.DataFrame()

    return pd.DataFrame()

# ============================================================================
# SAVE RESULTS FUNCTIONS
# ============================================================================
def save_drug_data(df, prot_code, formats=['csv', 'json'], filters=None):
    """
    Save drug data to CSV and/or JSON files with optional filtering
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing drug information
    prot_code : str
        Protein code (for filename)
    formats : list
        List of output formats ('csv', 'json')
    filters : dict, optional
        Dictionary of filters to apply before saving
        Example: 
        filters = {
            'require_smiles': True,           # Only keep drugs with SMILES
            'drug_types': ['small molecule'],  # Only specific drug types
            'groups': ['approved'],            # Only approved drugs
            'exclude_without_cas': True,       # Exclude drugs without CAS number
            'custom_filter': function          # Custom function that returns boolean mask
        }
    
    Returns:
    --------
    list : List of saved filenames
    """
    db_directory = "drugbank_selection"
    os.makedirs(db_directory, exist_ok=True)
    saved_files = []
    
    if df.empty:
        print("No data to save")
        return saved_files
    
    # Make a copy to avoid modifying original
    filtered_df = df.copy()
    original_count = len(filtered_df)
    
    # Apply filters if provided
    if filters:
        print(f"\nApplying filters to {original_count} compounds...")
        
        # Filter 1: Require SMILES
        if filters.get('require_smiles', False):
            before = len(filtered_df)
            filtered_df = filtered_df[filtered_df['SMILES'].str.len() > 0]
            print(f"  - Require SMILES: removed {before - len(filtered_df)} compounds without SMILES")
        
        # Filter 2: Drug types
        if filters.get('drug_types'):
            drug_types = filters['drug_types']
            if isinstance(drug_types, str):
                drug_types = [drug_types]
            before = len(filtered_df)
            filtered_df = filtered_df[filtered_df['Drug_Type'].isin(drug_types)]
            print(f"  - Drug types {drug_types}: removed {before - len(filtered_df)} compounds")
        
        # Filter 3: Groups/Approval status
        if filters.get('groups'):
            groups = filters['groups']
            if isinstance(groups, str):
                groups = [groups]
            before = len(filtered_df)
            # Groups column contains semicolon-separated values
            mask = filtered_df['Groups'].apply(
                lambda x: any(group in str(x).split('; ') for group in groups) if pd.notna(x) else False
            )
            filtered_df = filtered_df[mask]
            print(f"  - Groups {groups}: removed {before - len(filtered_df)} compounds")
        
        # Filter 4: Exclude without CAS number
        if filters.get('exclude_without_cas', False):
            before = len(filtered_df)
            filtered_df = filtered_df[filtered_df['CAS_Number'].str.len() > 0]
            print(f"  - Require CAS number: removed {before - len(filtered_df)} compounds")
        
        # Filter 5: Custom filter function
        if filters.get('custom_filter'):
            custom_func = filters['custom_filter']
            before = len(filtered_df)
            filtered_df = filtered_df[custom_func(filtered_df)]
            print(f"  - Custom filter: removed {before - len(filtered_df)} compounds")
        
        # Filter 6: Remove by specific DrugBank IDs (blacklist)
        if filters.get('exclude_ids'):
            exclude_ids = filters['exclude_ids']
            before = len(filtered_df)
            filtered_df = filtered_df[~filtered_df['Primary_ID'].isin(exclude_ids)]
            print(f"  - Exclude specific IDs: removed {before - len(filtered_df)} compounds")
        
        # Filter 7: Keep only specific DrugBank IDs (whitelist)
        if filters.get('include_only_ids'):
            include_ids = filters['include_only_ids']
            before = len(filtered_df)
            filtered_df = filtered_df[filtered_df['Primary_ID'].isin(include_ids)]
            print(f"  - Include only specific IDs: removed {before - len(filtered_df)} compounds")
        
        # Filter 8: Molecular weight range (if you add MW to extraction)
        if filters.get('mw_range') and 'Molecular_Weight' in filtered_df.columns:
            min_mw, max_mw = filters['mw_range']
            before = len(filtered_df)
            filtered_df = filtered_df[
                (filtered_df['Molecular_Weight'] >= min_mw) & 
                (filtered_df['Molecular_Weight'] <= max_mw)
            ]
            print(f"  - MW range {min_mw}-{max_mw}: removed {before - len(filtered_df)} compounds")
        
        print(f"  → {len(filtered_df)} compounds remaining after filtering")
    
    if filtered_df.empty:
        print("\n⚠ No compounds left after filtering!")
        return saved_files
    
    base_filename = f"drugs_targeting_{prot_code}"
    
    # Add filter info to filename if filters were applied
    if filters:
        filter_suffix = "_filtered"
        base_filename = f"drugs_targeting_{prot_code}{filter_suffix}"
    
    # Save filtered data
    if 'csv' in formats:
        csv_file = os.path.join(db_directory,f"{base_filename}.csv")
        filtered_df.to_csv(csv_file, index=False)
        print(f"✓ CSV file saved: {csv_file}")
        saved_files.append(csv_file)
    
    if 'json' in formats:
        json_file =  os.path.join(db_directory,f"{base_filename}.json")
        filtered_df.to_json(json_file, orient='records', indent=2)
        print(f"✓ JSON file saved: {json_file}")
        saved_files.append(json_file)
    
    # Save simplified SMILES file (filtered)
    if 'Primary_ID' in filtered_df.columns and 'Name' in filtered_df.columns and 'SMILES' in filtered_df.columns:
        simple_df = filtered_df[['Primary_ID', 'Name', 'SMILES']].copy()
        simple_csv =  os.path.join(db_directory,f"{base_filename}_smiles_only.csv")
        simple_df.to_csv(simple_csv, index=False)
        print(f"✓ Simplified SMILES file saved: {simple_csv}")
        saved_files.append(simple_csv)
    
    # Save filtered DrugBank IDs list
    if 'Primary_ID' in filtered_df.columns:
        id_file =  os.path.join(db_directory,f"{base_filename}_ids.txt")
        with open(id_file, 'w') as f:
            for drug_id in sorted(filtered_df['Primary_ID'].tolist()):
                f.write(f"{drug_id}\n")
        print(f"✓ Filtered DrugBank IDs list saved: {id_file}")
        saved_files.append(id_file)
    
    # Also save the unfiltered data if you want to compare
    if filters and filters.get('save_unfiltered', False):
        unfiltered_base = f"drugs_targeting_{prot_code}_unfiltered"
        if 'csv' in formats:
            df.to_csv(os.path.join(db_directory,f"{unfiltered_base}.csv"), index=False)
            print(f"✓ Unfiltered CSV saved: {unfiltered_base}.csv")
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"FILTERED DRUG INFORMATION FOR PROTEIN: {prot_code}")
    print(f"{'='*60}")
    print(f"Original compounds: {original_count}")
    print(f"Filtered compounds: {len(filtered_df)}")
    print(f"Removed: {original_count - len(filtered_df)} compounds")
    
    if len(filtered_df) > 0:
        print(f"\nFiltered compounds statistics:")
        print(f"  Drugs with SMILES: {len(filtered_df[filtered_df['SMILES'].str.len() > 0])}")
        print(f"  Drug types: {filtered_df['Drug_Type'].value_counts().to_dict()}")
        
        print(f"\nFirst 5 filtered drugs:")
        for idx, row in filtered_df.head().iterrows():
            print(f"  {row['Primary_ID']}: {row['Name']}")
            if row['SMILES']:
                smiles_preview = row['SMILES'][:80] + "..." if len(row['SMILES']) > 80 else row['SMILES']
                print(f"    SMILES: {smiles_preview}")
    
    return saved_files
# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def main():
    # Define paths
    p_data = 'D:\\TESI\\Script_finale\\Screening\\input\\p_data.txt'
    db_path = 'D:\\TESI\\Script_finale\\Screening\\input\\db_path'
    
    # Get protein code from user
    prot_code = input("Insert the Standard UNIPROT code of the target: ")  # e.g., P56817
    
    # Check if input files exist
    if not os.path.exists(p_data):
        print(f"Error: Credentials file not found at {p_data}")
        return
    
    if not os.path.exists(db_path):
        print(f"Error: Database path file not found at {db_path}")
        return
    
    print("\n" + "="*60)
    print("STARTING DRUG SCREENING PIPELINE")
    print("="*60)
    
    # STEP 1: Download databases (only if needed)
    print("\n--- STEP 1: DATABASE DOWNLOAD ---")
    db_files = download_databases(db_path, p_data, force_download=False)
    
    if not db_files:
        print("Error: No databases were downloaded or found")
        return
    
    # Identify database files
    target_zip_path = None
    full_db_path = None
    
    for db_type, file_path in db_files:
        if db_type == "target_polypeptides":
            target_zip_path = file_path
            print(f"\n✓ Target mapping file: {os.path.basename(file_path)}")
        elif db_type == "full_database":
            full_db_path = file_path
            print(f"✓ Full database file: {os.path.basename(file_path)}")
    
    if not target_zip_path:
        print("\n⚠ Warning: Target polypeptide database not found")
        print("  Make sure 'target-all-polypeptide-ids' is in your db_path file")
        
    if not full_db_path:
        print("\n⚠ Warning: Full database not found")
        print("  Make sure 'all-full-database' is in your db_path file")
        return
    
    # STEP 2: Extract DrugBank IDs from target database
    print("\n--- STEP 2: EXTRACTING DRUGBANK IDs FROM TARGET DATABASE ---")
    found, drugbank_ids = extract_drugbank_ids_from_target_db(target_zip_path, prot_code)
    
    if not found or not drugbank_ids:
        print(f"\n✗ No DrugBank IDs found for protein {prot_code}")
        return
    
    # STEP 3: Extract drug information (including SMILES) from full database
    print("\n--- STEP 3: EXTRACTING DRUG INFORMATION FROM FULL DATABASE ---")
    print(f"Searching for {len(drugbank_ids)} DrugBank IDs in XML...")
    
    # Read the full database zip file into memory
    with open(full_db_path, 'rb') as f:
        zip_data = f.read()
    
    # Parse XML and extract drug information
    drugs_df = extract_drug_info_from_xml(zip_data, drugbank_ids)
    
    # STEP 4: Save results
    print("\n--- STEP 4: SAVING RESULTS ---")
    
    r_path='D:\\TESI\\Script_finale\\Screening\\input\\ranking.txt'
    with open(r_path, 'r') as f:
        whitelist = [line.strip() for line in f.readlines()]  

    filters_id = {
        'require_smiles': True,  
        'include_only_ids': whitelist
    }

    """
    Whitelist derived from the manual ranking, which was executed by analyzing the PDBs using
    MOE 2024, using the preliminary mode of calculating possible affinity to select only 10 compounds.
    In a production environment, only 'require_smiles' is given as an input.  
    """

    if not drugs_df.empty:
        saved_files = save_drug_data(
            drugs_df,
            prot_code,
            formats=['csv', 'json'],
            filters=filters_id 
            )
        
        print(f"\n✓ Pipeline completed successfully!")
        print(f"  Generated {len(saved_files)} output files:")
        for file in saved_files:
            print(f"  - {file}")
    else:
        print(f"\n✗ No drug information found in full database")
        print("\nTroubleshooting:")
        print("1. Check if the DrugBank IDs are correct")
        print("2. Verify the XML structure matches expected format")
        print("3. Try searching for a few specific IDs manually")
        
        # Still save the DrugBank IDs list
        id_file = f"drugbank_ids_{prot_code}.txt"
        with open(id_file, 'w') as f:
            for drug_id in sorted(drugbank_ids):
                f.write(f"{drug_id}\n")
        print(f"\n✓ DrugBank IDs list saved to {id_file}")

if __name__ == "__main__":
    main()