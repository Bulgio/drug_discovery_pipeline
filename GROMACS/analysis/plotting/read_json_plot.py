"""

"""

import matplotlib
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