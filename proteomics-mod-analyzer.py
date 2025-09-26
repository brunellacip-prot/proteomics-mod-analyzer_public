# Main script for MaxQuant (MQ) evidence output file analysis 

#------------------------------------------------------------#
#  Uses MaxQuant (MQ) evidence output file to:
#    - Computes any modification extent and its statistics 
#              within group/s of sample replicates (_REP*)
#    - Returns:
#        > Extent of modification calculated for each protein separately (SPEM)
#        > Extent of modification calculated taking into account all proteins collectively (GEM)
#------------------------------------------------------------#
# Import libraries
import pandas as pd
import os
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
#------------------------------------------------------------#
# STEP 1 : Indicate input file and output folder paths

#   input file : MQ evidence file path provided as str
#        The evidence file should be in Excel format (.xlsx) or tab-delimited text (.txt). 
#        Should contain "Experiment", "Proteins", "Modifications" and "MS/MS count" columns.
#   output_dir : output folder path provided as str

input_file = r"C:\Users\brune\deamidation\Input_files\combined\txt\evidence.txt" #add correct path
output_dir = r"C:\Users\brune\deamidation\Results directory" #add correct path
#------------------------------------------------------------#
# STEP 2 : Load the MQ evidence output file (Excel or TXT)

def load_evidence_file(file_path):
    """
    Load MaxQuant evidence file - supports both Excel (.xlsx) and Tab-delimited (.txt) formats
    """
    file_extension = os.path.splitext(file_path)[1].lower()
    
    if file_extension == '.xlsx':
        df = pd.read_excel(file_path)
        print("🔹 Excel file loaded successfully")
    elif file_extension == '.txt':
        df = pd.read_csv(file_path, sep='\t', low_memory=False)
        print("🔹 TXT file loaded successfully")
    else:
        raise ValueError("Unsupported file format. Please use .xlsx or .txt files.")
    
    return df

maxquant_dataframe = load_evidence_file(input_file)

#------------------------------------------------------------#
# STEP 3 : Define key input parameters

#   target_aa : provided as str in letter code ("N", "Q", "S", "N|Q", etc...)
#       if you are targeting two amino acids at the same time please use "|" (e.g. "N|Q")
#   target_modification : provided as str (exactly as reported in the "Modifications" column)
#   separated_modified_aa : provided as str (exactly as reported in the "Modified sequence" column)
#       e.g. "N(Deamidation (NQ))", "Q(Deamidation (NQ))", "S(Phospho (ST))", "T(Phospho (ST)), etc..."
config_file = r"C:\Users\brune\deamidation\Input_modifications_MQ.xlsx" #add correct path

# # Fill with labels of modifications to be combined, grouping them according to the example.
combined_labels = [
    (["N(Deamidation (NQ))", "Q(Deamidation (NQ))"], "N|Q(Deamidation (NQ))")] #add correct labels to combine

#   uniprot_id : uniprot IDs (e.g. P02668, P02662, etc...) provided as list
uniprot_id = []  # Leave empty [] to extract all

# Define colors and samples sorting into plots

color_map = {
    'Sample_1': "#8cc5e3",
    'Sample_2': "#1a80bb",
    'Sample_3': "#d8a6a6",
    'Blank': "#b8b8b8"
    } #add correct sample group names and desidered colors for each of them

samples_sorted = ['Sample_1', 'Sample_2', 'Sample_3', 'Blank'] #add correct sample group sorting

mod_color_map = {
    "N(Deamidation (NQ))": "#1f77b4",  # blue
    "Q(Deamidation (NQ))": "#ff7f0e",  # orange
    "S(Phospho (ST))": "#2ca02c",      # green
    "T(Phospho (ST))": "#d62728"       # red
} #add correct color mapping for modifications

def extract_uniprot_ids_from_proteins(df):
    """
    Extract UniProt IDs from the 'Proteins' column.
    Handles formats like 'sp|P81265|PIGR_BOVIN' and extracts 'P81265'
    """
    all_proteins = df['Proteins'].dropna().unique()
    uniprot_ids = set()
    
    for protein_entry in all_proteins:
        # Split by semicolon in case of multiple proteins
        protein_list = str(protein_entry).split(';')
        
        for protein in protein_list:
            protein = protein.strip()
            
            # Pattern 1: sp|P81265|PIGR_BOVIN -> extract P81265
            if '|' in protein:
                parts = protein.split('|')
                if len(parts) >= 2:
                    potential_id = parts[1]
                    # Check if it looks like a UniProt ID (6 chars, starts with letter/number)
                    if len(potential_id) == 6 and potential_id[0].isalnum():
                        uniprot_ids.add(potential_id)
            else:
                # Pattern 2: Direct UniProt ID (P81265)
                if len(protein) == 6 and protein[0].isalnum():
                    uniprot_ids.add(protein)
    
    return sorted(list(uniprot_ids))

# If uniprot_id is empty, extract all UniProt IDs from the data
if not uniprot_id:
    uniprot_id = extract_uniprot_ids_from_proteins(maxquant_dataframe)
   
    print(f"🔹 Extracted {len(uniprot_id)} UniProt IDs from data: {uniprot_id[:5]}{'...' if len(uniprot_id) > 5 else ''}")
#------------------------------------------------------------#
#  Extract target_aa, target_modification and separated_modified_aa from Input_modifications_MQ Excel file

# Read the Excel file into a DataFrame
input_mod_file_MQ = pd.read_excel(config_file)

# Extract columns as lists, filtering out NaN values
target_aa = input_mod_file_MQ["target_aa"].dropna().tolist()
target_modification = input_mod_file_MQ["target_modification"].dropna().tolist()
separated_modified_aa = input_mod_file_MQ["separated_modified_aa"].dropna().tolist()

# Print to verify
print("Target AA:", target_aa)
print("Target Modifications:", target_modification)
print("Separated Modified AA:", separated_modified_aa)
print("Combined Labels:", combined_labels)
print("🔹 Input parameters loaded successfully")

#------------------------------------------------------------#
# STEP 4 : Read, sort and organize the experiments in the MQ evidence file

# Read and sort the experiments in alphabetical order

def read_and_sort_experiments(df):
    '''
    read_and_sort_experiments

    Read experiments from maxquant evidence dataframe 
    and sort them in alphabetical order.
    
    '''
    return sorted(df["Experiment"].unique())

experiments_list = read_and_sort_experiments(maxquant_dataframe)

#------------------------------------------------------------#
# Group the experiments by sample -> Create sample groups that contain replicates

def group_experiments_by_sample(experiments_list):
    """
    group_experiments_by_sample
    
    Groups experiments by sample name (removes _REP suffix).
    
    Returns a dictionary where keys are group names and values are lists of replicate names
        e.g. {'CAS': ['CAS_REP1', 'CAS_REP2', 'CAS_REP3'], ...}
    
    """
    from collections import defaultdict
    sample_groups = defaultdict(list)

    for exp in experiments_list:
        sample_name = exp.split('_REP')[0]
        sample_groups[sample_name].append(exp)
    
    return dict(sample_groups)

sample_groups = group_experiments_by_sample(experiments_list)

# if you want to get elements from the "sample_groups" dictionary:
#   for example: CAS_replicates = sample_groups["CAS"]

# Filter and order the MQ evidence file by sample group and protein

def filter_and_order_by_sample_group_protein(df, sample_groups, uniprot_id):
    ''' 
    filter_and_order_by_sample_group_protein

    Filter and order DataFrame by sample group and protein.

    df : maxquant_dataframe (MQ output evidence file)
        Input dataframe with 'Proteins' and 'Experiment' columns
    sample_groups : dict
        Dictionary where keys are group names and values are lists of replicate names
        e.g. {'CAS': ['CAS_REP1', 'CAS_REP2', 'CAS_REP3'], ...}
    uniprot_id : list
        List of protein identifiers to search for
        e.g. ["P02668", "P02662", "P02666", "P02663"]
    
    Returns pandas.DataFrame
        Filtered and ordered DataFrame containing all matching combinations    
    '''
    # Store all filtered combinations
    filtered_data = {}

    # Process each sample group
    for group_name, replicates in sample_groups.items():
        
        # Create pattern for matching any replicate in this group
        replicate_pattern = '|'.join(re.escape(rep) for rep in replicates)

        # Process each protein
        for protein_id in uniprot_id:

            # Filter by protein first
            df_protein = df[df["Proteins"].str.contains(protein_id, na=False, regex=False)]

            if not df_protein.empty:
                # Filter by sample group (any replicate)
                df_filtered = df_protein[df_protein["Experiment"].str.contains(replicate_pattern, na=False, regex=True)]
                
                if not df_filtered.empty:
                    # Add group_name and protein_id columns for easier sorting
                    df_filtered = df_filtered.copy()
                    df_filtered['Sample_Group'] = group_name
                    df_filtered['Protein_ID'] = protein_id

                    # Store in dictionary with tuple key
                    filtered_data[(group_name, protein_id)] = df_filtered
                    
    # Combine all filtered data
    if filtered_data:
        combined_df = pd.concat(filtered_data.values(), ignore_index=True)
        
        # Order by Sample_Group first, then by Protein_ID
        ordered_df = combined_df.sort_values(by=['Sample_Group', 'Protein_ID'], ascending=[True, True])
        ordered_df = ordered_df.reset_index(drop=True)
        
        return ordered_df
    else:
        # Return empty DataFrame with same columns if no matches found
        return pd.DataFrame(columns=df.columns.tolist() + ['Sample_Group', 'Protein_ID'])
    
maxquant_dataframe_by_sample_group_protein = filter_and_order_by_sample_group_protein(
    maxquant_dataframe, 
    sample_groups,
    uniprot_id
)

#------------------------------------------------------------#
# STEP 5 : Compute the number of targetable sites by experiment and protein (SP)
#          Return a Multi-index dataframe

# Compute targetable sites by experiment and protein -> Return single-protein targetable sites (SP)
# Organize targetable sites (SP) values in a Matrix with MultiIndex rows and columns

def compute_targetable_sites_SP(df, target_aa, return_individual=False):

    '''
    compute_targetable_sites_SP

    Compute targetable sites (SP) number from MS/MS counts and targetable amino acid counts.
    
    df : maxquant_dataframe ordered by sample_group and protein
        Input dataframe with 'MS/MS count' and 'Sequence' columns
    target_aa : list or str
        Target amino acids as strings in letter code (e.g. ["N", "Q"] or "NQ")
    return_individual : bool, default False
        If True, returns dict with individual amino acid results
        If False, returns total sum across all amino acids
    
    Returns:
    float or dict
        If return_individual=False: Total sum of (MS/MS count * # targetable sites)
        If return_individual=True: Dict with amino acid as key and targetable_sum as value
    '''
    
    MSMS_count = df["MS/MS count"]
    results = {}
    
    # Calculate for each amino acid
    for aa in target_aa:
        target_aa_count = df["Sequence"].str.count(aa)
        product_target_aa_MSMS = target_aa_count * MSMS_count
        results[aa] = product_target_aa_MSMS.sum()
    
    # Return based on the flag
    if return_individual:
        return results
    else:
        return sum(results.values())
    
    #-------------------------------------#

#   Organize targetable sites (SP) values in a Matrix with MultiIndex rows and columns

def create_matrix_targetable_sites_SP(df, target_aa, sample_groups=None, uniprot_id=None):
    '''
    create_matrix_targetable_sites_SP

    Create a matrix of targetable sites (SP) with Experiments/Sample_Groups as rows 
    and Protein_IDs + target amino acids as columns.
    
    df : pandas.DataFrame
        ordered MaxQuant dataset with columns: 'Experiment', 'Sample_Group', 'Protein_ID', 'MS/MS count', 'Sequence'
    target_aa : list or str
        Target amino acids as strings in letter code (e.g. ["N", "Q"] or "NQ")
    sample_groups : dict, optional
        Dictionary of sample groups.
    uniprot_id :
        List of protein IDs.
    
    Returns: 
    pandas.DataFrame
        Matrix with MultiIndex rows (Experiment, Sample_Group) and 
        MultiIndex columns (Protein_ID, Target_Amino_Acid)
    '''
    # Initialize results dictionary
    results = {}
    
    # Process each combination
    for sample_group, experiments in sample_groups.items():
        for experiment in experiments:
            # Check if this experiment actually exists in the data
            exp_data = df[df['Experiment'].str.contains(experiment, na=False)]
            if exp_data.empty:
                continue
                
            for protein_id in uniprot_id:
                # Filter data for this specific combination
                filtered_df = df[
                    (df['Sample_Group'] == sample_group) & 
                    (df['Experiment'].str.contains(experiment, na=False)) &
                    (df['Protein_ID'] == protein_id)
                ]
                
                if not filtered_df.empty:
                    # Calculate targetable sites for each amino acid separately
                    aa_results = compute_targetable_sites_SP(filtered_df, target_aa, return_individual=True)
                    
                    # Store results for each amino acid
                    for aa, value in aa_results.items():
                        results[(experiment, sample_group, protein_id, aa)] = value
                else:
                    # Set to 0 for missing combinations
                    for aa in target_aa:
                        results[(experiment, sample_group, protein_id, aa)] = 0
    
    # Convert results to DataFrame format
    if not results:
        print("No data found for the specified combinations.")
        return pd.DataFrame()
    
    # Get all unique experiment-sample_group combinations
    exp_sample_combinations = set()
    for (exp, sg, prot, aa), value in results.items():
        exp_sample_combinations.add((exp, sg))
    
    exp_sample_combinations = sorted(list(exp_sample_combinations))
    
    # Create MultiIndex for columns (Protein_ID, Target_Amino_Acid)
    column_combinations = []
    for protein_id in sorted(uniprot_id):
        for aa in sorted(target_aa):
            column_combinations.append((protein_id, aa))
    
    # Build the matrix
    data_matrix = []
    row_indices = []
    
    for exp, sample_group in exp_sample_combinations:
        row_indices.append((exp, sample_group))
        row_data = []
        
        for protein_id, aa in column_combinations:
            key = (exp, sample_group, protein_id, aa)
            value = results.get(key, 0)
            row_data.append(value)
        
        data_matrix.append(row_data)
    
    # Create DataFrame with MultiIndex for both rows and columns
    row_multi_index = pd.MultiIndex.from_tuples(row_indices, names=['Experiment', 'Sample_Group'])
    col_multi_index = pd.MultiIndex.from_tuples(column_combinations, names=['Protein_ID', 'Target_Amino_Acid'])
    
    result_df = pd.DataFrame(data_matrix, 
                            index=row_multi_index, 
                            columns=col_multi_index)
    
    return result_df

matrix_targetable_sites_SP = create_matrix_targetable_sites_SP(
    maxquant_dataframe_by_sample_group_protein, 
    target_aa,
    sample_groups,
    uniprot_id
)
print("🔹 Targetable sites (SP) matrix created")

# ---------------------------------------------------
# Refer to "utils.py" to get info about the DataFrames structure
# ---------------------------------------------------

# STEP 6 : Compute the number of modified sites by experiment and protein (SP)
#          Return a Multi-index dataframe

# Compute modified sites by experiment and protein -> Return single-protein modified sites (SP)
# Organize modified sites (SP) values in a Matrix with MultiIndex rows and columns

def create_weighted_modification_matrix_with_combined(
    df,
    separated_modified_aa,
    sample_groups,
    uniprot_id,
    combined_labels
):
    '''
    Create a matrix with weighted counts of individual modifications and their grouped combined sums.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input DataFrame with columns:
        'Modified sequence', 'MS/MS count', 'Experiment', 'Sample_Group', 'Protein_ID'
    
    separated_modified_aa : list of str
        List of all individual modifications to quantify.

    sample_groups : dict
        Dictionary of sample group names -> list of experiment names
    
    uniprot_id : list of str
        List of protein IDs to include

    combined_labels : list of tuples, optional
        Each tuple should be (list_of_mods, combined_label_str), e.g.:
        [ (["N(Deamidation (NQ))", "Q(Deamidation (NQ))"], "N|Q(Deamidation (NQ))") ]
    
    Returns:
    --------
    pandas.DataFrame
        Matrix with MultiIndex rows (Experiment, Sample_Group) and 
        MultiIndex columns (Protein_ID, Target_Modification)
    '''
    
    df["MS/MS count"] = pd.to_numeric(df["MS/MS count"], errors="coerce").fillna(0)
    results = {}

    for sample_group, experiments in sample_groups.items():
        for experiment in experiments:
            for protein_id in uniprot_id:
                # Filter subset
                subset = df[
                    (df["Sample_Group"] == sample_group) &
                    (df["Experiment"] == experiment) &
                    (df["Protein_ID"] == protein_id)
                ]

                # Track totals for each mod
                mod_totals = {}

                # Individual mods
                for mod in separated_modified_aa:
                    escaped_mod = re.escape(mod)
                    mod_count = subset["Modified sequence"].str.count(escaped_mod)
                    weighted_count = mod_count * subset["MS/MS count"]
                    total_weighted = weighted_count.sum()
                    mod_totals[mod] = total_weighted
                    results[(experiment, sample_group, protein_id, mod)] = total_weighted

                # Combined mods
                if combined_labels:
                    for mod_list, label in combined_labels:
                        combined_total = sum(mod_totals.get(m, 0) for m in mod_list)
                        results[(experiment, sample_group, protein_id, label)] = combined_total

    # Build index
    row_index = sorted({(exp, grp) for (exp, grp, _, _) in results})
    col_index = sorted({(prot, mod) for (_, _, prot, mod) in results})

    row_multi = pd.MultiIndex.from_tuples(row_index, names=["Experiment", "Sample_Group"])
    col_multi = pd.MultiIndex.from_tuples(col_index, names=["Protein_ID", "Target_Modification"])

    # Build matrix
    data = []
    for row in row_index:
        row_data = []
        for col in col_index:
            value = results.get((*row, *col), 0)
            row_data.append(value)
        data.append(row_data)

    return pd.DataFrame(data, index=row_multi, columns=col_multi)

matrix_modified_sites_SP = create_weighted_modification_matrix_with_combined(
    maxquant_dataframe_by_sample_group_protein,
    separated_modified_aa,
    sample_groups,
    uniprot_id,
    combined_labels
)
print("🔹 Modified sites (SP) matrix created")

# ---------------------------------------------------
# STEP 7 : Compute the global number of targetable sites by experiment (G)
#          Return a Multi-index dataframe

# Define functions necessary to compute sum number of #targetable sites for each experiment and target_aa
# Sum number of targetable sites for each experiment and target_aa (G)
#   For a given experiment and target_aa,
#   Compute sum of targetable sites (SP) 
#   Return targetable sites (G)

def find_matching_amino_acids(available_aas, target_aa):
    '''
    Find amino acids that match the target using flexible matching strategies.
    
    Strategies:
    1. Exact match
    2. Case-insensitive match  
    3. Partial match (for combined patterns like 'N|Q')
    4. Individual amino acid extraction from patterns
    '''
    matching_aas = []
    
    for aa in target_aa:
        # Strategy 1: Exact match
        exact_matches = [col for col in available_aas if col == aa]
        if exact_matches:
            matching_aas.extend(exact_matches)
            continue
            
        # Strategy 2: Case-insensitive exact match
        case_matches = [col for col in available_aas if str(col).upper() == aa.upper()]
        if case_matches:
            matching_aas.extend(case_matches)
            continue
            
        # Strategy 3: Handle combined patterns like 'N|Q' or 'S|T'
        if '|' in aa:
            individual_aas = aa.split('|')
            for individual_aa in individual_aas:
                individual_matches = [col for col in available_aas 
                                    if str(col).upper() == individual_aa.strip().upper()]
                matching_aas.extend(individual_matches)
        else:
            # Strategy 4: Partial match (column contains the amino acid)
            partial_matches = [col for col in available_aas 
                             if aa.upper() in str(col).upper()]
            matching_aas.extend(partial_matches)
    
    # Remove duplicates while preserving order
    return list(dict.fromkeys(matching_aas))

#-------------------------------------#

def handle_multiindex_columns_structure(df, exp, target_aa):
    '''
    Handle case where columns are MultiIndex, likely (Target_Amino_Acid, Protein)
    '''
    
    # Verify this is a MultiIndex row structure  
    if not isinstance(df.index, pd.MultiIndex):
        raise ValueError("DataFrame must have MultiIndex rows")
    
    # Find experiment level in index
    exp_level = 0
    for i, level_name in enumerate(df.index.names):
        if level_name and 'experiment' in str(level_name).lower():
            exp_level = i
            break
    
    # Filter rows by experiment
    exp_mask = df.index.get_level_values(exp_level).isin(exp)
    df_exp_filtered = df.loc[exp_mask]
    
    if df_exp_filtered.empty:
        print(f"Warning: No data found for experiment(s): {exp}")
        return pd.DataFrame()
    
    # Find target amino acid level in columns (usually level 0)
    aa_level = 0
    for i, level_name in enumerate(df_exp_filtered.columns.names):
        if level_name and ('amino' in str(level_name).lower() or 'aa' in str(level_name).lower() or 'target' in str(level_name).lower()):
            aa_level = i
            break
    
    # Find matching amino acid columns
    available_aas = df_exp_filtered.columns.get_level_values(aa_level).unique()
    matching_aas = find_matching_amino_acids(available_aas, target_aa)
    
    if not matching_aas:
        print(f"Warning: No amino acids found matching: {target_aa}")
        print(f"Available amino acids: {available_aas.tolist()}")
        return pd.DataFrame()
    
    # Filter columns by amino acids
    aa_mask = df_exp_filtered.columns.get_level_values(aa_level).isin(matching_aas)
    df_filtered = df_exp_filtered.loc[:, aa_mask]
    
    # Sum across proteins (assuming proteins are in level 1 of columns)
    # Group by amino acid level and sum across protein level
    result = df_filtered.T.groupby(level=aa_level).sum().T
    
    return result

#-------------------------------------#

def handle_single_level_columns_structure(df, exp, target_aa):
    '''
    Handle case where columns are single-level Target_Amino_Acid names
    '''
    
    # Verify this is a MultiIndex row structure
    if not isinstance(df.index, pd.MultiIndex):
        raise ValueError("DataFrame must have MultiIndex rows")
    
    # Find experiment level in index
    exp_level = 0
    for i, level_name in enumerate(df.index.names):
        if level_name and 'experiment' in str(level_name).lower():
            exp_level = i
            break
    
    # Filter rows by experiment
    exp_mask = df.index.get_level_values(exp_level).isin(exp)
    df_exp_filtered = df.loc[exp_mask]
    
    if df_exp_filtered.empty:
        print(f"Warning: No data found for experiment(s): {exp}")
        return pd.DataFrame()
    
    # Find matching columns
    matching_cols = find_matching_amino_acids(df_exp_filtered.columns, target_aa)
    
    if not matching_cols:
        print(f"Warning: No columns found matching: {target_aa}")
        print(f"Available columns: {df_exp_filtered.columns.tolist()}")
        return pd.DataFrame()
    
    # Filter and return
    return df_exp_filtered[matching_cols]

#-------------------------------------#

def compute_targetable_sites_G(df_targetable_SP, exp, target_aa):
    '''
    compute_targetable_sites_G

    Filter by experiment and target amino acid, then sum across all proteins
    for each (experiment, sample_group, target_amino_acid) combination.
    
    Keeps sample groups separate but sums across proteins.

    Returns:
    targetable sites (G) for any experiment and target_aa

    '''
    # Make copy to avoid modifying original data
    df = df_targetable_SP.copy()
    
    # Ensure exp and target_aa are lists for consistent handling
    if isinstance(exp, str):
        exp = [exp]
    if isinstance(target_aa, str):
        target_aa = [target_aa]
    
    # Check if columns are MultiIndex
    if isinstance(df.columns, pd.MultiIndex):
        if len(df.columns.levels) > 1:
            return handle_multiindex_columns_structure(df, exp, target_aa)
    else:
        return handle_single_level_columns_structure(df, exp, target_aa)
    
targetable_sites_G = compute_targetable_sites_G(
    matrix_targetable_sites_SP,
    experiments_list,
    target_aa
)

# ---------------------------------------------------
# STEP 8 : Compute the global number of modified sites by experiment (G)
#          Return a Multi-index dataframe

# Sum number of modified sites for each experiment and target_modification (G)
#   For a given experiment and target_modification,
#   Compute sum of modified sites (SP) 
#   Return modified sites (G)

def sum_modified_sites_across_proteins(df_modified_sites_SP, exp, target_aa):
    '''
    Filter by experiment and target amino acid, then sum across all proteins
    for each (experiment, sample_group, target_amino_acid) combination.
    
    Keeps sample groups separate but sums across proteins.
    
    Expected structure variations:
    1. Columns: MultiIndex (Target_Amino_Acid, Protein)
    2. Columns: Target_Amino_Acid, with proteins as values in cells
    3. Separate protein dimension
    '''
    
    # Make copy to avoid modifying original data
    df = df_modified_sites_SP.copy()
    
    # Ensure exp and target_aa are lists for consistent handling
    if isinstance(exp, str):
        exp = [exp]
    if isinstance(target_aa, str):
        target_aa = [target_aa]
    
    # Check if columns are MultiIndex
    if isinstance(df.columns, pd.MultiIndex):
        if len(df.columns.levels) > 1:
            return handle_multiindex_columns_structure(df, exp, target_aa)
    else:
        return handle_single_level_columns_structure(df, exp, target_aa)
    
modified_sites_G = sum_modified_sites_across_proteins(
    matrix_modified_sites_SP,
    experiments_list,
    target_aa
)
print("🔹 Modified sites and targetable sites (G) calculated successfully")

# ---------------------------------------------------
# STEP 9 : Compute ratios #modified sites (G) / #targetable sites (G)
#          Return a Multi-index dataframe

#Automated Mapping Code
def infer_mod_to_aa_map(mod_columns):
    """
    Automatically map mod column names to target amino acids.

    Example:
        "S|T(Phospho (ST))" → "S|T"

    Parameters:
        mod_columns (list or Index): mod column names as strings

    Returns:
        dict: mapping of {mod_column: target_amino_acid}
    """
    mod_to_aa_map = {}

    for col in mod_columns:
        if "(" in col:
            target_aa = col.split("(")[0]  # e.g., "S|T"
            mod_to_aa_map[col] = target_aa.strip()
    
    return mod_to_aa_map

def modified_targetable_G_ratio(mod_counts_df, aa_counts_df):
    """
    Compute ratio of #modified sites / #targetable sites for each modification.
    """
    # Infer AA mapping
    mod_to_aa_map = infer_mod_to_aa_map(mod_counts_df.columns)

    # Align index (safe even if aa_counts_df has missing entries)
    aa_counts_df = aa_counts_df.reindex(mod_counts_df.index)

    # Normalize
    normalized_df = mod_counts_df.copy()
    for mod_col, aa_col in mod_to_aa_map.items():
        if mod_col in mod_counts_df.columns and aa_col in aa_counts_df.columns:
            normalized_df[mod_col] = mod_counts_df[mod_col] / aa_counts_df[aa_col]
        else:
            raise KeyError(f"Missing column: '{mod_col}' or '{aa_col}' in inputs.")

    return normalized_df

ratios_G = modified_targetable_G_ratio(modified_sites_G, targetable_sites_G)

# ---------------------------------------------------
# STEP 10 : Calculate Single-Protein Extent of Modification (SPEM)

# Compute ratios #modified sites (SP) / #targetable sites (SP)
# Return a Multi-index dataframe

def modified_targetable_SP_ratio(df_modified_SP, df_targetable_SP):
    """ 
    Compute ratio of modified to targetable sites.
    
    Parameters:
    - df_modified_SP: DataFrame with MultiIndex rows (Experiment, Sample_Group),
                      and MultiIndex columns (Protein_ID, Target_Modification)
    - df_targetable_SP: DataFrame with same rows and MultiIndex columns 
                        (Protein_ID, Target_Amino_Acid)
    
    Returns:
    - ratio_df: DataFrame with same structure as df_modified_SP
    """

    def extract_residue_key(mod_label):
        """Extract AA key like 'N|Q' from 'N(Deamidation (NQ))' or 'N|Q(Deamidation (NQ))'"""
        prefix_match = re.match(r'^([A-Z|]+)\(', mod_label)
        if prefix_match:
            return prefix_match.group(1)
        inner_match = re.search(r'\(([^()]*)\)', mod_label)
        if inner_match:
            aas = re.findall(r'[A-Z]', inner_match.group(1))
            return '|'.join(sorted(set(aas)))
        return None

    ratio_df = pd.DataFrame(index=df_modified_SP.index, columns=df_modified_SP.columns, dtype="float")

    for (protein_id, mod_label) in df_modified_SP.columns:
        aa_key = extract_residue_key(mod_label)
        if not aa_key:
            continue

        aa_col = (protein_id, aa_key)

        if aa_col in df_targetable_SP.columns:
            modified = df_modified_SP[(protein_id, mod_label)]
            targetable = df_targetable_SP[aa_col]

            ratio = modified / targetable.replace(0, np.nan)
            ratio = ratio.replace([np.inf, -np.inf], np.nan)

            ratio_df[(protein_id, mod_label)] = ratio
        else:
            ratio_df[(protein_id, mod_label)] = np.nan

    return ratio_df

ratios_SP = modified_targetable_SP_ratio(matrix_modified_sites_SP, matrix_targetable_sites_SP)

# ---------------------------------------------------
# Calculate Single-Protein Extent of Modification (SPEM)
#   Summarizes ratios_SP by Sample_Group
#   Compute Mean (SPEM), Variance, Std_Dev, Std_Err (SEM)
#   Return a MultiIndexed DataFrame

def summarize_ratios_SP_data(df):
    """
    Summarizes ratios_SP values by Sample_Group and computes mean, variance, std dev, and SEM
    for each (Protein_ID, Target_Modification) pair.

    Parameters:
    - df: DataFrame with:
        - Row MultiIndex: ['Experiment', 'Sample_Group']
        - Column MultiIndex: ['Protein_ID', 'Target_Modification']

    Returns:
    - summary_df: DataFrame with:
        - Index: Sample_Group
        - Columns: MultiIndex [(Protein_ID, Target_Modification), Statistic]
    """

    # Group by 'Sample_Group' level in the index (level=1)
    grouped = df.groupby(level='Sample_Group')

    # Aggregate using mean, var, std, and sem (standard error of the mean)
    summary = grouped.agg(['mean', 'var', 'std', ('sem', lambda x: x.std(ddof=1) / np.sqrt(x.count()))])

    # Rename the column levels for clarity
    summary.columns.set_names(['Protein_ID', 'Target_Modification', 'Statistic'], inplace=True)

    return summary

SPEM_and_statistics = summarize_ratios_SP_data(ratios_SP)

# ---------------------------------------------------
# STEP 11 : Calculate Global Extent of Modification (GEM)

#   Calculate Global Extent of Modification (GEM)
#   Summarizes ratios_G by Sample_Group
#   Compute Mean (GEM), Variance, Std_Dev, Std_Err (SEM)
#   Return a MultiIndexed DataFrame

def compute_GEM_and_statistics(flat_df):
    """
    compute_GEM_and_statistics

    Summarizes a flat ratio DataFrame by Sample_Group.

    Parameters:
    - flat_df: DataFrame with columns ['Experiment', 'Sample_Group', <ratio columns>]

    Returns:
    - summary_df: DataFrame indexed by Sample_Group, with MultiIndex columns:
                  (Modification_Ratio, Statistic)
    """
    ratio_cols = [col for col in flat_df.columns if col not in ['Experiment', 'Sample_Group']]

    results = {}

    for ratio_col in ratio_cols:
        grouped = flat_df.groupby('Sample_Group')[ratio_col]

        stats = {
            'Mean': grouped.mean(),
            'Variance': grouped.var(ddof=1),
            'Std_Dev': grouped.std(ddof=1),
            'Std_Err': grouped.sem(ddof=1)
        }

        for stat_name, series in stats.items():
            results[(ratio_col, stat_name)] = series

    summary_df = pd.DataFrame(results)
    summary_df.columns = pd.MultiIndex.from_tuples(summary_df.columns, names=['Modification_Ratio', 'Statistic'])
    summary_df.index.name = 'Sample_Group'

    return summary_df

GEM_and_statistics = compute_GEM_and_statistics(ratios_G)

# ---------------------------------------------------
# STEP 12 : Export to Excel

SPEM_and_statistics.to_excel(os.path.join(output_dir, "SPEM_and_statistics_results.xlsx"))
GEM_and_statistics.to_excel(os.path.join(output_dir,"GEM_and_statistics_results.xlsx"))

print("SPEM and GEM results exported to Excel")

# ---------------------------------------------------
# STEP 13 : Plot results

#   Plots bar charts for each Modification, bars grouped by Sample_Group, with one bar per Protein_ID

def plot_all_modifications_with_error_bars_grouped_by_protein(summary_df, 
                                                              color_map=None, 
                                                              figsize=(10, 6), 
                                                              output_dir=None,
                                                              file_format='png', 
                                                              dpi=300, 
                                                              show=False,
                                                              sample_groups_order=None):
    """
    Plots bar charts for each Modification, showing bars grouped by Sample_Group,
    with one bar per Protein_ID.

    Parameters:
    - summary_df: DataFrame with MultiIndex columns: (Protein_ID, Target_Modification, Statistic)
    - color_map: Optional dict mapping Protein_IDs to specific colors
    - figsize: Tuple for figure size
    - output_dir: Optional path to save plots
    - file_format: 'png', 'pdf', etc.
    - dpi: Image resolution
    - show: Whether to display plots interactively
    - sample_groups_order: Optional list specifying the order of sample groups (e.g., ["G1", "G2", "VEG", "PELB"])
    """

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Get available sample groups from the dataframe
    available_sample_groups = summary_df.index.tolist()
    
    # Determine sample group order
    if sample_groups_order is not None:
        # Use the provided order, but only include groups that exist in the data
        sample_groups_ordered = [g for g in sample_groups_order if g in available_sample_groups]
        # Add any remaining groups that weren't in the custom order (sorted alphabetically)
        remaining_groups = sorted([g for g in available_sample_groups if g not in sample_groups_order])
        sample_groups_ordered.extend(remaining_groups)
        
        # Warn if any requested groups are missing from the data
        missing_groups = [g for g in sample_groups_order if g not in available_sample_groups]
        if missing_groups:
            print(f"Warning: The following sample groups were not found in the data: {missing_groups}")
    else:
        # Default behavior: alphabetical order
        sample_groups_ordered = sorted(available_sample_groups)

    # Extract unique modifications and protein IDs
    mods = summary_df.columns.get_level_values('Target_Modification').unique()
    proteins = summary_df.columns.get_level_values('Protein_ID').unique()

    # Generate color map if not provided
    if color_map is None:
        palette = sns.color_palette("Set3", len(proteins))
        color_map = dict(zip(proteins, palette))

    for mod in mods:
        # Extract mean and sem for this modification
        # Use IndexSlice for better MultiIndex slicing
        idx = pd.IndexSlice
        try:
            mean_df = summary_df.loc[:, idx[:, mod, 'mean']]
            sem_df = summary_df.loc[:, idx[:, mod, 'sem']]
        except KeyError:
            print(f"Skipping {mod} — missing data.")
            continue

        # Drop modification and statistic levels, keep only Protein_ID as columns
        mean_df.columns = mean_df.columns.droplevel([1, 2])
        sem_df.columns = sem_df.columns.droplevel([1, 2])

        # Filter proteins that have data for this modification
        valid_proteins = [p for p in proteins if p in mean_df.columns]

        if not valid_proteins:
            print(f"Skipping {mod} — no valid proteins.")
            continue

        # Reorder rows by sample group order
        mean_df = mean_df.loc[sample_groups_ordered]
        sem_df = sem_df.loc[sample_groups_ordered]

        x = np.arange(len(sample_groups_ordered))
        num_proteins = len(valid_proteins)
        bar_width = 0.8 / num_proteins

        plt.figure(figsize=figsize)

        for i, protein in enumerate(valid_proteins):
            offsets = x + i * bar_width - ((num_proteins - 1) / 2) * bar_width
            bar_vals = mean_df[protein]
            bar_errs = sem_df[protein]

            plt.bar(offsets, bar_vals, width=bar_width, yerr=bar_errs,
                    capsize=4,
                    label=protein,
                    color=color_map.get(protein, 'gray'),
                    edgecolor='black',
                    alpha=0.85)

        plt.xticks(x, sample_groups_ordered, rotation=45, ha='right')
        plt.ylabel('Fraction of modified')
        plt.title(f'Extent of {mod}')
        plt.grid(axis='y', linestyle='--', alpha=0.6)
        plt.legend(title='Protein ID', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()

        def sanitize_filename(name):
            "Replace invalid filename characters with underscores"
            return re.sub(r'[<>:"/\\|?*\n]+', '_', str(name)).strip()

        if output_dir:
            safe_mod = sanitize_filename(mod)
            filename = f"SPEM_{safe_mod}.{file_format}"
            filepath = os.path.join(output_dir, filename)
            plt.savefig(filepath, format=file_format, dpi=dpi, bbox_inches='tight')
            plt.close()
            print(f"Saved: {filepath}")
        elif show:
            plt.show()
        else:
            plt.close()

#   Plots a heatmap of SPEM mean values, with all modifications in the same heatmap

def plot_SPEM_mean_heatmap(SPEM_df, figsize=(10, 6), cmap="viridis", annot=True, save_path=None, show=False, sample_groups_order=None):
    """
    Plots a heatmap of SPEM mean values from a DataFrame with MultiIndex columns:
    (Protein_ID, Target_Modification, Statistic).

    Parameters:
    - SPEM_df: pd.DataFrame — index: Sample_Group, columns: MultiIndex (Protein_ID, Target_Modification, Statistic)
    - figsize: tuple — size of the figure
    - cmap: str — seaborn colormap
    - annot: bool — annotate values
    - save_path: str or Path — optional, if provided saves the plot to this path
    - show: bool — whether to display plot interactively
    - sample_groups_order: Optional list specifying the order of sample groups (e.g., ["G1", "G2", "VEG", "PELB"])

    Returns:
    - None
    """
    # Select only "mean" statistics
    mean_df = SPEM_df.loc[:, SPEM_df.columns.get_level_values("Statistic") == "mean"]

    # Transpose for heatmap: rows = (Protein_ID, Target_Modification), cols = Sample_Group
    mean_df_T = mean_df.T

    # Drop "Statistic" level from index
    if "Statistic" in mean_df_T.index.names:
        mean_df_T.index = mean_df_T.index.droplevel("Statistic")

    # Create a DataFrame with index split for sorting
    idx_df = pd.DataFrame(mean_df_T.index.tolist(), columns=["Protein_ID", "Target_Modification"])

    # Sort by Target_Modification first, then by Protein_ID to keep order consistent
    idx_df = idx_df.sort_values(by=["Target_Modification", "Protein_ID"])

    # Reorder mean_df_T rows accordingly
    mean_df_T = mean_df_T.loc[idx_df.apply(tuple, axis=1)]

    # Create readable row labels (Protein_ID | Target_Modification)
    mean_df_T.index = [f"{pid} | {mod}" for pid, mod in mean_df_T.index]

    # Get available sample groups from the dataframe
    available_sample_groups = mean_df_T.columns.tolist()
    
    # Determine sample group order
    if sample_groups_order is not None:
        # Use the provided order, but only include groups that exist in the data
        sample_groups_ordered = [g for g in sample_groups_order if g in available_sample_groups]
        # Add any remaining groups that weren't in the custom order (sorted alphabetically)
        remaining_groups = sorted([g for g in available_sample_groups if g not in sample_groups_order])
        sample_groups_ordered.extend(remaining_groups)
        
        # Warn if any requested groups are missing from the data
        missing_groups = [g for g in sample_groups_order if g not in available_sample_groups]
        if missing_groups:
            print(f"Warning: The following sample groups were not found in the data: {missing_groups}")
    else:
        # Default behavior: alphabetical order
        sample_groups_ordered = sorted(available_sample_groups)
    
    # Reorder columns according to sample group order
    mean_df_T = mean_df_T[sample_groups_ordered]

    # Plot heatmap
    plt.figure(figsize=figsize)
    sns.heatmap(mean_df_T, cmap=cmap, annot=annot)

    plt.title("SPEM Mean Values Heatmap (Protein_ID | Target_Modification)")
    plt.xlabel("Sample Group")
    plt.ylabel("Protein_ID | Target_Modification")
    plt.tight_layout()

    # Save figure if path is provided
    if save_path:
        filename = "SPEM_heatmap.png"
        filepath = os.path.join(save_path, filename)
        plt.savefig(filepath, dpi=300)
        plt.close()
        print(f"Saved: {filepath}")
    else:
        if show:
            plt.show()
        else:
            plt.close()

def plot_SPEM_mean_heatmaps_by_protein(SPEM_df, 
                                       figsize=(10, 6), 
                                       cmap="viridis", 
                                       annot=True, 
                                       save_path=None,
                                       sample_groups_order=None):
    """
    Plots a separate heatmap for each Protein_ID showing the mean values across sample groups.

    Parameters:
    - SPEM_df: pd.DataFrame — index: Sample_Group, columns: MultiIndex (Protein_ID, Target_Modification, Statistic)
    - figsize: tuple — size of the figure
    - cmap: str — seaborn colormap
    - annot: bool — annotate values
    - save_path: str or Path — optional, directory to save plots (files are named per Protein_ID)
    - sample_groups_order: Optional list specifying the order of sample groups (e.g., ["G1", "G2", "VEG", "PELB"])
    
    Returns:
    - None
    """

    # Extract only the mean values
    mean_df = SPEM_df.loc[:, SPEM_df.columns.get_level_values("Statistic") == "mean"]

    # Drop the 'Statistic' level from the columns
    mean_df.columns = mean_df.columns.droplevel("Statistic")

    # Get list of all unique Protein_IDs
    protein_ids = mean_df.columns.get_level_values("Protein_ID").unique()

    # Get available sample groups from the dataframe
    available_sample_groups = mean_df.index.tolist()
    
    # Determine sample group order
    if sample_groups_order is not None:
        # Use the provided order, but only include groups that exist in the data
        sample_groups_ordered = [g for g in sample_groups_order if g in available_sample_groups]
        # Add any remaining groups that weren't in the custom order (sorted alphabetically)
        remaining_groups = sorted([g for g in available_sample_groups if g not in sample_groups_order])
        sample_groups_ordered.extend(remaining_groups)
        
        # Warn if any requested groups are missing from the data
        missing_groups = [g for g in sample_groups_order if g not in available_sample_groups]
        if missing_groups:
            print(f"Warning: The following sample groups were not found in the data: {missing_groups}")
    else:
        # Default behavior: alphabetical order
        sample_groups_ordered = sorted(available_sample_groups)

    # Reorder the DataFrame index (rows) according to sample_groups_ordered
    mean_df = mean_df.loc[sample_groups_ordered]

    # Create output directory if saving
    if save_path:
        os.makedirs(save_path, exist_ok=True)

    # Loop through each Protein_ID
    for protein in protein_ids:
        # Subset columns for this protein
        protein_df = mean_df.loc[:, mean_df.columns.get_level_values("Protein_ID") == protein]

        # Transpose: rows = Target_Modification, columns = Sample_Group
        protein_df_T = protein_df.T
        protein_df_T.index = protein_df_T.index.get_level_values("Target_Modification")

        # Plot
        plt.figure(figsize=figsize)
        sns.heatmap(protein_df_T, cmap=cmap, annot=annot, fmt=".2f", 
                    cbar_kws={"label": "Mean Target_Modification Extent"})

        plt.title(f"SPEM Mean Values Heatmap for {protein}")
        plt.xlabel("Sample Group")
        plt.ylabel("Target_Modification")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()

        # Save plot if path provided
        if save_path:
            filename = f"{protein}_SPEM_heatmap.png".replace("/", "_").replace("|", "_").replace(" ", "_")
            filepath = os.path.join(save_path, filename)
            plt.savefig(filepath, dpi=300)
            plt.close()
            print(f"Saved: {filepath}")
        else:
            plt.show()

# Plots one heatmap per Target_Modification of SPEM mean values 

def plot_SPEM_mean_heatmap_by_modification(SPEM_df, figsize=(8, 3), cmap="viridis", annot=True, save_path=None, show=False, sample_groups_order=None):
    """
    Plots one heatmap per Target_Modification of SPEM mean values from a DataFrame with MultiIndex columns:
    (Protein_ID, Target_Modification, Statistic).

    Parameters:
    - SPEM_df: pd.DataFrame — index: Sample_Group, columns: MultiIndex (Protein_ID, Target_Modification, Statistic)
    - figsize: tuple — size of the figure
    - cmap: str — seaborn colormap
    - annot: bool — annotate values
    - save_path: str or Path — optional directory to save plots
    - show: bool — whether to display plots interactively
    - sample_groups_order: Optional list specifying the order of sample groups (e.g., ["G1", "G2", "VEG", "PELB"])

    Returns:
    - None
    """
    # Select only "mean" statistics
    mean_df = SPEM_df.loc[:, SPEM_df.columns.get_level_values("Statistic") == "mean"]

    # Get all unique Target_Modifications
    target_mods = mean_df.columns.get_level_values("Target_Modification").unique()

    # Get available sample groups from the dataframe
    available_sample_groups = mean_df.index.tolist()
    
    # Determine sample group order
    if sample_groups_order is not None:
        # Use the provided order, but only include groups that exist in the data
        sample_groups_ordered = [g for g in sample_groups_order if g in available_sample_groups]
        # Add any remaining groups that weren't in the custom order (sorted alphabetically)
        remaining_groups = sorted([g for g in available_sample_groups if g not in sample_groups_order])
        sample_groups_ordered.extend(remaining_groups)
        
        # Warn if any requested groups are missing from the data
        missing_groups = [g for g in sample_groups_order if g not in available_sample_groups]
        if missing_groups:
            print(f"Warning: The following sample groups were not found in the data: {missing_groups}")
    else:
        # Default behavior: alphabetical order
        sample_groups_ordered = sorted(available_sample_groups)

    for mod in target_mods:
        # Filter columns for this Target_Modification (all proteins)
        mod_cols = [col for col in mean_df.columns if col[1] == mod]

        # Subset DataFrame and transpose
        mod_df = mean_df.loc[:, mod_cols].T

        # Index is MultiIndex (Protein_ID, Target_Modification) — drop Target_Modification level
        if "Target_Modification" in mod_df.index.names:
            mod_df.index = mod_df.index.droplevel("Target_Modification")

        # Reorder columns (sample groups) according to sample_groups_ordered
        mod_df = mod_df[sample_groups_ordered]

        # Convert MultiIndex rows to strings for labels
        mod_df.index = [str(idx) if not isinstance(idx, tuple) else " | ".join(map(str, idx)) for idx in mod_df.index]

        # Plot heatmap
        plt.figure(figsize=figsize)
        sns.heatmap(mod_df, cmap=cmap, annot=annot)

        plt.title(f"SPEM Mean Heatmap - Target Modification: {mod}")
        plt.xlabel("Sample Group")
        plt.ylabel("Protein_ID")
        plt.tight_layout()

        # Save or show
        if save_path:
            # Ensure directory exists
            if not os.path.exists(save_path):
                os.makedirs(save_path)

            # Clean mod string for filename
            safe_mod = "".join(c if c.isalnum() else "_" for c in mod)
            filename = f"SPEM_heatmap_{safe_mod}.png"
            filepath = os.path.join(save_path, filename)

            plt.savefig(filepath, dpi=300)
            plt.close()
            print(f"Saved: {filepath}")
        else:
            if show:
                plt.show()
            else:
                plt.close()

print("Plotting GEM results")

#   Plots barplots of GEM modifications grouped by sample group with error bars.

def plot_GEM_modifications_grouped(summary_df, output_dir, color_map=None, default_color='skyblue', sample_groups_order=None):
    """
    Plots barplots of GEM modifications grouped by sample group with error bars.
    
    Parameters:
    - summary_df: pd.DataFrame with MultiIndex columns (Modification_Ratio, Statistic)
    - output_dir: str, directory to save plots
    - color_map: dict, mapping sample group names to colors (optional)
    - default_color: str, fallback color if sample group not in color_map
    - sample_groups_order: Optional list specifying the order of sample groups (e.g., ["G1", "G2", "VEG", "PELB"])
    
    Returns:
    - None
    """
    # Defensive check for columns MultiIndex structure
    expected_levels = ["Modification_Ratio", "Statistic"]
    if not isinstance(summary_df.columns, pd.MultiIndex) or summary_df.columns.names != expected_levels:
        raise ValueError(f"Expected columns with MultiIndex levels {expected_levels}, but got: {summary_df.columns.names}")

    # Get unique modifications from 'Modification_Ratio' level
    modifications = summary_df.columns.get_level_values("Modification_Ratio").unique()

    # Get available sample groups from the dataframe
    available_sample_groups = summary_df.index.tolist()
    
    # Determine sample group order
    if sample_groups_order is not None:
        # Use the provided order, but only include groups that exist in the data
        sample_groups_ordered = [g for g in sample_groups_order if g in available_sample_groups]
        # Add any remaining groups that weren't in the custom order (sorted alphabetically)
        remaining_groups = sorted([g for g in available_sample_groups if g not in sample_groups_order])
        sample_groups_ordered.extend(remaining_groups)
        
        # Warn if any requested groups are missing from the data
        missing_groups = [g for g in sample_groups_order if g not in available_sample_groups]
        if missing_groups:
            print(f"Warning: The following sample groups were not found in the data: {missing_groups}")
    else:
        # Default behavior: alphabetical order
        sample_groups_ordered = sorted(available_sample_groups)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    for mod in modifications:
        plt.figure(figsize=(8, 6))
        
        # Extract mean and standard error columns for the modification
        try:
            mean_vals = summary_df[(mod, "Mean")].reindex(sample_groups_ordered)
            sem_vals = summary_df[(mod, "Std_Err")].reindex(sample_groups_ordered)
        except KeyError as e:
            print(f"Warning: Missing data for modification '{mod}': {e}")
            plt.close()
            continue
        
        # Determine colors for each sample group
        if color_map is not None:
            bar_colors = [color_map.get(sg, default_color) for sg in sample_groups_ordered]
        else:
            bar_colors = default_color
        
        # Plot barplot with error bars and colors
        x_pos = np.arange(len(sample_groups_ordered))
        plt.bar(x_pos, mean_vals, yerr=sem_vals, capsize=5, color=bar_colors, edgecolor='black')
        
        # Labels and title
        plt.xticks(x_pos, sample_groups_ordered, rotation=45, ha='right')
        plt.ylabel("Fraction of modified (Mean ± SEM)")
        plt.title(f"GEM Modification Extent for {mod}")
        plt.tight_layout()

        # Save figure
        safe_mod_name = mod.replace(" ", "_").replace("(", "").replace(")", "").replace("|", "_")
        save_path = os.path.join(output_dir, f"GEM_{safe_mod_name}.png")
        try:
            plt.savefig(save_path, dpi=300)
            print(f"Saved: {save_path}")
        except Exception as e:
            print(f"Error saving {save_path}: {e}")
        finally:
            plt.close()

#   Plots grouped barplots of related GEM modifications by sample group with error bars.

def plot_GEM_modifications_grouped_combined(summary_df, output_dir, combined_labels,
                                            modification_color_map=None, default_colors=None, sample_groups_order=None):
    """
    Plots grouped barplots of related GEM modifications by sample group with error bars.
    
    Parameters:
    - summary_df: pd.DataFrame with MultiIndex columns (Modification_Ratio, Statistic)
    - output_dir: str, directory to save plots
    - combined_labels: list of tuples (list_of_modifications, group_name)
    - modification_color_map: dict mapping modification name to color (optional)
    - default_colors: list of colors to cycle through if no color map provided (optional)
    - sample_groups_order: Optional list specifying the order of sample groups (e.g., ["G1", "G2", "VEG", "PELB"])
    
    Returns:
    - None
    """
    # Defensive check
    expected_levels = ["Modification_Ratio", "Statistic"]
    if not isinstance(summary_df.columns, pd.MultiIndex) or summary_df.columns.names != expected_levels:
        raise ValueError(f"Expected columns with MultiIndex levels {expected_levels}, but got: {summary_df.columns.names}")

    # Get available sample groups from the dataframe
    available_sample_groups = summary_df.index.tolist()
    
    # Determine sample group order
    if sample_groups_order is not None:
        # Use the provided order, but only include groups that exist in the data
        sample_groups_ordered = [g for g in sample_groups_order if g in available_sample_groups]
        # Add any remaining groups that weren't in the custom order (sorted alphabetically)
        remaining_groups = sorted([g for g in available_sample_groups if g not in sample_groups_order])
        sample_groups_ordered.extend(remaining_groups)
        
        # Warn if any requested groups are missing from the data
        missing_groups = [g for g in sample_groups_order if g not in available_sample_groups]
        if missing_groups:
            print(f"Warning: The following sample groups were not found in the data: {missing_groups}")
    else:
        # Default behavior: alphabetical order
        sample_groups_ordered = sorted(available_sample_groups)

    os.makedirs(output_dir, exist_ok=True)

    if default_colors is None:
        default_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']  # matplotlib default palette

    for group_mods, group_name in combined_labels:
        plt.figure(figsize=(10, 6))
        
        means = []
        sems = []
        labels = []

        for mod in group_mods:
            try:
                mean_vals = summary_df[(mod, "Mean")].reindex(sample_groups_ordered)
                sem_vals = summary_df[(mod, "Std_Err")].reindex(sample_groups_ordered)
            except KeyError as e:
                print(f"Warning: Missing data for modification '{mod}': {e}")
                mean_vals = pd.Series([np.nan] * len(sample_groups_ordered), index=sample_groups_ordered)
                sem_vals = pd.Series([np.nan] * len(sample_groups_ordered), index=sample_groups_ordered)

            means.append(mean_vals)
            sems.append(sem_vals)
            labels.append(mod)

        n_mods = len(group_mods)
        x = np.arange(len(sample_groups_ordered))
        width = 0.8 / n_mods  # total bar width divided by number of mods

        for i, (mean_vals, sem_vals, label) in enumerate(zip(means, sems, labels)):
            color = None
            if modification_color_map and label in modification_color_map:
                color = modification_color_map[label]
            else:
                color = default_colors[i % len(default_colors)]

            positions = x - 0.4 + i*width + width/2  # centered bars

            plt.bar(positions, mean_vals, width, yerr=sem_vals, capsize=5,
                    label=label, color=color, edgecolor='black')

        plt.xticks(x, sample_groups_ordered, rotation=45, ha='right')
        plt.ylabel("Fraction of modified (Mean ± SEM)")
        plt.title(f"GEM Modification Extent for {group_name}")
        plt.legend()
        plt.tight_layout()

        safe_group_name = group_name.replace(" ", "_").replace("(", "").replace(")", "").replace("|", "_")
        save_path = os.path.join(output_dir, f"GEM_{safe_group_name}_combined.png")

        try:
            plt.savefig(save_path, dpi=300)
            print(f"Saved: {save_path}")
        except Exception as e:
            print(f"Error saving {save_path}: {e}")
        finally:
            plt.close()

#Plots a heatmap of GEM mean values

def plot_GEM_mean_heatmap(GEM_df, 
                          figsize=(10, 5), 
                          cmap="viridis", 
                          annot=True, 
                          save_path=None,
                          show=False,
                          sample_groups_order=None):
    """
    Plots a heatmap of GEM mean values from a DataFrame with MultiIndex columns:
    (Modification, Statistic), with custom ordering of sample groups.

    Parameters:
    - GEM_df: pd.DataFrame — index: Sample_Group, columns: MultiIndex (Modification, Statistic)
    - figsize: tuple — size of the figure
    - cmap: str — seaborn colormap
    - annot: bool — annotate values
    - save_path: str or Path — optional, if provided saves the plot to this path
    - show: bool — whether to display plot interactively
    - sample_groups_order: Optional list specifying the order of sample groups (e.g., ["G1", "G2", "VEG", "PELB"])

    Returns:
    - None
    """

    # Extract mean values
    mean_df = GEM_df.loc[:, GEM_df.columns.get_level_values("Statistic") == "Mean"]

    # Drop the 'Statistic' level from columns
    mean_df.columns = mean_df.columns.droplevel("Statistic")

    # Transpose so rows = Modification, columns = Sample_Group
    mean_df_T = mean_df.T

    # Get available sample groups from the dataframe
    available_sample_groups = mean_df_T.columns.tolist()
    
    # Determine sample group order
    if sample_groups_order is not None:
        # Use the provided order, but only include groups that exist in the data
        sample_groups_ordered = [g for g in sample_groups_order if g in available_sample_groups]
        # Add any remaining groups that weren't in the custom order (sorted alphabetically)
        remaining_groups = sorted([g for g in available_sample_groups if g not in sample_groups_order])
        sample_groups_ordered.extend(remaining_groups)
        
        # Warn if any requested groups are missing from the data
        missing_groups = [g for g in sample_groups_order if g not in available_sample_groups]
        if missing_groups:
            print(f"Warning: The following sample groups were not found in the data: {missing_groups}")
    else:
        # Default behavior: alphabetical order
        sample_groups_ordered = sorted(available_sample_groups)

    # Reorder columns (sample groups)
    mean_df_T = mean_df_T[sample_groups_ordered]

    # Plot
    plt.figure(figsize=figsize)
    sns.heatmap(mean_df_T, cmap=cmap, annot=annot, fmt=".2f", 
                cbar_kws={"label": "Mean Modification Extent"})

    plt.title("GEM Mean Values Heatmap")
    plt.xlabel("Sample Group")
    plt.ylabel("Modification")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Save figure if path is provided
    if save_path:
        filename = "GEM_heatmap.png"
        filepath = os.path.join(save_path, filename)
        plt.savefig(filepath, dpi=300)
        plt.close()
        print(f"Saved: {filepath}")
    else:
        if show:
            plt.show()
        else:
            plt.close()

plot_all_modifications_with_error_bars_grouped_by_protein(
    SPEM_and_statistics,
    color_map=None, 
    figsize=(10, 6), 
    output_dir=output_dir,
    file_format='png', 
    dpi=300, 
    show=False,
    sample_groups_order = samples_sorted
    )

plot_SPEM_mean_heatmap(
    SPEM_and_statistics,
    figsize=(14, 12),
    cmap="viridis",
    annot=True,
    save_path=output_dir,
    show=False,
    sample_groups_order=samples_sorted
    )

plot_SPEM_mean_heatmaps_by_protein(
    SPEM_and_statistics,
    figsize=(10, 6),
    cmap="viridis",
    annot=True,
    save_path=output_dir,
    sample_groups_order=samples_sorted
    )

plot_SPEM_mean_heatmap_by_modification(
    SPEM_and_statistics,
    save_path=output_dir, 
    show=False, 
    sample_groups_order = samples_sorted
    )

plot_GEM_modifications_grouped(
    GEM_and_statistics,
    output_dir=output_dir,
    color_map = color_map,
    default_color= "skyblue",
    sample_groups_order=samples_sorted
    )

plot_GEM_modifications_grouped_combined(
    summary_df=GEM_and_statistics,
    output_dir=output_dir,
    combined_labels=combined_labels,
    modification_color_map=mod_color_map,
    sample_groups_order=samples_sorted
    )

plot_GEM_mean_heatmap(
    GEM_and_statistics,
    save_path=output_dir,
    sample_groups_order=samples_sorted
    )

print(f"Plots saved in: {output_dir}")