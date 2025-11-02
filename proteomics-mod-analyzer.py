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

input_file = "evidence_test.txt"
output_dir = "results_directory"

#------------------------------------------------------------#
# STEP 2 : Define key input parameters

#   target_aa : provided as str in letter code ("N", "Q", "S", "N|Q", etc...)
#       if you are targeting two amino acids at the same time please use "|" (e.g. "N|Q")
#   target_modification : provided as str (exactly as reported in the "Modifications" column)
#   separated_modified_aa : provided as str (exactly as reported in the "Modified sequence" column)
#       e.g. "N(Deamidation (NQ))", "Q(Deamidation (NQ))", "S(Phospho (ST))", "T(Phospho (ST)), etc..."

config_file = "Input_modifications_MQ_test.xlsx"

# # Fill with labels of modifications to be combined, grouping them according to the example.
combined_labels = [
    (["N(Deamidation (NQ))", "Q(Deamidation (NQ))"], "N|Q(Deamidation (NQ))")] #add correct labels to combine

#   uniprot_id : uniprot IDs (e.g. P02668, P02662, etc...) provided as list
uniprot_id = ["P02662", "P02663", "P02666", "P02668"]  # Leave empty [] to extract all

#------------------------------------------------------------#
# STEP 3 : Load the MQ evidence output file (Excel or TXT)

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

# The end