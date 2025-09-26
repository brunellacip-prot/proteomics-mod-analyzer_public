# proteomics-mod-analyzer
Script Name: proteomics-mod-analyzer.py

## Overview

This Python script analyzes MaxQuant evidence output files to compute protein modification extents within groups of sample replicates. It calculates two key metrics:

- **SPEM (Single-Protein Extent of Modification)**: Modification extent calculated for each protein separately
- **GEM (Global Extent of Modification)**: Modification extent calculated considering all proteins collectively

## Prerequisites

### Required Libraries
```python
import pandas as pd
import os
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
```

### Input File Requirements
- **Format**: Excel (.xlsx) or tab-delimited text (.txt)
- **Required columns**: 
  - `Experiment`
  - `Proteins` 
  - `Modifications`
  - `MS/MS count`
  - `Sequence`
  - `Modified sequence`

## Configuration

### Step 1: Set File Paths
```python
input_file = r"path/to/your/evidence_file.xlsx"
output_dir = r"path/to/your/output/directory"
```

### Step 2: Define Analysis Parameters

#### Protein Selection
```python
# Specify UniProt IDs or leave empty to analyze all proteins
uniprot_id = ["P02668", "P02662", "P02666", "P02663"]  # Leave empty [] for all proteins
```

#### Target Amino Acids and Modifications
```python
target_aa = ["N", "Q", "M", "S", "T", "N|Q", "S|T"]
target_modification = ["Deamidation (NQ)", "Oxidation (M)", "Phospho (ST)"]
separated_modified_aa = ["N(Deamidation (NQ))", "Q(Deamidation (NQ))", 
                        "S(Phospho (ST))", "T(Phospho (ST))"]
```

#### Combined Modifications
```python
combined_labels = [
    (["N(Deamidation (NQ))", "Q(Deamidation (NQ))"], "N|Q(Deamidation (NQ))"),
    (["S(Phospho (ST))", "T(Phospho (ST))"], "S|T(Phospho (ST))")
]
```

## Script Workflow

### Steps 1-4: Data Loading and Organization
1. **Load evidence file** - Supports both Excel and text formats
2. **Extract UniProt IDs** - Automatically parses protein identifiers
3. **Define parameters** - Set target amino acids and modifications
4. **Organize experiments** - Groups experiments by sample (removes _REP suffix)

### Steps 5-6: Single-Protein Analysis (SP)
5. **Calculate targetable sites (SP)** - Counts potential modification sites per protein
6. **Calculate modified sites (SP)** - Counts actual modifications per protein

### Steps 7-8: Global Analysis (G)
7. **Calculate targetable sites (G)** - Sums across all proteins per experiment
8. **Calculate modified sites (G)** - Sums modified sites across all proteins

### Steps 9-11: Extent Calculations
9. **Calculate G ratios** - Modified/targetable sites globally
10. **Calculate SPEM** - Single-protein modification extents with statistics
11. **Calculate GEM** - Global modification extents with statistics

### Steps 12-13: Output and Visualization
12. **Export results** - Saves SPEM and GEM results to Excel files
13. **Generate plots** - Creates bar charts and heatmaps

## Key Functions

### Data Processing Functions
- `load_evidence_file()` - Loads MaxQuant evidence files
- `extract_uniprot_ids_from_proteins()` - Extracts UniProt IDs from protein column
- `filter_and_order_by_sample_group_protein()` - Filters data by sample groups and proteins

### Analysis Functions
- `compute_targetable_sites_SP()` - Calculates targetable amino acid sites
- `create_weighted_modification_matrix_with_combined()` - Creates modification count matrices
- `modified_targetable_SP_ratio()` - Computes modification ratios
- `summarize_ratios_SP_data()` - Calculates statistics for SPEM
- `compute_GEM_and_statistics()` - Calculates statistics for GEM

### Visualization Functions
- `plot_all_modifications_with_error_bars()` - Creates bar plots with error bars
- `plot_all_protein_modifications()` - Plots individual protein modifications
- `plot_SPEM_mean_heatmap()` - Creates SPEM heatmaps
- `plot_GEM_mean_heatmap()` - Creates GEM heatmaps

## Output Files

### Excel Files
- `SPEM_and_statistics_results.xlsx` - Single-protein modification extents
- `GEM_and_statistics_results.xlsx` - Global modification extents

### Plot Files
- Bar charts for each modification type
- Individual protein modification plots
- Heatmaps showing modification patterns across samples

## Data Structure

### Sample Naming Convention
The script expects experiment names with replicate suffixes:
- `SAMPLE_REP1`, `SAMPLE_REP2`, `SAMPLE_REP3`, etc.
- Automatically groups by removing `_REP*` suffix

### Statistical Outputs
For both SPEM and GEM calculations:
- **Mean**: Average modification extent
- **Variance**: Variance across replicates
- **Std_Dev**: Standard deviation
- **Std_Err (SEM)**: Standard error of the mean

## Customization Options

### Color Mapping for Plots
```python
colors = {
    'AZZURRITE': "#e2c0cdc5",
    'CAS': "#dca8c9",
    'CINABRO': "#960B0BD7",
    'CaCO3': "#93afb5",
    'MINIO': "#d64908"
}
```

### Plot Parameters
- Figure size, DPI, file format
- Color schemes for different sample groups
- Statistical error bar display options

## Common Modification Labels

The script supports various modification types commonly found in MaxQuant:
- `Deamidation (NQ)` - Deamidation of asparagine and glutamine
- `Oxidation (M)` - Methionine oxidation
- `Phospho (ST)` - Serine/threonine phosphorylation
- `Phospho (STY)` - Serine/threonine/tyrosine phosphorylation

## Usage Notes

1. Ensure your evidence file contains all required columns
2. Experiment names should follow the `SAMPLE_REP#` convention
3. The script automatically handles missing data and zero values
4. Output directory will be created if it doesn't exist
5. All plots can be saved in various formats (PNG, PDF, etc.)

## Troubleshooting

### Common Issues
- **File format errors**: Ensure file is .xlsx or .txt with proper delimiters
- **Missing columns**: Verify all required columns are present in input file
- **UniProt ID extraction**: Check if protein identifiers follow expected format
- **Memory issues**: Large datasets may require increased system memory

### Data Quality Checks
- The script includes warnings for missing experiments or proteins
- Zero values are handled appropriately in ratio calculations
- Statistical calculations account for missing replicates
