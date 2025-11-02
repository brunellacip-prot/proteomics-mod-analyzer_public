# MaxQuant Modification Extent Analysis Guide

## proteomics-mod-analyzer

Script Name: proteomics-mod-analyzer.py

## Overview

This script analyzes MaxQuant evidence files to calculate modification extents in proteomics data. It computes two key metrics:

- **SPEM (Single-Protein Extent of Modification)**: Modification extent calculated separately for each protein
- **GEM (Global Extent of Modification)**: Modification extent calculated across all proteins collectively

---

## Requirements

### Python Libraries

```python
pandas
numpy
matplotlib
seaborn
openpyxl  # For Excel file handling
```

Install with:

```bash
pip install pandas numpy matplotlib seaborn openpyxl
```

### Input Files

1. **Evidence File** (`.txt` or `.xlsx`)
   - MaxQuant output file
   - Must contain columns: `Experiment`, `Proteins`, `Modifications`, `MS/MS count`, `Modified sequence`, `Sequence`

2. **Configuration File** (`Input_modifications_MQ.xlsx`)
   - Excel file defining target modifications
   - Required columns:
     - `target_aa`: Target amino acids (e.g., "N", "Q", "S|T")
     - `target_modification`: Modification names from MaxQuant
     - `separated_modified_aa`: Modified amino acid format (e.g., "N(Deamidation (NQ))")

---

## Setup Instructions

### Step 1: Define File Paths

```python
input_file = "evidence_test.txt"  # Your MaxQuant evidence file
output_dir = "results_directory"   # Output folder for results
config_file = "Input_modifications_MQ_test.xlsx"  # Configuration file
```

### Step 2: Configure Analysis Parameters

#### Protein Selection

```python
# Option 1: Analyze specific proteins
uniprot_id = ["P02662", "P02663", "P02666", "P02668"]

# Option 2: Analyze all proteins (leave empty)
uniprot_id = []
```

#### Combined Modifications

Define modifications that should be grouped together:

```python
combined_labels = [
    (["N(Deamidation (NQ))", "Q(Deamidation (NQ))"], "N|Q(Deamidation (NQ))")
]
```

**Format**: `(list_of_individual_modifications, "combined_label")`

### Step 3: Experiment Naming Convention

Experiments must follow this naming pattern:

```
SAMPLE_NAME_REP1
SAMPLE_NAME_REP2
SAMPLE_NAME_REP3
```

**Example**:
- `CAS_REP1`, `CAS_REP2`, `CAS_REP3` → Sample group: `CAS`
- `WHEY_REP1`, `WHEY_REP2` → Sample group: `WHEY`

---

## Configuration File Format

Create an Excel file with these columns:

| target_aa | target_modification | separated_modified_aa |
|-----------|--------------------|-----------------------|
| N\|Q | Deamidation (NQ) | N(Deamidation (NQ)) |
| N\|Q | Deamidation (NQ) | Q(Deamidation (NQ)) |
| S\|T | Phospho (ST) | S(Phospho (ST)) |
| S\|T | Phospho (ST) | T(Phospho (ST)) |

**Notes**:
- Use `|` to separate multiple amino acids (e.g., `N|Q`, `S|T`)
- Modification names must exactly match MaxQuant output
- Modified AA format must match the `Modified sequence` column format

---

## Understanding the Analysis

### Workflow Overview

1. **Load Data**: Reads MaxQuant evidence file
2. **Group Experiments**: Organizes replicates by sample name
3. **Filter Data**: Extracts data for specified proteins
4. **Calculate Targetable Sites**: Counts amino acids that can be modified
5. **Calculate Modified Sites**: Counts observed modifications
6. **Compute Ratios**: Calculates modification extents
7. **Statistical Summary**: Computes mean, variance, std dev, and SEM

### Key Concepts

#### Targetable Sites (SP)

Number of amino acids in a protein sequence that could potentially be modified, weighted by MS/MS counts:

```
Targetable Sites = Σ(#amino_acids × MS/MS_count)
```

#### Modified Sites (SP)

Number of observed modifications, weighted by MS/MS counts:

```
Modified Sites = Σ(#modifications × MS/MS_count)
```

#### SPEM (Single-Protein Extent)

Modification extent for each protein individually:

```
SPEM = Modified Sites (SP) / Targetable Sites (SP)
```

#### GEM (Global Extent)

Modification extent across all proteins combined:

```
GEM = Σ Modified Sites (all proteins) / Σ Targetable Sites (all proteins)
```

---

## Output Files

The script generates two Excel files in the output directory:

### 1. `SPEM_and_statistics_results.xlsx`

**Structure**:
- **Index**: Sample groups
- **Columns**: Multi-level index
  - Level 1: Protein ID
  - Level 2: Target modification
  - Level 3: Statistic (mean, var, std, sem)

**Example**:

| Sample_Group | P02662 - N\|Q(Deam) - Mean | P02662 - N\|Q(Deam) - Std_Err | ... |
|--------------|---------------------------|------------------------------|-----|
| CAS | 0.234 | 0.012 | ... |
| WHEY | 0.189 | 0.015 | ... |

### 2. `GEM_and_statistics_results.xlsx`

**Structure**:
- **Index**: Sample groups
- **Columns**: Multi-level index
  - Level 1: Modification ratio
  - Level 2: Statistic (Mean, Variance, Std_Dev, Std_Err)

**Example**:

| Sample_Group | N\|Q - Mean | N\|Q - Std_Err | S\|T - Mean | ... |
|--------------|------------|---------------|------------|-----|
| CAS | 0.210 | 0.008 | 0.156 | ... |
| WHEY | 0.178 | 0.011 | 0.142 | ... |

---

## Common Issues & Troubleshooting

### Issue: "No data found for experiment(s)"

**Solution**: Check that experiment names in the evidence file match the naming convention (`SAMPLE_REP#`)

### Issue: "No columns found matching target_aa"

**Solution**: Verify that `target_aa` values in the config file match amino acids in your data

### Issue: Missing modifications in results

**Solution**: Ensure modification names in config file exactly match MaxQuant's `Modifications` column

### Issue: Empty SPEM/GEM results

**Solution**: 
- Check that UniProt IDs are correct
- Verify protein naming format in the `Proteins` column
- Ensure modifications exist in your data

### Issue: Division by zero warnings

**Solution**: Normal behavior when targetable sites = 0; these are handled as NaN values

---

## Best Practices

1. **Backup Your Data**: Keep original evidence files unchanged
2. **Test First**: Run with a small subset of proteins initially
3. **Check Config File**: Verify all modification labels match your data exactly
4. **Review Intermediate Outputs**: Check the printed messages during execution
5. **Validate Results**: Spot-check a few calculations manually

---

## Example Use Case

**Scenario**: Analyzing deamidation in milk proteins

```python
# Files
input_file = "milk_proteins_evidence.txt"
output_dir = "deamidation_results"
config_file = "deamidation_config.xlsx"

# Target proteins
uniprot_id = ["P02662", "P02663", "P02666", "P02668"]  # αs1, αs2, β, κ-casein

# Combined modifications
combined_labels = [
    (["N(Deamidation (NQ))", "Q(Deamidation (NQ))"], "N|Q(Deamidation (NQ))")
]
```

This will analyze deamidation extent in four casein proteins across your sample groups.

---

## Quick Start Checklist

- [ ] Install required Python libraries
- [ ] Prepare MaxQuant evidence file (`.txt` or `.xlsx`)
- [ ] Create configuration Excel file with modification details
- [ ] Set file paths in script
- [ ] Define target proteins (or leave empty for all)
- [ ] Configure combined modifications if needed
- [ ] Verify experiment naming follows `SAMPLE_REP#` convention
- [ ] Create output directory
- [ ] Run script and check console output
- [ ] Review generated Excel files

---

## Support

For issues or questions:
1. Check the printed messages during script execution
2. Verify input file formats match requirements
3. Review the configuration file structure
4. Check that experiment naming follows conventions

---

## Version Notes

- Supports both `.xlsx` and `.txt` evidence files
- Handles multiple protein ID formats (`sp|P02662|CASB_BOVIN` or `P02662`)
- Automatically extracts all proteins if `uniprot_id = []`
- Handles missing data with appropriate NaN values
- Console feedback with 🔹 indicators for progress tracking