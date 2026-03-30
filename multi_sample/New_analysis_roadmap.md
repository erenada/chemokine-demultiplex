# New Chemokine Receptor Analysis Pipeline Roadmap

**Author:** Eren Ada, PhD  
**Date:** 10/22/2025  
**Project:** ARGES Tumor Chemokine Receptor Analysis - New Dataset Processing

## Project Overview

This analysis processes a new dataset of chemokine receptor screening data with:
- 6 pre-separated FASTQ files (3 GFP+ replicates, 3 Negative control replicates)
- 21 chemokine receptor targets (expanded from original 19)
- ARGES tumor model experimental design
- Files already organized by condition (no sample demultiplexing required)

## New Dataset Characteristics

### Input Files Structure
```
New_analysis/
├── new_data/
│   ├── new_barcodes.csv (21 chemokine targets)
│   └── new_Fastq/
│       ├── CKR-ARGES-Tumor-GFP_1.fastq
│       ├── CKR-ARGES-Tumor-GFP_2.fastq  
│       ├── CKR-ARGES-Tumor-GFP_3.fastq
│       ├── CKR-ARGES-Tumor-Neg_1.fastq
│       ├── CKR-ARGES-Tumor-Neg_2.fastq
│       └── CKR-ARGES-Tumor-Neg_3.fastq
```

### Key Differences from Original Dataset
- **Sample Organization:** Files pre-separated by condition (GFP vs Negative)
- **Target Expansion:** 21 chemokine targets (added Gpr183, tTA)
- **Experimental Model:** ARGES tumor samples vs original screening model
- **Replication:** 3 biological replicates per condition

## Modified Pipeline Workflow

### Phase 0: Setup and Environment Preparation
- [ ] **0.1** Create analysis directory structure
  - [ ] Set up scripts/, reference_files/, results/, qc_reports/, logs/ directories
  - [ ] Copy and adapt pipeline scripts for new data structure
  - [ ] Create sample annotation mapping file

- [ ] **0.2** Environment activation
  - [ ] Activate demultiplex_pipeline conda environment
  - [ ] Verify all required tools and libraries are available

### Phase 1: Data Preparation and Quality Assessment
- [ ] **1.1** Multi-file data verification
  - [ ] Verify integrity of all 6 FASTQ files
  - [ ] Validate new barcode file format (21 targets)
  - [ ] Generate read statistics for each file individually
  - [ ] Create combined statistics summary

- [ ] **1.2** Generate reference files
  - [ ] Create reverse complement sequences for 21 chemokine barcodes
  - [ ] Generate sample mapping file for 6 input files
  - [ ] Create searchable barcode databases

- [ ] **1.3** Quality control assessment
  - [ ] Run FastQC on all 6 input files
  - [ ] Generate read length distributions per file
  - [ ] Assess overall data quality metrics
  - [ ] Compare file characteristics across replicates

**Expected Outputs:**
- `reference_files/chemokine_barcodes_with_rc.csv` (21 targets)
- `reference_files/sample_mapping.csv` (6 files mapping)
- `qc_reports/phase1_multi_file_assessment.png`
- Individual and combined quality reports

### Phase 2: SKIPPED - Sample Demultiplexing
**Reason:** Input files are already separated by experimental condition

### Phase 3: Chemokine Barcode Quantification (Modified)
- [ ] **3.1** Prepare expanded target barcode library
  - [ ] Process 21 chemokine barcodes (vs original 19)
  - [ ] Generate reverse complement versions for all targets
  - [ ] Create comprehensive searchable database
  - [ ] Validate new targets: Gpr183, tTA

- [ ] **3.2** Multi-file barcode quantification
  - [ ] Process each of 6 FASTQ files independently
  - [ ] Search for all 21 chemokine barcodes in each file
  - [ ] Handle bidirectional reads (forward and reverse complement)
  - [ ] Apply fuzzy matching with appropriate mismatch tolerance
  - [ ] Generate per-file barcode count matrices

- [ ] **3.3** Generate combined count matrix
  - [ ] Create 6×21 samples × chemokines count matrix
  - [ ] Calculate per-sample detection rates
  - [ ] Apply normalization methods (CPM, TPM, Log2)
  - [ ] Export analysis-ready formats

**Expected Outputs:**
- `results/raw_count_matrix_6x21.csv`
- `results/cpm_normalized_counts.csv`
- `results/tpm_normalized_counts.csv` 
- `results/log2_cpm_normalized_counts.csv`
- Per-sample barcode detection summaries

### Phase 4: Quality Control and Validation (Enhanced)
- [ ] **4.1** Multi-file pipeline metrics
  - [ ] Calculate per-file barcode detection efficiency
  - [ ] Assess detection consistency across replicates
  - [ ] Compare GFP vs Negative detection patterns
  - [ ] Identify potential batch effects

- [ ] **4.2** Biological validation
  - [ ] Compare GFP+ vs Negative samples within replicates
  - [ ] Assess replicate consistency (correlation analysis)
  - [ ] Identify condition-specific chemokine patterns
  - [ ] Flag outliers or technical issues

- [ ] **4.3** Enhanced QC reporting
  - [ ] Multi-file processing summary statistics
  - [ ] Cross-replicate correlation matrices
  - [ ] Condition comparison visualizations
  - [ ] Technical validation metrics

### Phase 5: Output Generation and Analysis-Ready Data
- [ ] **5.1** Primary analysis outputs
  - [ ] Final 6×21 count matrix with proper sample annotations
  - [ ] Multiple normalization versions for different analyses
  - [ ] Sample metadata and experimental design file
  - [ ] Comprehensive detection statistics

- [ ] **5.2** Analysis-ready datasets
  - [ ] Format for statistical comparison (GFP vs Negative)
  - [ ] Include proper experimental metadata
  - [ ] Generate R/Python-compatible data objects
  - [ ] Prepare for downstream differential analysis

## Required Script Modifications

### Scripts to Adapt from Original Pipeline
1. **`phase1_data_preparation.py`**
   - Modify to handle 6 input files instead of 1
   - Update for 21 chemokine targets
   - Add multi-file quality assessment

2. **`phase3_barcode_quantification.py`**
   - Adapt for direct FASTQ file processing (skip demultiplexing)
   - Update for 21 target barcodes
   - Add per-file processing logic

3. **`phase4_quality_control.py`**
   - Enhance for multi-file analysis
   - Add replicate consistency checks
   - Include condition comparison metrics

### New Files to Create
1. **`sample_mapping.csv`**
   ```
   Sample_no,Sample_annotation,Condition,Replicate,fastq_file
   1,ARGES-Tumor-GFP-Rep1,GFP,1,CKR-ARGES-Tumor-GFP_1.fastq
   2,ARGES-Tumor-GFP-Rep2,GFP,2,CKR-ARGES-Tumor-GFP_2.fastq
   3,ARGES-Tumor-GFP-Rep3,GFP,3,CKR-ARGES-Tumor-GFP_3.fastq
   4,ARGES-Tumor-Neg-Rep1,Negative,1,CKR-ARGES-Tumor-Neg_1.fastq
   5,ARGES-Tumor-Neg-Rep2,Negative,2,CKR-ARGES-Tumor-Neg_2.fastq
   6,ARGES-Tumor-Neg-Rep3,Negative,3,CKR-ARGES-Tumor-Neg_3.fastq
   ```

## Expected Final Directory Structure

```
New_analysis/
├── New_analysis_roadmap.md
├── new_data/
│   ├── new_barcodes.csv
│   └── new_Fastq/ (6 FASTQ files)
├── reference_files/
│   ├── chemokine_barcodes_with_rc.csv (21 targets)
│   └── sample_mapping.csv
├── scripts/
│   ├── phase1_data_preparation_multi.py
│   ├── phase3_barcode_quantification_multi.py
│   └── phase4_quality_control_enhanced.py
├── results/
│   ├── raw_count_matrix_6x21.csv
│   ├── cpm_normalized_counts.csv
│   ├── tpm_normalized_counts.csv
│   ├── log2_cpm_normalized_counts.csv
│   └── detection_summary_by_sample.csv
├── qc_reports/
│   ├── phase1_multi_file_assessment.png
│   ├── phase1_combined_report.txt
│   ├── phase3_detection_efficiency.png
│   ├── phase4_replicate_correlations.png
│   └── fastqc/ (6 individual FastQC reports)
└── logs/
    ├── phase1_multi_preparation.log
    ├── phase3_multi_quantification.log
    └── phase4_enhanced_qc.log
```

## Success Criteria

- [ ] All 6 FASTQ files successfully processed
- [ ] >90% chemokine barcode detection across all 21 targets
- [ ] High replicate consistency within conditions (r > 0.8)
- [ ] Clear differentiation between GFP and Negative conditions
- [ ] Comprehensive quality control documentation
- [ ] Analysis-ready 6×21 count matrix generated

## Key Advantages of New Dataset Structure

1. **Pre-separated samples:** Eliminates demultiplexing step and associated errors
2. **Biological replicates:** Enables proper statistical analysis
3. **Expanded target set:** Additional chemokine receptors for comprehensive analysis
4. **Controlled conditions:** Clear GFP+ vs Negative comparison

## Timeline Estimate

- **Phase 0-1:** 1 day (setup and data preparation)
- **Phase 3:** 1-2 days (barcode quantification across 6 files)
- **Phase 4-5:** 1 day (QC and final outputs)
- **Total:** 3-4 days

## Next Steps

1. Set up directory structure and copy/adapt required scripts
2. Create sample mapping file for the 6 input files
3. Begin Phase 1 with multi-file data preparation
4. Proceed through modified pipeline workflow

## Notes

- Pipeline optimized for pre-separated FASTQ files
- Enhanced for multi-replicate experimental design
- Expanded target set requires updated barcode databases
- Quality control enhanced for biological replicate analysis
- Ready for downstream statistical comparison between conditions
