# Chemokine Receptor Barcode Demultiplexing Pipeline

**Author:** Eren Ada, PhD  
**Last updated:** 03/30/2026

---

## Overview

This repository contains two Python pipelines for processing and quantifying chemokine receptor barcodes from Illumina sequencing data. Both pipelines use approximate barcode matching via [edlib](https://github.com/Martinsos/edlib) to handle sequencing errors, and support both forward and reverse-complement barcode orientations.

| Pipeline | Location | Use case |
|---|---|---|
| Main pipeline (4-phase) | `scripts/` | Single pooled FASTQ — demultiplexes samples, then quantifies barcodes |
| Multi-sample pipeline | `multi_sample/` | Per-sample FASTQ files already separated — skips demultiplexing, quantifies directly |

---

## Repository Structure

```
chemokine-demultiplex/
├── scripts/                              # Main 4-phase pipeline
│   ├── phase1_data_preparation.py
│   ├── phase2_sample_demultiplexing.py
│   ├── phase3_barcode_quantification.py
│   └── phase4_quality_control.py
├── multi_sample/                         # Multi-sample pipeline
│   ├── README.md                         # Detailed usage instructions
│   └── scripts/
│       ├── phase1_data_preparation_multi.py
│       └── phase3_barcode_quantification_multi.py
├── environment.yml                       # Conda environment
├── setup_environment.sh                  # Environment setup helper
└── .gitignore
```

---

## Setup

```bash
conda env create -f environment.yml
conda activate demultiplex_pipeline
```

---

## Which Pipeline to Use

**Use `scripts/` (main pipeline) if:**
- All samples are pooled in a single FASTQ file.
- The pipeline needs to demultiplex reads into per-sample files first (Phase 2), then quantify barcodes.

**Use `multi_sample/` if:**
- Each sample already has its own FASTQ file (e.g., standard Illumina index-based demultiplexing was done upstream).
- You only need barcode quantification across samples.

See [`multi_sample/README.md`](multi_sample/README.md) for detailed setup, input file formats, and expected outputs for the multi-sample pipeline.

---

## Dependencies

- Python 3.9+
- pandas, numpy
- matplotlib, seaborn
- biopython
- pyfastx
- edlib
- tqdm

All dependencies are defined in `environment.yml`.
