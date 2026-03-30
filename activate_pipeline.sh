#!/bin/bash
# Quick activation script for the demultiplexing pipeline environment
# Author: Eren Ada, PhD
# Date: 08/07/2025

echo "🧬 Activating Chemokine Demultiplexing Pipeline Environment"
echo "============================================================"

# Activate the conda environment
source /Users/eren/miniconda3/bin/activate demultiplex_pipeline

echo "Environment: demultiplex_pipeline"
echo ""
echo "Available Bioinformatics Tools:"
echo "  ✓ cutadapt v5.1 - Adapter trimming and barcode extraction"
echo "  ✓ samtools v1.21 - Sequence file manipulation"
echo "  ✓ fastqc v0.12.1 - Quality control reports"
echo "  ✓ seqkit v2.10.0 - FASTQ manipulation and statistics"
echo ""
echo "Python Libraries:"
echo "  ✓ biopython v1.85 - Sequence analysis"
echo "  ✓ pyfastx v2.2.0 - Fast FASTQ parsing"
echo "  ✓ pandas v2.3.1 - Data manipulation"
echo "  ✓ numpy v2.0.2 - Numerical computing"
echo "  ✓ edlib - Sequence alignment for fuzzy matching"
echo ""
echo "Project Structure:"
echo "  📁 scripts/ - Custom analysis scripts"
echo "  📁 reference_files/ - Barcode reference files"
echo "  📁 demultiplexed_samples/ - Sample-specific FASTQ files"
echo "  📁 results/ - Count matrices and analysis outputs"
echo "  📁 qc_reports/ - Quality control reports"
echo "  📁 logs/ - Pipeline execution logs"
echo ""
echo "🚀 Ready to run the demultiplexing pipeline!"
echo ""
echo "Next steps:"
echo "  1. Proceed with Phase 1: Data preparation"
echo "  2. Generate reverse complement barcodes"
echo "  3. Start sample demultiplexing"