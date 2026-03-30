#!/bin/bash

# Chemokine Receptor Demultiplexing Pipeline - Environment Setup
# Author: Eren Ada, PhD
# Date: 08/07/2025

set -e  # Exit on any error

echo "=================================================="
echo "Setting up Chemokine Demultiplexing Pipeline Environment"
echo "=================================================="

# Environment name
ENV_NAME="demultiplex_pipeline"

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda is not installed or not in PATH"
    echo "Please install Miniconda or Anaconda first"
    exit 1
fi

echo "✓ Conda found: $(conda --version)"

# Remove existing environment if it exists
echo "Checking for existing environment..."
if conda env list | grep -q "^${ENV_NAME}"; then
    echo "Environment '${ENV_NAME}' already exists. Removing..."
    conda env remove -n ${ENV_NAME} -y
fi

echo "Creating new conda environment: ${ENV_NAME}"
conda create -n ${ENV_NAME} python=3.9 -y

echo "Activating environment..."
source activate ${ENV_NAME}

echo "=================================================="
echo "Installing Bioinformatics Tools"
echo "=================================================="

# Install bioinformatics tools from bioconda
echo "Installing bioinformatics tools from bioconda..."
conda install -c bioconda -c conda-forge -y \
    seqkit \
    cutadapt \
    bioawk \
    samtools \
    fastqc \
    gawk

echo "=================================================="
echo "Installing Python Packages"
echo "=================================================="

# Install Python packages
echo "Installing Python packages..."
pip install \
    biopython \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    scipy \
    plotly \
    pyfastx \
    edlib \
    regex \
    tqdm \
    jupyter

echo "=================================================="
echo "Verifying Installation"
echo "=================================================="

echo "Checking installed tools..."

# Check bioinformatics tools
tools=("seqkit" "cutadapt" "bioawk" "samtools" "fastqc")
for tool in "${tools[@]}"; do
    if command -v ${tool} &> /dev/null; then
        version=$(${tool} --version 2>&1 | head -n1 || echo "version unknown")
        echo "✓ ${tool}: ${version}"
    else
        echo "✗ ${tool}: NOT FOUND"
    fi
done

# Check Python packages
echo "Checking Python packages..."
python -c "
import sys
packages = [
    'Bio', 'pandas', 'numpy', 'matplotlib', 'seaborn', 
    'scipy', 'plotly', 'pyfastx', 'edlib', 'regex', 'tqdm'
]

for pkg in packages:
    try:
        __import__(pkg)
        print(f'✓ {pkg}: Available')
    except ImportError:
        print(f'✗ {pkg}: NOT FOUND')
"

echo "=================================================="
echo "Creating Project Directory Structure"
echo "=================================================="

# Create directory structure
mkdir -p scripts
mkdir -p reference_files
mkdir -p demultiplexed_samples
mkdir -p results
mkdir -p qc_reports
mkdir -p logs

echo "Created directory structure:"
echo "  - scripts/          (custom Python scripts)"
echo "  - reference_files/  (barcode references)"
echo "  - demultiplexed_samples/ (sample-specific FASTQ files)"
echo "  - results/          (count matrices and analysis results)"
echo "  - qc_reports/       (quality control reports)"
echo "  - logs/             (pipeline logs)"

echo "=================================================="
echo "Creating Environment Activation Script"
echo "=================================================="

cat > activate_pipeline.sh << 'EOF'
#!/bin/bash
# Quick activation script for the demultiplexing pipeline environment

echo "Activating demultiplex_pipeline environment..."
source activate demultiplex_pipeline

echo "Environment activated!"
echo "Available tools:"
echo "  - seqkit (FASTQ manipulation)"
echo "  - cutadapt (adapter trimming)"
echo "  - fastqc (quality control)"
echo "  - bioawk (sequence processing)"
echo "  - samtools (sequence file manipulation)"
echo ""
echo "Python packages:"
echo "  - biopython, pandas, numpy, matplotlib"
echo "  - pyfastx (fast FASTQ parsing)"
echo "  - edlib (sequence alignment)"
echo ""
echo "Ready to run demultiplexing pipeline!"
EOF

chmod +x activate_pipeline.sh

echo "=================================================="
echo "Installation Summary"
echo "=================================================="

echo "Environment setup complete!"
echo ""
echo "To activate the environment:"
echo "  conda activate ${ENV_NAME}"
echo "  OR"
echo "  ./activate_pipeline.sh"
echo ""
echo "Next steps:"
echo "  1. Run 'conda activate ${ENV_NAME}'"
echo "  2. Proceed with Phase 1 of the pipeline"
echo ""
echo "Environment location: $(conda env list | grep ${ENV_NAME} | awk '{print $2}')"

# Deactivate environment
conda deactivate

echo "=================================================="
echo "Setup Complete!"
echo "=================================================="