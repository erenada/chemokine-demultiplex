#!/usr/bin/env python3

"""
Phase 1: Data Preparation and Quality Assessment
Chemokine Receptor Demultiplexing Pipeline

Author: Eren Ada, PhD
Date: 08/07/2025

This script performs:
1. FASTQ file verification and basic statistics
2. Generation of reverse complement barcode sequences
3. Initial quality control analysis
4. Preparation of reference files
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Seq import Seq
from Bio import SeqIO
import pyfastx
import subprocess
import logging
from datetime import datetime
import json

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/phase1_preparation.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class DataPreparation:
    def __init__(self, fastq_file, demux_barcodes_file, chemokine_barcodes_file):
        """
        Initialize data preparation pipeline
        
        Args:
            fastq_file: Path to input FASTQ file
            demux_barcodes_file: Path to sample demultiplexing barcodes CSV
            chemokine_barcodes_file: Path to chemokine target barcodes CSV
        """
        self.fastq_file = fastq_file
        self.demux_barcodes_file = demux_barcodes_file
        self.chemokine_barcodes_file = chemokine_barcodes_file
        
        # Output directories
        self.output_dirs = {
            'reference_files': 'reference_files',
            'qc_reports': 'qc_reports',
            'logs': 'logs'
        }
        
        # Results storage
        self.stats = {}
        self.demux_barcodes = None
        self.chemokine_barcodes = None
        
        logger.info("Data Preparation Pipeline initialized")
        logger.info(f"Input FASTQ: {fastq_file}")
        logger.info(f"Demux barcodes: {demux_barcodes_file}")
        logger.info(f"Chemokine barcodes: {chemokine_barcodes_file}")
    
    def verify_files(self):
        """Verify input files exist and are readable"""
        logger.info("=== Step 1: File Verification ===")
        
        files_to_check = [
            self.fastq_file,
            self.demux_barcodes_file,
            self.chemokine_barcodes_file
        ]
        
        for file_path in files_to_check:
            if not os.path.exists(file_path):
                logger.error(f"File not found: {file_path}")
                sys.exit(1)
            
            file_size = os.path.getsize(file_path)
            logger.info(f"✓ {file_path} ({file_size:,} bytes)")
        
        # Check FASTQ format
        try:
            # Test reading first few records
            fq = pyfastx.Fastq(self.fastq_file)
            first_read = next(iter(fq))
            logger.info(f"✓ FASTQ format valid, first read: {first_read.name}")
        except Exception as e:
            logger.error(f"FASTQ format error: {e}")
            sys.exit(1)
        
        logger.info("All input files verified successfully")
    
    def analyze_fastq_stats(self):
        """Generate comprehensive FASTQ statistics"""
        logger.info("=== Step 2: FASTQ Analysis ===")
        
        # Basic statistics using pyfastx
        fq = pyfastx.Fastq(self.fastq_file)
        
        self.stats['total_reads'] = len(fq)
        self.stats['total_bases'] = fq.size
        
        logger.info(f"Total reads: {self.stats['total_reads']:,}")
        logger.info(f"Total bases: {self.stats['total_bases']:,}")
        
        # Read length distribution
        logger.info("Analyzing read length distribution...")
        read_lengths = []
        quality_scores = []
        
        # Sample reads for analysis (use all if <100k reads, otherwise sample)
        sample_size = min(self.stats['total_reads'], 50000)
        sample_step = max(1, self.stats['total_reads'] // sample_size)
        
        sampled_reads = 0
        for i, read in enumerate(fq):
            if i % sample_step == 0:
                read_lengths.append(len(read.seq))
                # Calculate average quality score for this read
                qual_scores = [ord(q) - 33 for q in read.qual]  # Convert to numeric
                quality_scores.append(np.mean(qual_scores))
                sampled_reads += 1
                
                if sampled_reads >= sample_size:
                    break
        
        # Calculate statistics
        self.stats['read_length'] = {
            'min': min(read_lengths),
            'max': max(read_lengths),
            'mean': np.mean(read_lengths),
            'median': np.median(read_lengths),
            'std': np.std(read_lengths)
        }
        
        self.stats['quality'] = {
            'mean': np.mean(quality_scores),
            'median': np.median(quality_scores),
            'std': np.std(quality_scores)
        }
        
        self.stats['read_lengths'] = read_lengths
        self.stats['quality_scores'] = quality_scores
        
        logger.info(f"Read length - Min: {self.stats['read_length']['min']}, "
                   f"Max: {self.stats['read_length']['max']}, "
                   f"Mean: {self.stats['read_length']['mean']:.1f}")
        logger.info(f"Quality score - Mean: {self.stats['quality']['mean']:.1f}")
        
        # Run FastQC for detailed quality analysis
        self.run_fastqc()
    
    def run_fastqc(self):
        """Run FastQC quality control"""
        logger.info("Running FastQC analysis...")
        
        fastqc_output = os.path.join(self.output_dirs['qc_reports'], 'fastqc')
        os.makedirs(fastqc_output, exist_ok=True)
        
        cmd = [
            'fastqc',
            self.fastq_file,
            '-o', fastqc_output,
            '--quiet'
        ]
        
        try:
            subprocess.run(cmd, check=True)
            logger.info(f"✓ FastQC completed, results in {fastqc_output}/")
        except subprocess.CalledProcessError as e:
            logger.warning(f"FastQC failed: {e}")
    
    def load_barcode_files(self):
        """Load and validate barcode files"""
        logger.info("=== Step 3: Loading Barcode Files ===")
        
        # Load demultiplexing barcodes
        try:
            self.demux_barcodes = pd.read_csv(self.demux_barcodes_file)
            logger.info(f"✓ Loaded {len(self.demux_barcodes)} sample demux barcodes")
            logger.info(f"Samples: {list(self.demux_barcodes['Sample_annotation'])}")
        except Exception as e:
            logger.error(f"Error loading demux barcodes: {e}")
            sys.exit(1)
        
        # Load chemokine barcodes
        try:
            self.chemokine_barcodes = pd.read_csv(self.chemokine_barcodes_file)
            logger.info(f"✓ Loaded {len(self.chemokine_barcodes)} chemokine target barcodes")
            logger.info(f"Targets: {list(self.chemokine_barcodes['Chemokine'])}")
        except Exception as e:
            logger.error(f"Error loading chemokine barcodes: {e}")
            sys.exit(1)
        
        # Validate barcode sequences
        self.validate_barcodes()
    
    def validate_barcodes(self):
        """Validate barcode sequences"""
        logger.info("Validating barcode sequences...")
        
        # Check demux barcodes
        demux_lengths = self.demux_barcodes['demultiplex_barcode'].str.len()
        logger.info(f"Demux barcode lengths: {demux_lengths.unique()}")
        
        # Check chemokine barcodes
        chemokine_lengths = self.chemokine_barcodes['Barcode'].str.len()
        logger.info(f"Chemokine barcode lengths: {chemokine_lengths.unique()}")
        
        # Check for invalid characters
        valid_chars = set('ATCG')
        for i, row in self.demux_barcodes.iterrows():
            barcode = row['demultiplex_barcode']
            if not set(barcode).issubset(valid_chars):
                logger.warning(f"Invalid characters in demux barcode {row['Sample_annotation']}: {barcode}")
        
        for i, row in self.chemokine_barcodes.iterrows():
            barcode = row['Barcode']
            if not set(barcode).issubset(valid_chars):
                logger.warning(f"Invalid characters in chemokine barcode {row['Chemokine']}: {barcode}")
        
        logger.info("Barcode validation completed")
    
    def generate_reverse_complements(self):
        """Generate reverse complement sequences for all barcodes"""
        logger.info("=== Step 4: Generating Reverse Complement Barcodes ===")
        
        # Generate RC for demux barcodes
        demux_rc = self.demux_barcodes.copy()
        demux_rc['demultiplex_barcode_rc'] = demux_rc['demultiplex_barcode'].apply(
            lambda x: str(Seq(x).reverse_complement())
        )
        
        # Generate RC for chemokine barcodes
        chemokine_rc = self.chemokine_barcodes.copy()
        chemokine_rc['Barcode_rc'] = chemokine_rc['Barcode'].apply(
            lambda x: str(Seq(x).reverse_complement())
        )
        
        # Save to reference files
        demux_rc_file = os.path.join(self.output_dirs['reference_files'], 'sample_barcodes_with_rc.csv')
        chemokine_rc_file = os.path.join(self.output_dirs['reference_files'], 'chemokine_barcodes_with_rc.csv')
        
        demux_rc.to_csv(demux_rc_file, index=False)
        chemokine_rc.to_csv(chemokine_rc_file, index=False)
        
        logger.info(f"✓ Saved demux barcodes with RC to: {demux_rc_file}")
        logger.info(f"✓ Saved chemokine barcodes with RC to: {chemokine_rc_file}")
        
        # Log some examples
        logger.info("Examples of reverse complement generation:")
        for i in range(min(3, len(demux_rc))):
            orig = demux_rc.iloc[i]['demultiplex_barcode']
            rc = demux_rc.iloc[i]['demultiplex_barcode_rc']
            sample = demux_rc.iloc[i]['Sample_annotation']
            logger.info(f"  {sample}: {orig} -> {rc}")
    
    def create_qc_plots(self):
        """Create quality control plots"""
        logger.info("=== Step 5: Creating QC Plots ===")
        
        # Set style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Phase 1: Data Quality Assessment', fontsize=16, fontweight='bold')
        
        # Plot 1: Read length distribution
        axes[0, 0].hist(self.stats['read_lengths'], bins=50, alpha=0.7, edgecolor='black')
        axes[0, 0].set_xlabel('Read Length (bp)')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Read Length Distribution')
        axes[0, 0].axvline(self.stats['read_length']['mean'], color='red', linestyle='--', 
                          label=f"Mean: {self.stats['read_length']['mean']:.1f} bp")
        axes[0, 0].legend()
        
        # Plot 2: Quality score distribution
        axes[0, 1].hist(self.stats['quality_scores'], bins=50, alpha=0.7, edgecolor='black')
        axes[0, 1].set_xlabel('Average Quality Score')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Quality Score Distribution')
        axes[0, 1].axvline(self.stats['quality']['mean'], color='red', linestyle='--',
                          label=f"Mean: {self.stats['quality']['mean']:.1f}")
        axes[0, 1].legend()
        
        # Plot 3: Dataset summary (expanded to fill available space)
        axes[1, 0].axis('off')
        summary_text = f"""
Dataset Summary:

Total Reads: {self.stats['total_reads']:,}
Total Bases: {self.stats['total_bases']:,}

Read Length:
  Mean: {self.stats['read_length']['mean']:.1f} bp
  Range: {self.stats['read_length']['min']}-{self.stats['read_length']['max']} bp

Quality:
  Mean Score: {self.stats['quality']['mean']:.1f}

Samples: {len(self.demux_barcodes)} (8bp barcodes)
Chemokine Targets: {len(self.chemokine_barcodes)} (20bp barcodes)

Sample Barcodes:
{chr(10).join([f"  {row['Sample_annotation']}: {row['demultiplex_barcode']}" for _, row in self.demux_barcodes.head(6).iterrows()])}

Processing Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}
        """
        
        axes[1, 0].text(0.05, 0.95, summary_text, transform=axes[1, 0].transAxes,
                        fontsize=11, verticalalignment='top', fontfamily='monospace',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.5))
        
        # Plot 4: Dataset summary (right panel)
        axes[1, 1].axis('off')
        
        # Add barcode design information
        design_text = f"""
Barcode Design Specifications:

Sample Demultiplexing:
  Length: 8 base pairs
  Count: {len(self.demux_barcodes)} unique barcodes
  
Chemokine Targeting:
  Length: 20 base pairs  
  Count: {len(self.chemokine_barcodes)} unique barcodes

Expected Read Structure:
  [Sample BC][Target BC][Other sequences]
  
Quality Metrics:
  All barcodes validated for:
  - Correct length
  - Valid nucleotides (A,T,G,C)
  - Reverse complements generated
        """
        
        axes[1, 1].text(0.05, 0.95, design_text, transform=axes[1, 1].transAxes,
                        fontsize=11, verticalalignment='top', fontfamily='monospace',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.3))
        
        plt.tight_layout()
        
        # Save plot
        plot_file = os.path.join(self.output_dirs['qc_reports'], 'phase1_quality_assessment.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        logger.info(f"✓ QC plots saved to: {plot_file}")
        
        plt.close()
    
    def save_summary_report(self):
        """Save comprehensive summary report"""
        logger.info("=== Step 6: Generating Summary Report ===")
        
        # Prepare summary data
        summary = {
            'pipeline_info': {
                'phase': 'Phase 1 - Data Preparation',
                'date': datetime.now().isoformat(),
                'input_file': self.fastq_file
            },
            'dataset_stats': self.stats,
            'sample_info': {
                'total_samples': len(self.demux_barcodes),
                'sample_list': self.demux_barcodes[['Sample_no', 'Sample_annotation', 'demultiplex_barcode']].to_dict('records')
            },
            'target_info': {
                'total_targets': len(self.chemokine_barcodes),
                'target_list': self.chemokine_barcodes[['Chemokine', 'Barcode']].to_dict('records')
            }
        }
        
        # Save JSON summary
        json_file = os.path.join(self.output_dirs['qc_reports'], 'phase1_summary.json')
        with open(json_file, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        
        logger.info(f"✓ Summary report saved to: {json_file}")
        
        # Create human-readable report
        report_file = os.path.join(self.output_dirs['qc_reports'], 'phase1_report.txt')
        with open(report_file, 'w') as f:
            f.write("CHEMOKINE RECEPTOR DEMULTIPLEXING PIPELINE\n")
            f.write("Phase 1: Data Preparation and Quality Assessment\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Input File: {self.fastq_file}\n\n")
            
            f.write("DATASET OVERVIEW\n")
            f.write("-" * 30 + "\n")
            f.write(f"Total Reads: {self.stats['total_reads']:,}\n")
            f.write(f"Total Bases: {self.stats['total_bases']:,}\n")
            f.write(f"Average Read Length: {self.stats['read_length']['mean']:.1f} bp\n")
            f.write(f"Read Length Range: {self.stats['read_length']['min']}-{self.stats['read_length']['max']} bp\n")
            f.write(f"Average Quality Score: {self.stats['quality']['mean']:.1f}\n\n")
            
            f.write("EXPERIMENTAL DESIGN\n")
            f.write("-" * 30 + "\n")
            f.write(f"Total Samples: {len(self.demux_barcodes)}\n")
            f.write("Sample Details:\n")
            for _, row in self.demux_barcodes.iterrows():
                f.write(f"  {row['Sample_no']}: {row['Sample_annotation']} ({row['demultiplex_barcode']})\n")
            
            f.write(f"\nTotal Chemokine Targets: {len(self.chemokine_barcodes)}\n")
            f.write("Target List: " + ", ".join(self.chemokine_barcodes['Chemokine']) + "\n\n")
            
            f.write("FILES GENERATED\n")
            f.write("-" * 30 + "\n")
            f.write("- reference_files/sample_barcodes_with_rc.csv\n")
            f.write("- reference_files/chemokine_barcodes_with_rc.csv\n")
            f.write("- qc_reports/phase1_quality_assessment.png\n")
            f.write("- qc_reports/fastqc/ (FastQC results)\n")
            f.write("- logs/phase1_preparation.log\n\n")
            
            f.write("NEXT STEPS\n")
            f.write("-" * 30 + "\n")
            f.write("1. Review quality control results\n")
            f.write("2. Proceed to Phase 2: Sample Demultiplexing\n")
            f.write("3. Use generated reference files for barcode matching\n")
        
        logger.info(f"✓ Human-readable report saved to: {report_file}")
    
    def run_pipeline(self):
        """Execute the complete Phase 1 pipeline"""
        logger.info("Starting Phase 1: Data Preparation and Quality Assessment")
        
        # Create output directories
        for dir_name in self.output_dirs.values():
            os.makedirs(dir_name, exist_ok=True)
        
        try:
            # Execute pipeline steps
            self.verify_files()
            self.analyze_fastq_stats()
            self.load_barcode_files()
            self.generate_reverse_complements()
            self.create_qc_plots()
            self.save_summary_report()
            
            logger.info("=" * 60)
            logger.info("Phase 1 completed successfully!")
            logger.info("=" * 60)
            logger.info("Ready to proceed to Phase 2: Sample Demultiplexing")
            
            return True
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            return False

def main():
    """Main execution function"""
    print("Chemokine Receptor Demultiplexing Pipeline")
    print("Phase 1: Data Preparation and Quality Assessment")
    print("=" * 60)
    
    # File paths
    fastq_file = "FCHKLT_1_EL1.fastq"
    demux_barcodes = "demultiplex_barcode.csv"
    chemokine_barcodes = "chemokine_barcodes.csv"
    
    # Check if files exist
    for file in [fastq_file, demux_barcodes, chemokine_barcodes]:
        if not os.path.exists(file):
            print(f"Error: File not found - {file}")
            print("Please ensure all required files are in the current directory.")
            sys.exit(1)
    
    # Initialize and run pipeline
    pipeline = DataPreparation(fastq_file, demux_barcodes, chemokine_barcodes)
    success = pipeline.run_pipeline()
    
    if success:
        print("\nPhase 1 completed successfully.")
        print("\nGenerated files:")
        print("- reference_files/sample_barcodes_with_rc.csv")
        print("- reference_files/chemokine_barcodes_with_rc.csv") 
        print("- qc_reports/phase1_quality_assessment.png")
        print("- qc_reports/phase1_report.txt")
        print("- logs/phase1_preparation.log")
        print("\nReady to proceed to Phase 2: Sample Demultiplexing")
    else:
        print("\nPhase 1 failed. Check logs for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()