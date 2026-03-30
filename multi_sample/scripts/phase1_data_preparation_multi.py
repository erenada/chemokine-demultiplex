#!/usr/bin/env python3

"""
Phase 1: Multi-File Data Preparation and Quality Assessment
New Chemokine Receptor Analysis Pipeline

Author: Eren Ada, PhD
Date: 10/22/2025

This script performs:
1. Multi-FASTQ file verification and statistics
2. Generation of reverse complement barcode sequences for 21 chemokine targets
3. Quality control analysis across all files
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
        logging.FileHandler('logs/phase1_multi_preparation.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class MultiFileDataPreparation:
    def __init__(self, sample_mapping_file, chemokine_barcodes_file, fastq_dir, output_base_dir=None):
        """
        Initialize multi-file data preparation pipeline
        
        Args:
            sample_mapping_file: Path to sample mapping CSV
            chemokine_barcodes_file: Path to chemokine target barcodes CSV
            fastq_dir: Directory containing FASTQ files
            output_base_dir: Base directory for all outputs (optional)
        """
        self.sample_mapping_file = sample_mapping_file
        self.chemokine_barcodes_file = chemokine_barcodes_file
        self.fastq_dir = fastq_dir
        
        # Set base output directory
        if output_base_dir is None:
            output_base_dir = '.'
        
        # Output directories
        self.output_dirs = {
            'reference_files': os.path.join(output_base_dir, 'reference_files'),
            'qc_reports': os.path.join(output_base_dir, 'qc_reports'),
            'logs': os.path.join(output_base_dir, 'logs')
        }
        
        # Results storage
        self.sample_mapping = None
        self.chemokine_barcodes = None
        self.file_stats = {}
        self.combined_stats = {}
        
        logger.info("Multi-File Data Preparation Pipeline initialized")
        logger.info(f"Sample mapping: {sample_mapping_file}")
        logger.info(f"Chemokine barcodes: {chemokine_barcodes_file}")
        logger.info(f"FASTQ directory: {fastq_dir}")
    
    def load_sample_mapping(self):
        """Load sample mapping file"""
        logger.info("=== Step 1: Loading Sample Mapping ===")
        
        try:
            self.sample_mapping = pd.read_csv(self.sample_mapping_file)
            logger.info(f"Loaded {len(self.sample_mapping)} samples")
            
            for _, row in self.sample_mapping.iterrows():
                logger.info(f"  Sample {row['Sample_no']}: {row['Sample_annotation']} ({row['Condition']} Rep{row['Replicate']})")
            
            return True
        except Exception as e:
            logger.error(f"Error loading sample mapping: {e}")
            return False
    
    def verify_fastq_files(self):
        """Verify all FASTQ files exist and are readable"""
        logger.info("=== Step 2: Verifying FASTQ Files ===")
        
        all_valid = True
        
        for _, row in self.sample_mapping.iterrows():
            fastq_file = os.path.join(self.fastq_dir, row['fastq_file'])
            
            if not os.path.exists(fastq_file):
                logger.error(f"File not found: {fastq_file}")
                all_valid = False
                continue
            
            file_size = os.path.getsize(fastq_file)
            logger.info(f"  {row['fastq_file']} ({file_size:,} bytes)")
            
            # Test FASTQ format
            try:
                fq = pyfastx.Fastq(fastq_file)
                first_read = next(iter(fq))
                logger.info(f"    FASTQ format valid, first read: {first_read.name}")
            except Exception as e:
                logger.error(f"    FASTQ format error: {e}")
                all_valid = False
        
        if all_valid:
            logger.info("All FASTQ files verified successfully")
        else:
            logger.error("Some FASTQ files failed verification")
        
        return all_valid
    
    def analyze_individual_fastq(self, fastq_file, sample_name):
        """Analyze a single FASTQ file"""
        logger.info(f"Analyzing {sample_name}...")
        
        fq = pyfastx.Fastq(fastq_file)
        
        stats = {}
        stats['total_reads'] = len(fq)
        stats['total_bases'] = fq.size
        
        # Read length and quality distribution
        read_lengths = []
        quality_scores = []
        
        sample_size = min(stats['total_reads'], 10000)
        sample_step = max(1, stats['total_reads'] // sample_size)
        
        sampled_reads = 0
        for i, read in enumerate(fq):
            if i % sample_step == 0:
                read_lengths.append(len(read.seq))
                qual_scores = [ord(q) - 33 for q in read.qual]
                quality_scores.append(np.mean(qual_scores))
                sampled_reads += 1
                
                if sampled_reads >= sample_size:
                    break
        
        stats['read_length'] = {
            'min': min(read_lengths) if read_lengths else 0,
            'max': max(read_lengths) if read_lengths else 0,
            'mean': np.mean(read_lengths) if read_lengths else 0,
            'median': np.median(read_lengths) if read_lengths else 0,
            'std': np.std(read_lengths) if read_lengths else 0
        }
        
        stats['quality'] = {
            'mean': np.mean(quality_scores) if quality_scores else 0,
            'median': np.median(quality_scores) if quality_scores else 0,
            'std': np.std(quality_scores) if quality_scores else 0
        }
        
        stats['read_lengths'] = read_lengths
        stats['quality_scores'] = quality_scores
        
        logger.info(f"  Total reads: {stats['total_reads']:,}")
        logger.info(f"  Mean length: {stats['read_length']['mean']:.1f} bp")
        logger.info(f"  Mean quality: {stats['quality']['mean']:.1f}")
        
        return stats
    
    def analyze_all_fastq_files(self):
        """Analyze all FASTQ files"""
        logger.info("=== Step 3: Analyzing All FASTQ Files ===")
        
        for _, row in self.sample_mapping.iterrows():
            fastq_file = os.path.join(self.fastq_dir, row['fastq_file'])
            sample_name = row['Sample_annotation']
            
            stats = self.analyze_individual_fastq(fastq_file, sample_name)
            self.file_stats[sample_name] = stats
        
        # Calculate combined statistics
        total_reads = sum([s['total_reads'] for s in self.file_stats.values()])
        total_bases = sum([s['total_bases'] for s in self.file_stats.values()])
        
        all_lengths = []
        all_qualities = []
        for stats in self.file_stats.values():
            all_lengths.extend(stats['read_lengths'])
            all_qualities.extend(stats['quality_scores'])
        
        self.combined_stats = {
            'total_reads': total_reads,
            'total_bases': total_bases,
            'read_length': {
                'min': min(all_lengths) if all_lengths else 0,
                'max': max(all_lengths) if all_lengths else 0,
                'mean': np.mean(all_lengths) if all_lengths else 0,
                'median': np.median(all_lengths) if all_lengths else 0,
                'std': np.std(all_lengths) if all_lengths else 0
            },
            'quality': {
                'mean': np.mean(all_qualities) if all_qualities else 0,
                'median': np.median(all_qualities) if all_qualities else 0,
                'std': np.std(all_qualities) if all_qualities else 0
            }
        }
        
        logger.info("\nCombined Statistics:")
        logger.info(f"  Total reads across all files: {total_reads:,}")
        logger.info(f"  Total bases: {total_bases:,}")
        logger.info(f"  Overall mean read length: {self.combined_stats['read_length']['mean']:.1f} bp")
        logger.info(f"  Overall mean quality: {self.combined_stats['quality']['mean']:.1f}")
    
    def load_chemokine_barcodes(self):
        """Load and validate chemokine barcode file"""
        logger.info("=== Step 4: Loading Chemokine Barcodes ===")
        
        try:
            # Load CSV without assuming header, then assign column names
            # encoding='utf-8-sig' handles BOM if present
            self.chemokine_barcodes = pd.read_csv(self.chemokine_barcodes_file, header=None, encoding='utf-8-sig')
            self.chemokine_barcodes.columns = ['Chemokine', 'Barcode']
            
            logger.info(f"Loaded {len(self.chemokine_barcodes)} chemokine target barcodes")
            
            # Display target list
            targets = list(self.chemokine_barcodes['Chemokine'])
            logger.info(f"Targets: {', '.join(targets)}")
            
            # Validate barcode sequences
            barcode_lengths = self.chemokine_barcodes['Barcode'].str.len()
            logger.info(f"Barcode lengths: {barcode_lengths.unique()}")
            
            # Check for valid characters
            valid_chars = set('ATCG')
            for i, row in self.chemokine_barcodes.iterrows():
                barcode = row['Barcode']
                if not set(barcode).issubset(valid_chars):
                    logger.warning(f"Invalid characters in barcode {row['Chemokine']}: {barcode}")
            
            logger.info("Barcode validation completed")
            return True
            
        except Exception as e:
            logger.error(f"Error loading chemokine barcodes: {e}")
            return False
    
    def generate_reverse_complements(self):
        """Generate reverse complement sequences for chemokine barcodes"""
        logger.info("=== Step 5: Generating Reverse Complement Barcodes ===")
        
        chemokine_rc = self.chemokine_barcodes.copy()
        
        chemokine_rc['Barcode_rc'] = chemokine_rc['Barcode'].apply(
            lambda x: str(Seq(x).reverse_complement())
        )
        
        # Save to reference file
        output_file = os.path.join(self.output_dirs['reference_files'], 'chemokine_barcodes_with_rc.csv')
        chemokine_rc.to_csv(output_file, index=False)
        
        logger.info(f"Saved chemokine barcodes with RC to: {output_file}")
        
        # Log examples
        logger.info("Examples of reverse complement generation:")
        for i in range(min(5, len(chemokine_rc))):
            name = chemokine_rc.iloc[i]['Chemokine']
            orig = chemokine_rc.iloc[i]['Barcode']
            rc = chemokine_rc.iloc[i]['Barcode_rc']
            logger.info(f"  {name}: {orig} -> {rc}")
    
    def create_qc_plots(self):
        """Create comprehensive quality control plots"""
        logger.info("=== Step 6: Creating QC Plots ===")
        
        plt.style.use('default')
        sns.set_palette("husl")
        
        # Create figure with multiple subplots
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Phase 1: Multi-File Data Quality Assessment', fontsize=16, fontweight='bold')
        
        # Plot 1: Per-file read counts
        samples = list(self.file_stats.keys())
        read_counts = [self.file_stats[s]['total_reads'] for s in samples]
        
        axes[0, 0].bar(range(len(samples)), read_counts, alpha=0.7, edgecolor='black')
        axes[0, 0].set_xlabel('Sample')
        axes[0, 0].set_ylabel('Read Count')
        axes[0, 0].set_title('Reads per Sample')
        axes[0, 0].set_xticks(range(len(samples)))
        axes[0, 0].set_xticklabels([s.split('-')[-1] for s in samples], rotation=45, ha='right')
        
        # Plot 2: Read length distributions by file
        for sample, stats in self.file_stats.items():
            axes[0, 1].hist(stats['read_lengths'], bins=30, alpha=0.5, label=sample.split('-')[-1])
        axes[0, 1].set_xlabel('Read Length (bp)')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Read Length Distributions')
        axes[0, 1].legend(fontsize=8)
        
        # Plot 3: Quality score distributions by file
        for sample, stats in self.file_stats.items():
            axes[0, 2].hist(stats['quality_scores'], bins=30, alpha=0.5, label=sample.split('-')[-1])
        axes[0, 2].set_xlabel('Average Quality Score')
        axes[0, 2].set_ylabel('Frequency')
        axes[0, 2].set_title('Quality Score Distributions')
        axes[0, 2].legend(fontsize=8)
        
        # Plot 4: Combined dataset summary
        axes[1, 0].axis('off')
        summary_text = f"""
Combined Dataset Summary:

Total Files: {len(self.file_stats)}
Total Reads: {self.combined_stats['total_reads']:,}
Total Bases: {self.combined_stats['total_bases']:,}

Read Length:
  Mean: {self.combined_stats['read_length']['mean']:.1f} bp
  Range: {self.combined_stats['read_length']['min']}-{self.combined_stats['read_length']['max']} bp

Quality:
  Mean Score: {self.combined_stats['quality']['mean']:.1f}

Chemokine Targets: {len(self.chemokine_barcodes)}

Analysis Date: {datetime.now().strftime('%m/%d/%Y')}
        """
        
        axes[1, 0].text(0.05, 0.95, summary_text, transform=axes[1, 0].transAxes,
                        fontsize=10, verticalalignment='top', fontfamily='monospace',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.5))
        
        # Plot 5: Per-sample statistics table
        axes[1, 1].axis('off')
        
        table_data = []
        for sample in samples:
            stats = self.file_stats[sample]
            table_data.append([
                sample.split('-')[-1],
                f"{stats['total_reads']:,}",
                f"{stats['read_length']['mean']:.0f}",
                f"{stats['quality']['mean']:.1f}"
            ])
        
        table = axes[1, 1].table(cellText=table_data,
                                colLabels=['Sample', 'Reads', 'Avg Len', 'Avg Qual'],
                                cellLoc='center',
                                loc='center',
                                bbox=[0, 0.3, 1, 0.6])
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2)
        axes[1, 1].set_title('Per-Sample Statistics', pad=20)
        
        # Plot 6: Experimental design
        axes[1, 2].axis('off')
        design_text = f"""
Experimental Design:

Conditions:
  - GFP: 3 replicates
  - Negative: 3 replicates

Chemokine Targets: {len(self.chemokine_barcodes)}
  {', '.join(list(self.chemokine_barcodes['Chemokine'])[:10])}
  ... and {len(self.chemokine_barcodes) - 10} more

Barcode Specifications:
  Length: 20 base pairs
  All reverse complements generated
  
Pipeline Status:
  Phase 1: Complete
  Next: Phase 3 Quantification
        """
        
        axes[1, 2].text(0.05, 0.95, design_text, transform=axes[1, 2].transAxes,
                        fontsize=9, verticalalignment='top', fontfamily='monospace',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.3))
        
        plt.tight_layout()
        
        # Save plot
        plot_file = os.path.join(self.output_dirs['qc_reports'], 'phase1_multi_file_assessment.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        logger.info(f"QC plots saved to: {plot_file}")
        
        plt.close()
    
    def save_summary_report(self):
        """Save comprehensive summary report"""
        logger.info("=== Step 7: Generating Summary Report ===")
        
        # JSON summary
        summary = {
            'pipeline_info': {
                'phase': 'Phase 1 - Multi-File Data Preparation',
                'date': datetime.now().isoformat(),
                'total_files': len(self.file_stats)
            },
            'combined_stats': self.combined_stats,
            'per_file_stats': {k: {
                'total_reads': v['total_reads'],
                'total_bases': v['total_bases'],
                'mean_length': v['read_length']['mean'],
                'mean_quality': v['quality']['mean']
            } for k, v in self.file_stats.items()},
            'target_info': {
                'total_targets': len(self.chemokine_barcodes),
                'target_list': list(self.chemokine_barcodes['Chemokine'])
            }
        }
        
        json_file = os.path.join(self.output_dirs['qc_reports'], 'phase1_summary.json')
        with open(json_file, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        
        logger.info(f"Summary report saved to: {json_file}")
        
        # Human-readable report
        report_file = os.path.join(self.output_dirs['qc_reports'], 'phase1_combined_report.txt')
        with open(report_file, 'w') as f:
            f.write("NEW CHEMOKINE RECEPTOR ANALYSIS PIPELINE\n")
            f.write("Phase 1: Multi-File Data Preparation and Quality Assessment\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Analysis Date: {datetime.now().strftime('%m/%d/%Y')}\n")
            f.write(f"Total Files Processed: {len(self.file_stats)}\n\n")
            
            f.write("COMBINED DATASET OVERVIEW\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total Reads: {self.combined_stats['total_reads']:,}\n")
            f.write(f"Total Bases: {self.combined_stats['total_bases']:,}\n")
            f.write(f"Average Read Length: {self.combined_stats['read_length']['mean']:.1f} bp\n")
            f.write(f"Read Length Range: {self.combined_stats['read_length']['min']}-{self.combined_stats['read_length']['max']} bp\n")
            f.write(f"Average Quality Score: {self.combined_stats['quality']['mean']:.1f}\n\n")
            
            f.write("PER-FILE STATISTICS\n")
            f.write("-" * 40 + "\n")
            for sample, stats in self.file_stats.items():
                f.write(f"\n{sample}:\n")
                f.write(f"  Reads: {stats['total_reads']:,}\n")
                f.write(f"  Mean Length: {stats['read_length']['mean']:.1f} bp\n")
                f.write(f"  Mean Quality: {stats['quality']['mean']:.1f}\n")
            
            f.write(f"\nCHEMOKINE TARGETS\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total Targets: {len(self.chemokine_barcodes)}\n")
            f.write("Target List:\n")
            for target in self.chemokine_barcodes['Chemokine']:
                f.write(f"  - {target}\n")
            
            f.write("\nFILES GENERATED\n")
            f.write("-" * 40 + "\n")
            f.write("- reference_files/chemokine_barcodes_with_rc.csv\n")
            f.write("- qc_reports/phase1_multi_file_assessment.png\n")
            f.write("- qc_reports/phase1_summary.json\n")
            f.write("- logs/phase1_multi_preparation.log\n\n")
            
            f.write("NEXT STEPS\n")
            f.write("-" * 40 + "\n")
            f.write("1. Review quality control results\n")
            f.write("2. Proceed to Phase 3: Barcode Quantification\n")
            f.write("   (Phase 2 skipped - files pre-demultiplexed)\n")
        
        logger.info(f"Human-readable report saved to: {report_file}")
    
    def run_pipeline(self):
        """Execute the complete Phase 1 pipeline"""
        logger.info("Starting Phase 1: Multi-File Data Preparation")
        
        # Create output directories
        for dir_name in self.output_dirs.values():
            os.makedirs(dir_name, exist_ok=True)
        
        try:
            # Execute pipeline steps
            if not self.load_sample_mapping():
                return False
            
            if not self.verify_fastq_files():
                return False
            
            self.analyze_all_fastq_files()
            
            if not self.load_chemokine_barcodes():
                return False
            
            self.generate_reverse_complements()
            self.create_qc_plots()
            self.save_summary_report()
            
            logger.info("=" * 70)
            logger.info("Phase 1 completed successfully!")
            logger.info("=" * 70)
            logger.info("Ready to proceed to Phase 3: Barcode Quantification")
            
            return True
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False

def main():
    """Main execution function"""
    print("New Chemokine Receptor Analysis Pipeline")
    print("Phase 1: Multi-File Data Preparation and Quality Assessment")
    print("=" * 70)
    
    # File paths
    sample_mapping = "New_analysis/reference_files/sample_mapping.csv"
    chemokine_barcodes = "New_analysis/new_data/new_barcodes.csv"  # 2-column input file
    fastq_dir = "New_analysis/data_nov14/SYFFKS_fastq"
    output_base_dir = "New_analysis/data_nov14"
    
    # Check if files exist
    if not os.path.exists(sample_mapping):
        print(f"Error: File not found - {sample_mapping}")
        sys.exit(1)
    
    if not os.path.exists(chemokine_barcodes):
        print(f"Error: File not found - {chemokine_barcodes}")
        sys.exit(1)
    
    if not os.path.exists(fastq_dir):
        print(f"Error: Directory not found - {fastq_dir}")
        sys.exit(1)
    
    # Initialize and run pipeline
    pipeline = MultiFileDataPreparation(sample_mapping, chemokine_barcodes, fastq_dir, output_base_dir)
    success = pipeline.run_pipeline()
    
    if success:
        print("\nPhase 1 completed successfully.")
        print("\nGenerated files:")
        print(f"- {output_base_dir}/reference_files/chemokine_barcodes_with_rc.csv")
        print(f"- {output_base_dir}/qc_reports/phase1_multi_file_assessment.png")
        print(f"- {output_base_dir}/qc_reports/phase1_combined_report.txt")
        print(f"- {output_base_dir}/qc_reports/phase1_summary.json")
        print(f"- {output_base_dir}/logs/phase1_multi_preparation.log")
        print("\nReady to proceed to Phase 3: Barcode Quantification")
    else:
        print("\nPhase 1 failed. Check logs for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()

