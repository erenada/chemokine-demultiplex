#!/usr/bin/env python3

"""
Phase 2: Sample Demultiplexing
Chemokine Receptor Demultiplexing Pipeline

Author: Eren Ada, PhD
Date: 08/07/2025

This script performs:
1. Load sample barcodes with reverse complements
2. Search for barcodes in reads (both orientations)
3. Assign reads to samples with mismatch tolerance
4. Generate sample-specific FASTQ files
5. Quality control and statistics reporting
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Seq import Seq
import pyfastx
import edlib
import regex
import logging
from datetime import datetime
import json
from collections import defaultdict, Counter
from tqdm import tqdm

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/phase2_demultiplexing.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class SampleDemultiplexer:
    def __init__(self, fastq_file, barcode_file, max_mismatches=1, min_quality=20):
        """
        Initialize sample demultiplexing pipeline
        
        Args:
            fastq_file: Path to input FASTQ file
            barcode_file: Path to sample barcodes with reverse complements CSV
            max_mismatches: Maximum mismatches allowed in barcode matching
            min_quality: Minimum average quality score for reads to process
        """
        self.fastq_file = fastq_file
        self.barcode_file = barcode_file
        self.max_mismatches = max_mismatches
        self.min_quality = min_quality
        
        # Output directories
        self.output_dirs = {
            'demultiplexed_samples': 'demultiplexed_samples',
            'qc_reports': 'qc_reports',
            'logs': 'logs'
        }
        
        # Results storage
        self.barcodes = None
        self.barcode_patterns = {}
        self.assignment_stats = defaultdict(int)
        self.quality_stats = defaultdict(list)
        self.mismatch_stats = defaultdict(int)
        
        logger.info("Sample Demultiplexer initialized")
        logger.info(f"Input FASTQ: {fastq_file}")
        logger.info(f"Barcode file: {barcode_file}")
        logger.info(f"Max mismatches: {max_mismatches}")
        logger.info(f"Min quality: {min_quality}")
    
    def load_barcodes(self):
        """Load sample barcodes with reverse complements"""
        logger.info("=== Step 1: Loading Sample Barcodes ===")
        
        try:
            self.barcodes = pd.read_csv(self.barcode_file)
            logger.info(f"✓ Loaded {len(self.barcodes)} sample barcodes")
            
            # Prepare barcode patterns for matching
            for _, row in self.barcodes.iterrows():
                sample_id = row['Sample_annotation']
                forward_bc = row['demultiplex_barcode']
                reverse_bc = row['demultiplex_barcode_rc']
                
                self.barcode_patterns[sample_id] = {
                    'forward': forward_bc,
                    'reverse': reverse_bc,
                    'sample_no': row['Sample_no']
                }
                
                logger.info(f"  {sample_id}: {forward_bc} / {reverse_bc}")
            
        except Exception as e:
            logger.error(f"Error loading barcodes: {e}")
            sys.exit(1)
    
    def find_barcode_in_sequence(self, sequence, quality_scores=None):
        """
        Find best barcode match in a sequence using fuzzy matching
        
        Args:
            sequence: DNA sequence string
            quality_scores: List of quality scores (optional)
            
        Returns:
            tuple: (sample_id, orientation, position, mismatches, barcode_quality)
        """
        best_match = None
        best_score = float('inf')
        
        for sample_id, patterns in self.barcode_patterns.items():
            # Try forward barcode
            for orientation, barcode in [('forward', patterns['forward']), 
                                       ('reverse', patterns['reverse'])]:
                
                # Use edlib for approximate matching
                result = edlib.align(barcode, sequence, mode="HW", 
                                   task="locations", k=self.max_mismatches)
                
                if result['editDistance'] != -1 and result['editDistance'] <= self.max_mismatches:
                    edit_distance = result['editDistance']
                    position = result['locations'][0][0] if result['locations'] else 0
                    
                    # Calculate barcode region quality if available
                    barcode_quality = None
                    if quality_scores and position + len(barcode) <= len(quality_scores):
                        barcode_quals = quality_scores[position:position + len(barcode)]
                        barcode_quality = np.mean([ord(q) - 33 for q in barcode_quals])
                    
                    # Prefer matches with fewer mismatches, then by position (closer to start)
                    score = edit_distance * 1000 + position
                    
                    if score < best_score:
                        best_score = score
                        best_match = (sample_id, orientation, position, edit_distance, barcode_quality)
        
        return best_match
    
    def calculate_read_quality(self, quality_string):
        """Calculate average quality score for a read"""
        if not quality_string:
            return 0
        return np.mean([ord(q) - 33 for q in quality_string])
    
    def demultiplex_reads(self):
        """Demultiplex all reads and write to sample-specific files"""
        logger.info("=== Step 2: Demultiplexing Reads ===")
        
        # Open input FASTQ
        fq = pyfastx.Fastq(self.fastq_file)
        total_reads = len(fq)
        
        # Open output files for each sample
        output_files = {}
        sample_writers = {}
        
        # Create output files
        for sample_id in self.barcode_patterns.keys():
            safe_name = sample_id.replace('/', '_').replace(' ', '_')
            output_file = os.path.join(self.output_dirs['demultiplexed_samples'], f"{safe_name}.fastq")
            output_files[sample_id] = output_file
            sample_writers[sample_id] = open(output_file, 'w')
            logger.info(f"  Created: {output_file}")
        
        # Create unassigned reads file
        unassigned_file = os.path.join(self.output_dirs['demultiplexed_samples'], "unassigned.fastq")
        unassigned_writer = open(unassigned_file, 'w')
        logger.info(f"  Created: {unassigned_file}")
        
        # Process reads with progress bar
        logger.info(f"Processing {total_reads:,} reads...")
        
        processed_reads = 0
        low_quality_reads = 0
        
        with tqdm(total=total_reads, desc="Demultiplexing") as pbar:
            for read in fq:
                processed_reads += 1
                pbar.update(1)
                
                read_name = read.name
                sequence = read.seq
                quality = read.qual
                
                # Calculate read quality
                avg_quality = self.calculate_read_quality(quality)
                
                # Skip low quality reads
                if avg_quality < self.min_quality:
                    low_quality_reads += 1
                    self.assignment_stats['low_quality'] += 1
                    continue
                
                # Find barcode match
                match_result = self.find_barcode_in_sequence(sequence, quality)
                
                if match_result:
                    sample_id, orientation, position, mismatches, barcode_quality = match_result
                    
                    # Write to sample file
                    sample_writers[sample_id].write(f"@{read_name}\n{sequence}\n+\n{quality}\n")
                    
                    # Update statistics
                    self.assignment_stats[sample_id] += 1
                    self.quality_stats[sample_id].append(avg_quality)
                    self.mismatch_stats[f"{sample_id}_{mismatches}mm"] += 1
                    
                    # Log detailed match info for first few reads
                    if processed_reads <= 10:
                        logger.info(f"  Read {processed_reads}: {read_name} -> {sample_id} "
                                  f"({orientation}, pos={position}, mm={mismatches})")
                
                else:
                    # Write to unassigned file
                    unassigned_writer.write(f"@{read_name}\n{sequence}\n+\n{quality}\n")
                    self.assignment_stats['unassigned'] += 1
        
        # Close all output files
        for writer in sample_writers.values():
            writer.close()
        unassigned_writer.close()
        
        logger.info(f"✓ Processed {processed_reads:,} reads")
        logger.info(f"  Low quality reads skipped: {low_quality_reads:,}")
        
        # Log assignment summary
        logger.info("Demultiplexing Summary:")
        for sample_id, count in self.assignment_stats.items():
            percentage = (count / processed_reads) * 100
            logger.info(f"  {sample_id}: {count:,} reads ({percentage:.1f}%)")
    
    def generate_statistics(self):
        """Generate comprehensive demultiplexing statistics"""
        logger.info("=== Step 3: Generating Statistics ===")
        
        # Calculate overall statistics
        # Split metrics into post-QC (assigned + unassigned) and pre-QC (includes low_quality)
        total_assigned = sum(count for sample, count in self.assignment_stats.items()
                             if sample not in ['unassigned', 'low_quality'])
        unassigned = self.assignment_stats['unassigned']
        low_quality = self.assignment_stats['low_quality']
        post_qc_total = total_assigned + unassigned
        pre_qc_total = post_qc_total + low_quality
        
        # Headline metric: post-QC assignment efficiency
        assignment_efficiency = (total_assigned / post_qc_total * 100.0) if post_qc_total > 0 else 0.0
        # Secondary metric: overall yield from all reads
        overall_yield = (total_assigned / pre_qc_total * 100.0) if pre_qc_total > 0 else 0.0
        
        logger.info(f"Assignment efficiency (post-QC): {assignment_efficiency:.1f}%")
        logger.info(f"Overall yield (including low-quality): {overall_yield:.1f}%")
        logger.info(f"Total assigned reads: {total_assigned:,}")
        logger.info(f"Unassigned reads: {unassigned:,}")
        logger.info(f"Low-quality reads: {low_quality:,}")
        
        # Prepare detailed statistics
        stats_data = []
        for sample_id, patterns in self.barcode_patterns.items():
            count = self.assignment_stats.get(sample_id, 0)
            avg_quality = np.mean(self.quality_stats[sample_id]) if self.quality_stats[sample_id] else 0
            
            # Percentages use post-QC denominator for consistency with headline metric
            stats_data.append({
                'Sample_ID': sample_id,
                'Sample_No': patterns['sample_no'],
                'Forward_Barcode': patterns['forward'],
                'Reverse_Barcode': patterns['reverse'],
                'Read_Count': count,
                'Percentage': (count / post_qc_total) * 100 if post_qc_total > 0 else 0.0,
                'Avg_Quality': avg_quality,
                'Output_File': f"{sample_id.replace('/', '_').replace(' ', '_')}.fastq"
            })
        
        # Add unassigned reads
        stats_data.append({
            'Sample_ID': 'Unassigned',
            'Sample_No': 'N/A',
            'Forward_Barcode': 'N/A',
            'Reverse_Barcode': 'N/A',
            'Read_Count': self.assignment_stats['unassigned'],
            'Percentage': (unassigned / post_qc_total) * 100 if post_qc_total > 0 else 0.0,
            'Avg_Quality': 0,
            'Output_File': 'unassigned.fastq'
        })
        
        # Save statistics
        stats_df = pd.DataFrame(stats_data)
        stats_file = os.path.join(self.output_dirs['qc_reports'], 'phase2_demultiplexing_stats.csv')
        stats_df.to_csv(stats_file, index=False)
        logger.info(f"✓ Statistics saved to: {stats_file}")
        
        # Return both the post-QC assignment efficiency and overall yield for reporting
        stats_df.attrs['overall_yield'] = overall_yield
        return stats_df, assignment_efficiency
    
    def create_qc_plots(self, stats_df, assignment_efficiency):
        """Create quality control plots for demultiplexing"""
        logger.info("=== Step 4: Creating QC Plots ===")
        
        # Set style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Phase 2: Sample Demultiplexing Quality Control', fontsize=16, fontweight='bold')
        
        # Plot 1: Read assignment counts
        sample_data = stats_df[stats_df['Sample_ID'] != 'Unassigned']
        
        bars = axes[0, 0].bar(range(len(sample_data)), sample_data['Read_Count'])
        axes[0, 0].set_xlabel('Sample')
        axes[0, 0].set_ylabel('Read Count')
        axes[0, 0].set_title('Reads per Sample')
        axes[0, 0].set_xticks(range(len(sample_data)))
        axes[0, 0].set_xticklabels(sample_data['Sample_ID'], rotation=45, ha='right')
        
        # Add value labels on bars
        for i, bar in enumerate(bars):
            height = bar.get_height()
            axes[0, 0].text(bar.get_x() + bar.get_width()/2., height + 50,
                           f'{int(height):,}', ha='center', va='bottom')
        
        # Plot 2: Assignment percentages (pie chart)
        pie_data = stats_df[['Sample_ID', 'Read_Count']].copy()
        pie_data = pie_data[pie_data['Read_Count'] > 0]
        
        colors = plt.cm.Set3(np.linspace(0, 1, len(pie_data)))
        wedges, texts, autotexts = axes[0, 1].pie(pie_data['Read_Count'], 
                                                 labels=pie_data['Sample_ID'],
                                                 autopct='%1.1f%%',
                                                 colors=colors)
        axes[0, 1].set_title('Read Distribution')
        
        # Plot 3: Quality scores by sample
        quality_data = []
        sample_labels = []
        for sample_id in sample_data['Sample_ID']:
            if sample_id in self.quality_stats and self.quality_stats[sample_id]:
                quality_data.append(self.quality_stats[sample_id])
                sample_labels.append(sample_id)
        
        if quality_data:
            bp = axes[0, 2].boxplot(quality_data, tick_labels=sample_labels)
            axes[0, 2].set_xlabel('Sample')
            axes[0, 2].set_ylabel('Average Quality Score')
            axes[0, 2].set_title('Quality Distribution by Sample')
            axes[0, 2].tick_params(axis='x', rotation=45)
        
        # Plot 4: Mismatch distribution
        mismatch_counts = defaultdict(int)
        for key, count in self.mismatch_stats.items():
            if '_' in key and 'mm' in key:
                mismatches = key.split('_')[-1]
                mismatch_counts[mismatches] += count
        
        if mismatch_counts:
            mm_labels = sorted(mismatch_counts.keys())
            mm_values = [mismatch_counts[label] for label in mm_labels]
            
            axes[1, 0].bar(mm_labels, mm_values)
            axes[1, 0].set_xlabel('Mismatches')
            axes[1, 0].set_ylabel('Read Count')
            axes[1, 0].set_title('Barcode Mismatch Distribution')
        
        # Plot 5: Assignment efficiency metrics (post-QC)
        metrics = {
            'Assigned (post-QC)': assignment_efficiency,
            'Unassigned (post-QC)': 100.0 - assignment_efficiency
        }
        
        colors = ['lightgreen', 'lightcoral']
        bars = axes[1, 1].bar(metrics.keys(), metrics.values(), color=colors)
        axes[1, 1].set_ylabel('Percentage (%)')
        axes[1, 1].set_title('Assignment Efficiency')
        axes[1, 1].set_ylim(0, 100)
        
        # Add percentage labels
        for bar, value in zip(bars, metrics.values()):
            axes[1, 1].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1,
                           f'{value:.1f}%', ha='center', va='bottom')
        
        # Plot 6: Summary statistics
        axes[1, 2].axis('off')
        
        # Compute pre-/post-QC totals locally for clear reporting.
        # Post-QC total equals the sum of assigned + unassigned in stats_df.
        post_qc_total = int(stats_df['Read_Count'].sum())
        assigned_reads = int(sample_data['Read_Count'].sum())
        low_quality = int(self.assignment_stats.get('low_quality', 0))
        pre_qc_total = int(post_qc_total + low_quality)
        
        summary_text = f"""
Demultiplexing Summary:

Total Reads (pre-QC): {pre_qc_total:,}
Post-QC Reads: {post_qc_total:,}
Successfully Assigned: {assigned_reads:,}
Assignment Efficiency (post-QC): {assignment_efficiency:.1f}%
Overall Yield (incl. low-quality): {stats_df.attrs.get('overall_yield', 0.0):.1f}%
Unassigned Reads (post-QC): {int(stats_df[stats_df['Sample_ID'] == 'Unassigned']['Read_Count'].iloc[0]):,}
Low-quality Reads: {low_quality:,}

Sample Distribution (post-QC %):
{chr(10).join([f"  {row['Sample_ID']}: {row['Read_Count']:,} ({row['Percentage']:.1f}%)" 
               for _, row in sample_data.iterrows()])}

Max Mismatches: {self.max_mismatches}
Min Quality: {self.min_quality}

Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}
        """
        
        axes[1, 2].text(0.05, 0.95, summary_text, transform=axes[1, 2].transAxes,
                        fontsize=10, verticalalignment='top', fontfamily='monospace',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.5))
        
        plt.tight_layout()
        
        # Save plot
        plot_file = os.path.join(self.output_dirs['qc_reports'], 'phase2_demultiplexing_qc.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        logger.info(f"✓ QC plots saved to: {plot_file}")
        
        plt.close()
    
    def save_summary_report(self, stats_df, assignment_efficiency):
        """Save comprehensive demultiplexing report"""
        logger.info("=== Step 5: Generating Summary Report ===")
        
        # Prepare detailed summary
        summary_data = {
            'pipeline_info': {
                'phase': 'Phase 2 - Sample Demultiplexing',
                'date': datetime.now().isoformat(),
                'input_file': self.fastq_file,
                'parameters': {
                    'max_mismatches': self.max_mismatches,
                    'min_quality': self.min_quality
                }
            },
            'demultiplexing_results': {
                'assignment_efficiency': assignment_efficiency,
                'total_processed': int(stats_df['Read_Count'].sum()),
                'total_assigned': int(stats_df[stats_df['Sample_ID'] != 'Unassigned']['Read_Count'].sum()),
                'unassigned': int(stats_df[stats_df['Sample_ID'] == 'Unassigned']['Read_Count'].iloc[0])
            },
            'sample_statistics': stats_df.to_dict('records'),
            'mismatch_statistics': dict(self.mismatch_stats)
        }
        
        # Save JSON summary
        json_file = os.path.join(self.output_dirs['qc_reports'], 'phase2_summary.json')
        with open(json_file, 'w') as f:
            json.dump(summary_data, f, indent=2, default=str)
        logger.info(f"✓ JSON summary saved to: {json_file}")
        
        # Create human-readable report
        report_file = os.path.join(self.output_dirs['qc_reports'], 'phase2_report.txt')
        with open(report_file, 'w') as f:
            f.write("CHEMOKINE RECEPTOR DEMULTIPLEXING PIPELINE\n")
            f.write("Phase 2: Sample Demultiplexing Results\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Input File: {self.fastq_file}\n")
            f.write(f"Max Mismatches: {self.max_mismatches}\n")
            f.write(f"Min Quality Score: {self.min_quality}\n\n")
            
            f.write("DEMULTIPLEXING RESULTS\n")
            f.write("-" * 30 + "\n")
            f.write(f"Assignment Efficiency: {assignment_efficiency:.1f}%\n")
            f.write(f"Total Reads Processed: {stats_df['Read_Count'].sum():,}\n")
            f.write(f"Successfully Assigned: {stats_df[stats_df['Sample_ID'] != 'Unassigned']['Read_Count'].sum():,}\n")
            f.write(f"Unassigned Reads: {stats_df[stats_df['Sample_ID'] == 'Unassigned']['Read_Count'].iloc[0]:,}\n\n")
            
            f.write("SAMPLE BREAKDOWN\n")
            f.write("-" * 30 + "\n")
            sample_data = stats_df[stats_df['Sample_ID'] != 'Unassigned']
            for _, row in sample_data.iterrows():
                f.write(f"{row['Sample_ID']}: {row['Read_Count']:,} reads ({row['Percentage']:.1f}%)\n")
                f.write(f"  Barcode: {row['Forward_Barcode']} / {row['Reverse_Barcode']}\n")
                f.write(f"  Avg Quality: {row['Avg_Quality']:.1f}\n")
                f.write(f"  Output: {row['Output_File']}\n\n")
            
            f.write("OUTPUT FILES\n")
            f.write("-" * 30 + "\n")
            f.write("Sample-specific FASTQ files:\n")
            for _, row in sample_data.iterrows():
                f.write(f"  demultiplexed_samples/{row['Output_File']}\n")
            f.write("  demultiplexed_samples/unassigned.fastq\n\n")
            
            f.write("QUALITY CONTROL\n")
            f.write("-" * 30 + "\n")
            f.write("- phase2_demultiplexing_qc.png\n")
            f.write("- phase2_demultiplexing_stats.csv\n")
            f.write("- logs/phase2_demultiplexing.log\n\n")
            
            f.write("NEXT STEPS\n")
            f.write("-" * 30 + "\n")
            f.write("1. Review demultiplexing efficiency and sample balance\n")
            f.write("2. Proceed to Phase 3: Chemokine Barcode Quantification\n")
            f.write("3. Use sample-specific FASTQ files for target analysis\n")
        
        logger.info(f"✓ Human-readable report saved to: {report_file}")
    
    def run_pipeline(self):
        """Execute the complete Phase 2 pipeline"""
        logger.info("Starting Phase 2: Sample Demultiplexing")
        
        # Create output directories
        for dir_name in self.output_dirs.values():
            os.makedirs(dir_name, exist_ok=True)
        
        try:
            # Execute pipeline steps
            self.load_barcodes()
            self.demultiplex_reads()
            stats_df, assignment_efficiency = self.generate_statistics()
            self.create_qc_plots(stats_df, assignment_efficiency)
            self.save_summary_report(stats_df, assignment_efficiency)
            
            logger.info("=" * 60)
            logger.info("Phase 2 completed successfully!")
            logger.info(f"Assignment efficiency: {assignment_efficiency:.1f}%")
            logger.info("=" * 60)
            logger.info("Ready to proceed to Phase 3: Chemokine Barcode Quantification")
            
            return True, assignment_efficiency
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False, 0

def main():
    """Main execution function"""
    print("Chemokine Receptor Demultiplexing Pipeline")
    print("Phase 2: Sample Demultiplexing")
    print("=" * 60)
    
    # File paths
    fastq_file = "FCHKLT_1_EL1.fastq"
    barcode_file = "reference_files/sample_barcodes_with_rc.csv"
    
    # Check if files exist
    for file in [fastq_file, barcode_file]:
        if not os.path.exists(file):
            print(f"Error: File not found - {file}")
            print("Please ensure Phase 1 has been completed and files are available.")
            sys.exit(1)
    
    # Initialize and run pipeline
    demultiplexer = SampleDemultiplexer(
        fastq_file=fastq_file,
        barcode_file=barcode_file,
        max_mismatches=1,  # Allow 1 mismatch in barcode
        min_quality=20     # Minimum average quality score
    )
    
    success, efficiency = demultiplexer.run_pipeline()
    
    if success:
        print("\nPhase 2 completed successfully.")
        print(f"Assignment efficiency: {efficiency:.1f}%")
        print("\nGenerated files:")
        print("- demultiplexed_samples/*.fastq (sample-specific files)")
        print("- demultiplexed_samples/unassigned.fastq")
        print("- qc_reports/phase2_demultiplexing_qc.png")
        print("- qc_reports/phase2_demultiplexing_stats.csv")
        print("- qc_reports/phase2_report.txt")
        print("- logs/phase2_demultiplexing.log")
        print("\nReady to proceed to Phase 3: Chemokine Barcode Quantification")
    else:
        print("\nPhase 2 failed. Check logs for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()