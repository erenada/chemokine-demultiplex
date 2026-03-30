#!/usr/bin/env python3

"""
Phase 3: Multi-File Chemokine Barcode Quantification
New Chemokine Receptor Analysis Pipeline

Author: Eren Ada, PhD
Date: 10/22/2025

This script performs:
1. Load chemokine target barcodes with reverse complements  
2. Quantify barcode occurrences in each of 6 FASTQ files
3. Handle both forward and reverse orientations with fuzzy matching
4. Generate comprehensive count matrix (6 samples × 21 chemokines)
5. Apply normalization methods and quality control analysis
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
import logging
from datetime import datetime
import json
from collections import defaultdict
from tqdm import tqdm

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/phase3_multi_quantification.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class MultiBarcodeQuantifier:
    def __init__(self, sample_mapping_file, barcode_file, fastq_dir, 
                 max_mismatches=2, min_quality=20, output_base_dir=None):
        """
        Initialize multi-file barcode quantification pipeline
        
        Args:
            sample_mapping_file: Path to sample mapping CSV
            barcode_file: Path to chemokine barcodes with reverse complements CSV
            fastq_dir: Directory containing FASTQ files
            max_mismatches: Maximum mismatches allowed in barcode matching
            min_quality: Minimum average quality score for reads to process
            output_base_dir: Base directory for all outputs (optional)
        """
        self.sample_mapping_file = sample_mapping_file
        self.barcode_file = barcode_file
        self.fastq_dir = fastq_dir
        self.max_mismatches = max_mismatches
        self.min_quality = min_quality
        
        # Set base output directory
        if output_base_dir is None:
            output_base_dir = '.'
        
        # Output directories
        self.output_dirs = {
            'results': os.path.join(output_base_dir, 'results'),
            'qc_reports': os.path.join(output_base_dir, 'qc_reports'),
            'logs': os.path.join(output_base_dir, 'logs')
        }
        
        # Results storage
        self.sample_mapping = None
        self.target_barcodes = None
        self.barcode_patterns = {}
        self.count_matrix = None
        self.detection_stats = defaultdict(dict)
        self.quality_stats = defaultdict(dict)
        
        logger.info("Multi-File Barcode Quantifier initialized")
        logger.info(f"Sample mapping: {sample_mapping_file}")
        logger.info(f"Barcode file: {barcode_file}")
        logger.info(f"FASTQ directory: {fastq_dir}")
        logger.info(f"Max mismatches: {max_mismatches}")
        logger.info(f"Min quality: {min_quality}")
    
    def load_sample_mapping(self):
        """Load sample mapping file"""
        logger.info("=== Step 1: Loading Sample Mapping ===")
        
        try:
            self.sample_mapping = pd.read_csv(self.sample_mapping_file)
            logger.info(f"Loaded {len(self.sample_mapping)} samples")
            
            for _, row in self.sample_mapping.iterrows():
                logger.info(f"  {row['Sample_annotation']}: {row['fastq_file']}")
            
            return True
        except Exception as e:
            logger.error(f"Error loading sample mapping: {e}")
            return False
    
    def load_target_barcodes(self):
        """Load chemokine target barcodes with reverse complements"""
        logger.info("=== Step 2: Loading Chemokine Target Barcodes ===")
        
        try:
            self.target_barcodes = pd.read_csv(self.barcode_file)
            logger.info(f"Loaded {len(self.target_barcodes)} chemokine target barcodes")
            
            # Prepare barcode patterns for matching
            for _, row in self.target_barcodes.iterrows():
                chemokine = row['Chemokine']
                forward_bc = row['Barcode']
                reverse_bc = row['Barcode_rc']
                
                self.barcode_patterns[chemokine] = {
                    'forward': forward_bc,
                    'reverse': reverse_bc,
                    'length': len(forward_bc)
                }
            
            logger.info(f"Barcode patterns prepared for {len(self.barcode_patterns)} targets")
            return True
            
        except Exception as e:
            logger.error(f"Error loading target barcodes: {e}")
            return False
    
    def find_barcode_in_sequence(self, sequence, quality_scores=None):
        """
        Find best barcode match in a sequence using fuzzy matching
        
        Args:
            sequence: DNA sequence string
            quality_scores: Quality scores string (optional)
            
        Returns:
            list: [(chemokine, orientation, position, mismatches, barcode_quality), ...]
        """
        matches = []
        
        for chemokine, patterns in self.barcode_patterns.items():
            # Try both forward and reverse barcode
            for orientation, barcode in [('forward', patterns['forward']), 
                                       ('reverse', patterns['reverse'])]:
                
                # Use edlib for approximate matching
                result = edlib.align(barcode, sequence, mode="HW", 
                                   task="locations", k=self.max_mismatches)
                
                if result['editDistance'] != -1 and result['editDistance'] <= self.max_mismatches:
                    edit_distance = result['editDistance']
                    
                    # Get all match locations
                    for location in result['locations']:
                        position = location[0]
                        
                        # Calculate barcode region quality if available
                        barcode_quality = None
                        if quality_scores and position + len(barcode) <= len(quality_scores):
                            barcode_quals = quality_scores[position:position + len(barcode)]
                            barcode_quality = np.mean([ord(q) - 33 for q in barcode_quals])
                        
                        matches.append((chemokine, orientation, position, edit_distance, barcode_quality))
        
        return matches
    
    def calculate_read_quality(self, quality_string):
        """Calculate average quality score for a read"""
        if not quality_string:
            return 0
        return np.mean([ord(q) - 33 for q in quality_string])
    
    def quantify_sample(self, sample_row):
        """Quantify barcode occurrences in a single sample"""
        sample_name = sample_row['Sample_annotation']
        fastq_file = os.path.join(self.fastq_dir, sample_row['fastq_file'])
        
        logger.info(f"Processing sample: {sample_name}")
        logger.info(f"  File: {sample_row['fastq_file']}")
        
        # Initialize counts
        barcode_counts = defaultdict(int)
        orientation_counts = defaultdict(lambda: defaultdict(int))
        mismatch_counts = defaultdict(lambda: defaultdict(int))
        quality_scores = []
        
        # Process reads
        fq = pyfastx.Fastq(fastq_file)
        total_reads = len(fq)
        processed_reads = 0
        low_quality_reads = 0
        reads_with_barcodes = 0
        
        logger.info(f"  Processing {total_reads:,} reads...")
        
        for read in tqdm(fq, desc=f"  {sample_name}", leave=False, disable=False):
            processed_reads += 1
            
            sequence = read.seq
            quality = read.qual
            
            # Calculate read quality
            avg_quality = self.calculate_read_quality(quality)
            quality_scores.append(avg_quality)
            
            # Skip low quality reads
            if avg_quality < self.min_quality:
                low_quality_reads += 1
                continue
            
            # Find all candidate barcode matches
            matches = self.find_barcode_in_sequence(sequence, quality)
            
            if matches:
                reads_with_barcodes += 1
                
                # Deduplicate: at most one count per chemokine per read
                # Select best match by: 1) fewest mismatches, 2) earliest position, 3) forward orientation
                best_by_chemokine = {}
                ambiguous_chemokines = set()
                
                for chemokine, orientation, position, mismatches, barcode_qual in matches:
                    orientation_rank = 0 if orientation == 'forward' else 1
                    score = (mismatches, position, orientation_rank)
                    
                    if chemokine not in best_by_chemokine:
                        best_by_chemokine[chemokine] = {
                            'orientation': orientation,
                            'position': position,
                            'mismatches': mismatches,
                            'score': score
                        }
                    else:
                        current = best_by_chemokine[chemokine]
                        if score < current['score']:
                            best_by_chemokine[chemokine] = {
                                'orientation': orientation,
                                'position': position,
                                'mismatches': mismatches,
                                'score': score
                            }
                        elif score == current['score']:
                            # Ambiguous tie
                            ambiguous_chemokines.add(chemokine)
                
                # Count best match per chemokine, skip ambiguous
                for chemokine, info in best_by_chemokine.items():
                    if chemokine in ambiguous_chemokines:
                        continue
                    barcode_counts[chemokine] += 1
                    orientation_counts[chemokine][info['orientation']] += 1
                    mismatch_counts[chemokine][info['mismatches']] += 1
        
        # Calculate statistics
        detection_rate = (reads_with_barcodes / processed_reads) * 100 if processed_reads > 0 else 0
        avg_sample_quality = np.mean(quality_scores) if quality_scores else 0
        
        logger.info(f"  {sample_name}: {reads_with_barcodes:,}/{processed_reads:,} reads with barcodes ({detection_rate:.1f}%)")
        logger.info(f"    Average quality: {avg_sample_quality:.1f}")
        logger.info(f"    Low quality reads: {low_quality_reads:,}")
        
        # Store sample statistics
        self.detection_stats[sample_name] = {
            'total_reads': total_reads,
            'processed_reads': processed_reads,
            'reads_with_barcodes': reads_with_barcodes,
            'detection_rate': detection_rate,
            'low_quality_reads': low_quality_reads
        }
        
        self.quality_stats[sample_name] = {
            'average_quality': avg_sample_quality,
            'quality_distribution': quality_scores
        }
        
        return {
            'sample_name': sample_name,
            'condition': sample_row['Condition'],
            'barcode_counts': dict(barcode_counts),
            'orientation_counts': dict(orientation_counts),
            'mismatch_counts': dict(mismatch_counts)
        }
    
    def quantify_all_samples(self):
        """Quantify barcode occurrences across all samples"""
        logger.info("=== Step 3: Quantifying Barcodes Across All Samples ===")
        
        all_results = []
        
        for _, row in self.sample_mapping.iterrows():
            result = self.quantify_sample(row)
            all_results.append(result)
        
        logger.info("\nQuantification complete for all samples")
        return all_results
    
    def create_count_matrix(self, quantification_results):
        """Create count matrix from quantification results"""
        logger.info("=== Step 4: Creating Count Matrix ===")
        
        # Initialize matrix
        samples = [result['sample_name'] for result in quantification_results]
        chemokines = sorted(list(self.barcode_patterns.keys()))
        
        # Create count matrix
        count_data = []
        for result in quantification_results:
            sample_counts = []
            for chemokine in chemokines:
                count = result['barcode_counts'].get(chemokine, 0)
                sample_counts.append(count)
            count_data.append(sample_counts)
        
        # Create DataFrame
        self.count_matrix = pd.DataFrame(
            count_data, 
            index=samples, 
            columns=chemokines
        )
        
        # Calculate summary statistics
        total_counts = self.count_matrix.sum().sum()
        detected_barcodes = (self.count_matrix > 0).sum().sum()
        detection_rate = (detected_barcodes / (len(samples) * len(chemokines))) * 100
        
        logger.info(f"Count matrix created: {len(samples)} samples × {len(chemokines)} chemokines")
        logger.info(f"  Total barcode detections: {total_counts:,}")
        logger.info(f"  Barcode detection rate: {detection_rate:.1f}%")
        
        # Save raw count matrix
        count_file = os.path.join(self.output_dirs['results'], 'raw_count_matrix_6x21.csv')
        self.count_matrix.to_csv(count_file)
        logger.info(f"Raw count matrix saved to: {count_file}")
        
        return self.count_matrix
    
    def normalize_counts(self):
        """Apply various normalization methods"""
        logger.info("=== Step 5: Applying Normalization Methods ===")
        
        normalized_matrices = {}
        
        # 1. Counts Per Million (CPM)
        total_counts_per_sample = self.count_matrix.sum(axis=1)
        safe_totals = total_counts_per_sample.replace(0, np.nan)
        cpm_matrix = self.count_matrix.div(safe_totals, axis=0) * 1e6
        cpm_matrix = cpm_matrix.fillna(0.0)
        normalized_matrices['CPM'] = cpm_matrix
        
        # 2. Transcripts Per Million (TPM)
        rpk = self.count_matrix
        rpk_sum = rpk.sum(axis=1)
        safe_rpk_sum = rpk_sum.replace(0, np.nan)
        tpm_matrix = rpk.div(safe_rpk_sum, axis=0) * 1e6
        tpm_matrix = tpm_matrix.fillna(0.0)
        normalized_matrices['TPM'] = tpm_matrix
        
        # 3. Log2(CPM + 1)
        log_cpm_matrix = np.log2(cpm_matrix + 1.0)
        normalized_matrices['Log2_CPM'] = log_cpm_matrix
        
        # 4. Relative abundance (proportions)
        relative_matrix = self.count_matrix.div(safe_totals, axis=0)
        relative_matrix = relative_matrix.fillna(0.0)
        normalized_matrices['Relative'] = relative_matrix
        
        # Save normalized matrices
        for method, matrix in normalized_matrices.items():
            norm_file = os.path.join(self.output_dirs['results'], f'{method.lower()}_normalized_counts.csv')
            matrix.to_csv(norm_file)
            logger.info(f"  {method} normalized counts saved")
        
        return normalized_matrices
    
    def create_qc_plots(self, quantification_results, normalized_matrices):
        """Create comprehensive quality control plots"""
        logger.info("=== Step 6: Creating Quality Control Plots ===")
        
        plt.style.use('default')
        sns.set_palette("husl")
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Phase 3: Barcode Quantification Quality Control', fontsize=16, fontweight='bold')
        
        # Plot 1: Total counts per sample
        # Create better labels showing condition and replicate
        sample_labels_short = []
        for s in self.count_matrix.index:
            if 'GFP' in s:
                rep_num = s.split('Rep')[-1]
                sample_labels_short.append(f'GFP-Rep{rep_num}')
            elif 'Neg' in s:
                rep_num = s.split('Rep')[-1]
                sample_labels_short.append(f'Neg-Rep{rep_num}')
            else:
                sample_labels_short.append(s.split('-')[-1])
        
        sample_totals = self.count_matrix.sum(axis=1)
        bars = axes[0, 0].bar(range(len(sample_totals)), sample_totals.values, alpha=0.7, edgecolor='black')
        axes[0, 0].set_xlabel('Sample')
        axes[0, 0].set_ylabel('Total Barcode Counts')
        axes[0, 0].set_title('Total Counts per Sample')
        axes[0, 0].set_xticks(range(len(sample_totals)))
        axes[0, 0].set_xticklabels(sample_labels_short, rotation=45, ha='right')
        
        # Add value labels
        for i, bar in enumerate(bars):
            height = bar.get_height()
            axes[0, 0].text(bar.get_x() + bar.get_width()/2., height,
                           f'{int(height):,}', ha='center', va='bottom', fontsize=8)
        
        # Plot 2: Detection rate per chemokine
        detection_rates = (self.count_matrix > 0).sum(axis=0) / len(self.count_matrix) * 100
        bars = axes[0, 1].bar(range(len(detection_rates)), detection_rates.values, alpha=0.7, edgecolor='black')
        axes[0, 1].set_xlabel('Chemokine Target')
        axes[0, 1].set_ylabel('Detection Rate (%)')
        axes[0, 1].set_title('Detection Rate per Chemokine')
        axes[0, 1].set_xticks(range(len(detection_rates)))
        axes[0, 1].set_xticklabels(detection_rates.index, rotation=90, ha='right', fontsize=7)
        axes[0, 1].set_ylim(0, 100)
        
        # Plot 3: Normalized expression heatmap
        # Create better labels showing condition and replicate
        sample_labels = []
        for s in normalized_matrices['Log2_CPM'].index:
            if 'GFP' in s:
                rep_num = s.split('Rep')[-1]
                sample_labels.append(f'GFP-Rep{rep_num}')
            elif 'Neg' in s:
                rep_num = s.split('Rep')[-1]
                sample_labels.append(f'Neg-Rep{rep_num}')
            else:
                sample_labels.append(s.split('-')[-1])
        
        sns.heatmap(normalized_matrices['Log2_CPM'].T, ax=axes[0, 2], cmap='RdYlBu_r', 
                    center=0, robust=True, cbar_kws={'label': 'Log2(CPM+1)'})
        axes[0, 2].set_title('Normalized Expression Heatmap')
        axes[0, 2].set_xlabel('Sample')
        axes[0, 2].set_ylabel('Chemokine Target')
        axes[0, 2].set_xticklabels(sample_labels, rotation=45, ha='right', fontsize=8)
        axes[0, 2].set_yticklabels(axes[0, 2].get_yticklabels(), fontsize=7)
        
        # Plot 4: Sample correlation matrix
        correlation_matrix = normalized_matrices['Log2_CPM'].T.corr(method='spearman')
        sns.heatmap(correlation_matrix, ax=axes[1, 0], cmap='coolwarm', center=0,
                    vmin=-1, vmax=1, annot=True, fmt='.2f', cbar_kws={'label': 'Spearman r'})
        axes[1, 0].set_title('Sample Correlation Matrix')
        axes[1, 0].set_xticklabels(sample_labels, rotation=45, ha='right', fontsize=8)
        axes[1, 0].set_yticklabels(sample_labels, rotation=0, fontsize=8)
        
        # Plot 5: Detection statistics
        axes[1, 1].axis('off')
        
        total_barcodes = (self.count_matrix > 0).sum().sum()
        total_possible = len(self.count_matrix) * len(self.count_matrix.columns)
        overall_detection = (total_barcodes / total_possible) * 100
        
        stats_text = f"""
Detection Statistics:

Total Samples: {len(self.count_matrix)}
Total Chemokine Targets: {len(self.count_matrix.columns)}

Total Barcode Detections: {self.count_matrix.sum().sum():,}
Overall Detection Rate: {overall_detection:.1f}%

Per-Sample Detection Rates:
{chr(10).join([f"  {s.split('-')[-1]}: {self.detection_stats[s]['detection_rate']:.1f}%" for s in self.count_matrix.index])}

Top 5 Detected Chemokines:
{chr(10).join([f"  {chem}: {count:,}" for chem, count in self.count_matrix.sum().nlargest(5).items()])}
        """
        
        axes[1, 1].text(0.05, 0.95, stats_text, transform=axes[1, 1].transAxes,
                        fontsize=9, verticalalignment='top', fontfamily='monospace',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.5))
        
        # Plot 6: Condition comparison
        axes[1, 2].axis('off')
        
        # Calculate average expression by condition
        gfp_samples = [s for s in self.count_matrix.index if 'GFP' in s]
        neg_samples = [s for s in self.count_matrix.index if 'Neg' in s]
        
        if gfp_samples and neg_samples:
            gfp_avg = self.count_matrix.loc[gfp_samples].mean()
            neg_avg = self.count_matrix.loc[neg_samples].mean()
            
            comparison_text = f"""
Condition Comparison:

GFP Samples (n={len(gfp_samples)}):
  Total Counts: {self.count_matrix.loc[gfp_samples].sum().sum():,}
  Avg per Sample: {self.count_matrix.loc[gfp_samples].sum(axis=1).mean():.0f}

Negative Samples (n={len(neg_samples)}):
  Total Counts: {self.count_matrix.loc[neg_samples].sum().sum():,}
  Avg per Sample: {self.count_matrix.loc[neg_samples].sum(axis=1).mean():.0f}

Top Enriched in GFP:
{chr(10).join([f"  {chem}: {gfp_avg[chem]:.0f}" for chem in gfp_avg.nlargest(5).index])}

Top Enriched in Negative:
{chr(10).join([f"  {chem}: {neg_avg[chem]:.0f}" for chem in neg_avg.nlargest(5).index])}
            """
        else:
            comparison_text = "\nCondition comparison not available"
        
        axes[1, 2].text(0.05, 0.95, comparison_text, transform=axes[1, 2].transAxes,
                        fontsize=9, verticalalignment='top', fontfamily='monospace',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.3))
        
        plt.tight_layout()
        
        # Save plot
        plot_file = os.path.join(self.output_dirs['qc_reports'], 'phase3_quantification_qc.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        logger.info(f"QC plots saved to: {plot_file}")
        
        plt.close()
    
    def save_summary_report(self, quantification_results):
        """Save comprehensive summary report"""
        logger.info("=== Step 7: Generating Summary Report ===")
        
        # JSON summary
        summary = {
            'pipeline_info': {
                'phase': 'Phase 3 - Multi-File Barcode Quantification',
                'date': datetime.now().isoformat(),
                'total_samples': len(self.sample_mapping),
                'total_targets': len(self.target_barcodes)
            },
            'quantification_stats': {
                'total_barcode_detections': int(self.count_matrix.sum().sum()),
                'detection_rate': float((self.count_matrix > 0).sum().sum() / 
                                      (len(self.count_matrix) * len(self.count_matrix.columns)) * 100)
            },
            'per_sample_stats': {k: {
                'total_reads': v['total_reads'],
                'reads_with_barcodes': v['reads_with_barcodes'],
                'detection_rate': v['detection_rate']
            } for k, v in self.detection_stats.items()}
        }
        
        json_file = os.path.join(self.output_dirs['qc_reports'], 'phase3_summary.json')
        with open(json_file, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        
        logger.info(f"Summary report saved to: {json_file}")
        
        # Human-readable report
        report_file = os.path.join(self.output_dirs['qc_reports'], 'phase3_quantification_report.txt')
        with open(report_file, 'w') as f:
            f.write("NEW CHEMOKINE RECEPTOR ANALYSIS PIPELINE\n")
            f.write("Phase 3: Multi-File Barcode Quantification\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Analysis Date: {datetime.now().strftime('%m/%d/%Y')}\n")
            f.write(f"Samples Processed: {len(self.sample_mapping)}\n")
            f.write(f"Chemokine Targets: {len(self.target_barcodes)}\n\n")
            
            f.write("QUANTIFICATION RESULTS\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total Barcode Detections: {self.count_matrix.sum().sum():,}\n")
            f.write(f"Overall Detection Rate: {(self.count_matrix > 0).sum().sum() / (len(self.count_matrix) * len(self.count_matrix.columns)) * 100:.1f}%\n\n")
            
            f.write("PER-SAMPLE STATISTICS\n")
            f.write("-" * 40 + "\n")
            for sample in self.count_matrix.index:
                stats = self.detection_stats[sample]
                f.write(f"\n{sample}:\n")
                f.write(f"  Total Reads: {stats['total_reads']:,}\n")
                f.write(f"  Reads with Barcodes: {stats['reads_with_barcodes']:,}\n")
                f.write(f"  Detection Rate: {stats['detection_rate']:.1f}%\n")
                f.write(f"  Total Barcode Counts: {self.count_matrix.loc[sample].sum():,}\n")
            
            f.write("\nTOP DETECTED CHEMOKINES\n")
            f.write("-" * 40 + "\n")
            for chemokine, count in self.count_matrix.sum().nlargest(10).items():
                f.write(f"  {chemokine}: {count:,}\n")
            
            f.write("\nFILES GENERATED\n")
            f.write("-" * 40 + "\n")
            f.write("- results/raw_count_matrix_6x21.csv\n")
            f.write("- results/cpm_normalized_counts.csv\n")
            f.write("- results/tpm_normalized_counts.csv\n")
            f.write("- results/log2_cpm_normalized_counts.csv\n")
            f.write("- results/relative_normalized_counts.csv\n")
            f.write("- qc_reports/phase3_quantification_qc.png\n")
            f.write("- qc_reports/phase3_summary.json\n")
            f.write("- logs/phase3_multi_quantification.log\n\n")
            
            f.write("NEXT STEPS\n")
            f.write("-" * 40 + "\n")
            f.write("1. Review count matrix and QC plots\n")
            f.write("2. Analyze condition-specific patterns (GFP vs Negative)\n")
            f.write("3. Perform statistical differential expression analysis\n")
        
        logger.info(f"Human-readable report saved to: {report_file}")
    
    def run_pipeline(self):
        """Execute the complete Phase 3 pipeline"""
        logger.info("Starting Phase 3: Multi-File Barcode Quantification")
        
        # Create output directories
        for dir_name in self.output_dirs.values():
            os.makedirs(dir_name, exist_ok=True)
        
        try:
            # Execute pipeline steps
            if not self.load_sample_mapping():
                return False
            
            if not self.load_target_barcodes():
                return False
            
            quantification_results = self.quantify_all_samples()
            self.create_count_matrix(quantification_results)
            normalized_matrices = self.normalize_counts()
            self.create_qc_plots(quantification_results, normalized_matrices)
            self.save_summary_report(quantification_results)
            
            logger.info("=" * 70)
            logger.info("Phase 3 completed successfully!")
            logger.info("=" * 70)
            
            return True
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False

def main():
    """Main execution function"""
    print("New Chemokine Receptor Analysis Pipeline")
    print("Phase 3: Multi-File Barcode Quantification")
    print("=" * 70)
    
    # File paths
    sample_mapping = "New_analysis/reference_files/sample_mapping.csv"
    barcode_file = "New_analysis/data_nov14/reference_files/chemokine_barcodes_with_rc.csv"
    fastq_dir = "New_analysis/data_nov14/SYFFKS_fastq"
    output_base_dir = "New_analysis/data_nov14"
    
    # Check if files exist
    for file in [sample_mapping, barcode_file]:
        if not os.path.exists(file):
            print(f"Error: File not found - {file}")
            sys.exit(1)
    
    if not os.path.exists(fastq_dir):
        print(f"Error: Directory not found - {fastq_dir}")
        sys.exit(1)
    
    # Initialize and run pipeline
    pipeline = MultiBarcodeQuantifier(
        sample_mapping_file=sample_mapping,
        barcode_file=barcode_file,
        fastq_dir=fastq_dir,
        max_mismatches=2,
        min_quality=20,
        output_base_dir=output_base_dir
    )
    
    success = pipeline.run_pipeline()
    
    if success:
        print("\nPhase 3 completed successfully.")
        print("\nGenerated files:")
        print(f"- {output_base_dir}/results/raw_count_matrix_6x21.csv")
        print(f"- {output_base_dir}/results/cpm_normalized_counts.csv")
        print(f"- {output_base_dir}/results/tpm_normalized_counts.csv")
        print(f"- {output_base_dir}/results/log2_cpm_normalized_counts.csv")
        print(f"- {output_base_dir}/results/relative_normalized_counts.csv")
        print(f"- {output_base_dir}/qc_reports/phase3_quantification_qc.png")
        print(f"- {output_base_dir}/qc_reports/phase3_quantification_report.txt")
        print(f"- {output_base_dir}/qc_reports/phase3_summary.json")
        print("\nAnalysis complete. Count matrix ready for statistical analysis.")
    else:
        print("\nPhase 3 failed. Check logs for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()

