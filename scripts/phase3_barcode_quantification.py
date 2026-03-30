#!/usr/bin/env python3

"""
Phase 3: Chemokine Barcode Quantification
Chemokine Receptor Demultiplexing Pipeline

Author: Eren Ada, PhD
Date: 08/07/2025

This script performs:
1. Load chemokine target barcodes with reverse complements
2. Quantify barcode occurrences in each sample-specific FASTQ file
3. Handle both forward and reverse orientations with fuzzy matching
4. Generate comprehensive count matrix (samples × chemokines)
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
import regex
import logging
from datetime import datetime
import json
from collections import defaultdict, Counter
from tqdm import tqdm
import glob

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/phase3_quantification.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class BarcodeQuantifier:
    def __init__(self, barcode_file, sample_dir, max_mismatches=2, min_quality=20):
        """
        Initialize barcode quantification pipeline
        
        Args:
            barcode_file: Path to chemokine barcodes with reverse complements CSV
            sample_dir: Directory containing sample-specific FASTQ files
            max_mismatches: Maximum mismatches allowed in barcode matching
            min_quality: Minimum average quality score for reads to process
        """
        self.barcode_file = barcode_file
        self.sample_dir = sample_dir
        self.max_mismatches = max_mismatches
        self.min_quality = min_quality
        
        # Output directories
        self.output_dirs = {
            'results': 'results',
            'qc_reports': 'qc_reports',
            'logs': 'logs'
        }
        
        # Results storage
        self.target_barcodes = None
        self.barcode_patterns = {}
        self.sample_files = []
        self.count_matrix = None
        self.detection_stats = defaultdict(dict)
        self.quality_stats = defaultdict(dict)
        
        logger.info("Barcode Quantifier initialized")
        logger.info(f"Barcode file: {barcode_file}")
        logger.info(f"Sample directory: {sample_dir}")
        logger.info(f"Max mismatches: {max_mismatches}")
        logger.info(f"Min quality: {min_quality}")
    
    def load_target_barcodes(self):
        """Load chemokine target barcodes with reverse complements"""
        logger.info("=== Step 1: Loading Chemokine Target Barcodes ===")
        
        try:
            self.target_barcodes = pd.read_csv(self.barcode_file)
            logger.info(f"✓ Loaded {len(self.target_barcodes)} chemokine target barcodes")
            
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
                
                logger.info(f"  {chemokine}: {forward_bc} / {reverse_bc}")
            
        except Exception as e:
            logger.error(f"Error loading target barcodes: {e}")
            sys.exit(1)
    
    def find_sample_files(self):
        """Find all sample-specific FASTQ files"""
        logger.info("=== Step 2: Finding Sample FASTQ Files ===")
        
        # Find all FASTQ files except unassigned
        pattern = os.path.join(self.sample_dir, "*.fastq")
        all_files = glob.glob(pattern)
        
        # Filter out unassigned file
        self.sample_files = [f for f in all_files if 'unassigned' not in os.path.basename(f)]
        self.sample_files.sort()
        
        if not self.sample_files:
            logger.error(f"No sample FASTQ files found in {self.sample_dir}")
            sys.exit(1)
        
        logger.info(f"Found {len(self.sample_files)} sample files:")
        for sample_file in self.sample_files:
            sample_name = os.path.basename(sample_file).replace('.fastq', '')
            fq = pyfastx.Fastq(sample_file)
            read_count = len(fq)
            logger.info(f"  {sample_name}: {read_count:,} reads")
    
    def find_barcode_in_sequence(self, sequence, quality_scores=None):
        """
        Find best barcode match in a sequence using fuzzy matching
        
        Args:
            sequence: DNA sequence string
            quality_scores: List of quality scores (optional)
            
        Returns:
            list: [(chemokine, orientation, position, mismatches, barcode_quality), ...]
        """
        matches = []
        
        for chemokine, patterns in self.barcode_patterns.items():
            # Try forward barcode
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
    
    def quantify_sample(self, sample_file):
        """Quantify barcode occurrences in a single sample"""
        sample_name = os.path.basename(sample_file).replace('.fastq', '')
        logger.info(f"Processing sample: {sample_name}")
        
        # Initialize counts
        barcode_counts = defaultdict(int)
        orientation_counts = defaultdict(lambda: defaultdict(int))
        mismatch_counts = defaultdict(lambda: defaultdict(int))
        quality_scores = []
        detection_positions = defaultdict(list)
        
        # Process reads
        fq = pyfastx.Fastq(sample_file)
        total_reads = len(fq)
        processed_reads = 0
        low_quality_reads = 0
        reads_with_barcodes = 0
        
        logger.info(f"  Processing {total_reads:,} reads...")
        
        for read in tqdm(fq, desc=f"  {sample_name}", leave=False):
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
            
            # Find all candidate barcode matches (may include multiple per chemokine)
            matches = self.find_barcode_in_sequence(sequence, quality)
            
            if matches:
                reads_with_barcodes += 1
                
                # Deduplicate to at most one count per chemokine per read.
                # Select the best match for each chemokine using a deterministic score:
                #   1) fewer mismatches (primary), 2) earlier position (secondary),
                #   3) prefer forward orientation over reverse (tertiary) to break ties.
                # If a tie persists (identical score), mark as ambiguous and skip counting for that chemokine.
                best_by_chemokine = {}
                ambiguous_chemokines = set()
                for chemokine, orientation, position, mismatches, barcode_qual in matches:
                    # Orientation rank: 0 for forward, 1 for reverse (forward preferred)
                    orientation_rank = 0 if orientation == 'forward' else 1
                    score = (mismatches, position, orientation_rank)
                    if chemokine not in best_by_chemokine:
                        best_by_chemokine[chemokine] = {
                            'orientation': orientation,
                            'position': position,
                            'mismatches': mismatches,
                            'barcode_qual': barcode_qual,
                            'score': score
                        }
                    else:
                        current = best_by_chemokine[chemokine]
                        if score < current['score']:
                            best_by_chemokine[chemokine] = {
                                'orientation': orientation,
                                'position': position,
                                'mismatches': mismatches,
                                'barcode_qual': barcode_qual,
                                'score': score
                            }
                        elif score == current['score']:
                            # Ambiguous: two matches are equally good for this chemokine in this read
                            ambiguous_chemokines.add(chemokine)
                
                # Apply counts: one per chemokine per read, skipping ambiguous ties
                for chemokine, info in best_by_chemokine.items():
                    if chemokine in ambiguous_chemokines:
                        continue
                    barcode_counts[chemokine] += 1
                    orientation_counts[chemokine][info['orientation']] += 1
                    mismatch_counts[chemokine][info['mismatches']] += 1
                    detection_positions[chemokine].append(info['position'])
        
        # Calculate statistics
        detection_rate = (reads_with_barcodes / processed_reads) * 100 if processed_reads > 0 else 0
        avg_sample_quality = np.mean(quality_scores) if quality_scores else 0
        
        logger.info(f"  ✓ {sample_name}: {reads_with_barcodes:,}/{processed_reads:,} reads with barcodes ({detection_rate:.1f}%)")
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
            'barcode_counts': dict(barcode_counts),
            'orientation_counts': dict(orientation_counts),
            'mismatch_counts': dict(mismatch_counts),
            'detection_positions': dict(detection_positions)
        }
    
    def quantify_all_samples(self):
        """Quantify barcode occurrences across all samples"""
        logger.info("=== Step 3: Quantifying Barcodes Across All Samples ===")
        
        all_results = []
        
        for sample_file in self.sample_files:
            result = self.quantify_sample(sample_file)
            all_results.append(result)
        
        return all_results
    
    def create_count_matrix(self, quantification_results):
        """Create count matrix from quantification results"""
        logger.info("=== Step 4: Creating Count Matrix ===")
        
        # Initialize matrix
        samples = [result['sample_name'] for result in quantification_results]
        chemokines = list(self.barcode_patterns.keys())
        
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
        
        logger.info(f"✓ Count matrix created: {len(samples)} samples × {len(chemokines)} chemokines")
        logger.info(f"  Total barcode detections: {total_counts:,}")
        logger.info(f"  Barcode detection rate: {detection_rate:.1f}%")
        
        # Save raw count matrix
        count_file = os.path.join(self.output_dirs['results'], 'raw_count_matrix.csv')
        self.count_matrix.to_csv(count_file)
        logger.info(f"✓ Raw count matrix saved to: {count_file}")
        
        return self.count_matrix
    
    def normalize_counts(self):
        """Apply various normalization methods"""
        logger.info("=== Step 5: Applying Normalization Methods ===")
        
        normalized_matrices = {}
        
        # Guard divisions by zero in per-sample totals by replacing zeros with NaN
        # and filling NaN with 0 after division. This prevents inf/NaN artifacts.
        
        # 1. Counts Per Million (CPM)
        total_counts_per_sample = self.count_matrix.sum(axis=1)
        safe_totals = total_counts_per_sample.replace(0, np.nan)
        cpm_matrix = self.count_matrix.div(safe_totals, axis=0) * 1e6
        cpm_matrix = cpm_matrix.fillna(0.0)
        normalized_matrices['CPM'] = cpm_matrix
        
        # 2. Transcripts Per Million (TPM) - treating each barcode as equal length
        # Here RPK equals raw counts; we normalize by per-sample sum, again guarding zeros
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
            logger.info(f"✓ {method} normalized counts saved to: {norm_file}")
        
        return normalized_matrices
    
    def create_qc_plots(self, quantification_results, normalized_matrices):
        """Create comprehensive quality control plots"""
        logger.info("=== Step 6: Creating Quality Control Plots ===")
        
        # Set style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # Create figure with subplots
        fig, axes = plt.subplots(3, 3, figsize=(20, 18))
        fig.suptitle('Phase 3: Chemokine Barcode Quantification Quality Control', fontsize=16, fontweight='bold')
        
        # Plot 1: Total counts per sample
        sample_totals = self.count_matrix.sum(axis=1)
        bars = axes[0, 0].bar(range(len(sample_totals)), sample_totals.values)
        axes[0, 0].set_xlabel('Sample')
        axes[0, 0].set_ylabel('Total Barcode Counts')
        axes[0, 0].set_title('Total Counts per Sample')
        axes[0, 0].set_xticks(range(len(sample_totals)))
        axes[0, 0].set_xticklabels(sample_totals.index, rotation=45, ha='right')
        
        # Add value labels
        for i, bar in enumerate(bars):
            height = bar.get_height()
            axes[0, 0].text(bar.get_x() + bar.get_width()/2., height + 10,
                           f'{int(height):,}', ha='center', va='bottom')
        
        # Plot 2: Detection rate per chemokine (moved from plot 3)
        detection_rates = (self.count_matrix > 0).sum(axis=0) / len(self.count_matrix) * 100
        bars = axes[0, 1].bar(range(len(detection_rates)), detection_rates.values)
        axes[0, 1].set_xlabel('Chemokine Target')
        axes[0, 1].set_ylabel('Detection Rate (%)')
        axes[0, 1].set_title('Detection Rate per Chemokine')
        axes[0, 1].set_xticks(range(len(detection_rates)))
        axes[0, 1].set_xticklabels(detection_rates.index, rotation=45, ha='right')
        axes[0, 1].set_ylim(0, 100)
        
        # Plot 3: Hierarchically clustered normalized expression heatmap with condition grouping
        # Extract condition information from sample names (e.g., Sp1-GFP -> GFP)
        sample_names = normalized_matrices['Log2_CPM'].index.tolist()
        conditions = [name.split('-')[1] if '-' in name else 'Unknown' for name in sample_names]
        
        # Separate samples by condition
        gfp_samples = [s for s in sample_names if 'GFP' in s]
        rfp_samples = [s for s in sample_names if 'RFP' in s]
        
        # Hierarchical clustering within each condition using correlation distance
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import pdist
        
        ordered_samples = []
        
        # Cluster GFP samples
        if len(gfp_samples) > 1:
            gfp_data = normalized_matrices['Log2_CPM'].loc[gfp_samples]
            gfp_distances = pdist(gfp_data, metric='correlation')
            gfp_linkage = linkage(gfp_distances, method='average')
            gfp_order = leaves_list(gfp_linkage)
            ordered_gfp = [gfp_samples[i] for i in gfp_order]
        else:
            ordered_gfp = gfp_samples
        
        # Cluster RFP samples  
        if len(rfp_samples) > 1:
            rfp_data = normalized_matrices['Log2_CPM'].loc[rfp_samples]
            rfp_distances = pdist(rfp_data, metric='correlation')
            rfp_linkage = linkage(rfp_distances, method='average')
            rfp_order = leaves_list(rfp_linkage)
            ordered_rfp = [rfp_samples[i] for i in rfp_order]
        else:
            ordered_rfp = rfp_samples
        
        # Combine ordered samples (GFP first, then RFP)
        ordered_samples = ordered_gfp + ordered_rfp
        
        # ALSO CLUSTER CHEMOKINES by expression similarity
        # This will group chemokines with similar expression patterns together
        chemokine_data = normalized_matrices['Log2_CPM'].T  # chemokines x samples
        
        if len(chemokine_data) > 1:
            # Calculate distances between chemokines based on expression patterns
            chemokine_distances = pdist(chemokine_data.values, metric='correlation')
            chemokine_linkage = linkage(chemokine_distances, method='average')
            chemokine_order = leaves_list(chemokine_linkage)
            
            # Get ordered chemokine names
            ordered_chemokines = [chemokine_data.index[i] for i in chemokine_order]
        else:
            ordered_chemokines = chemokine_data.index.tolist()
        
        # Create color map for the ordered samples
        ordered_colors = ['lightgreen' if 'GFP' in s else 'lightcoral' if 'RFP' in s else 'gray' 
                         for s in ordered_samples]
        
        # Reorder the matrix: ordered samples (columns) AND ordered chemokines (rows)
        ordered_matrix = normalized_matrices['Log2_CPM'].loc[ordered_samples, ordered_chemokines]
        
        # Create the doubly clustered heatmap (both samples and chemokines clustered)
        sns.heatmap(ordered_matrix.T, ax=axes[0, 2], cmap='RdYlBu_r', center=0, robust=True,
                    cbar_kws={'label': 'Log2(CPM+1)'})
        axes[0, 2].set_title('Normalized Expression Heatmap\n(Samples & Chemokines Clustered)')
        axes[0, 2].set_xlabel('Sample (GFP | RFP - Clustered)')
        axes[0, 2].set_ylabel('Chemokine Target (Clustered by Expression)')
        
        # Add condition color bar at the bottom
        for i, color in enumerate(ordered_colors):
            axes[0, 2].add_patch(plt.Rectangle((i, -0.5), 1, 0.3, facecolor=color, clip_on=False))
        
        # Plot 4: Clustered sample correlation matrix
        # Use Spearman correlation which is more robust for expression data
        # and doesn't assume normality
        correlation_matrix = normalized_matrices['Log2_CPM'].T.corr(method='spearman')
        
        # Use the same sample ordering as the expression heatmap for consistency
        ordered_correlation = correlation_matrix.loc[ordered_samples, ordered_samples]
        
        # Create clustered correlation heatmap
        sns.heatmap(ordered_correlation, ax=axes[1, 0], annot=True, cmap='RdYlBu_r',
                    center=0, square=True, cbar_kws={'label': 'Spearman Correlation'},
                    vmin=-1, vmax=1)
        axes[1, 0].set_title('Sample Correlation Matrix\n(Grouped by Condition)')
        axes[1, 0].set_xlabel('Sample (GFP | RFP)')
        axes[1, 0].set_ylabel('Sample (GFP | RFP)')
        
        # Add condition color bars for both axes
        for i, color in enumerate(ordered_colors):
            # Bottom color bar
            axes[1, 0].add_patch(plt.Rectangle((i, -0.5), 1, 0.3, facecolor=color, clip_on=False))
            # Left color bar  
            axes[1, 0].add_patch(plt.Rectangle((-0.5, i), 0.3, 1, facecolor=color, clip_on=False))
        
        # Plot 5: Top expressed chemokines
        top_chemokines = self.count_matrix.sum(axis=0).sort_values(ascending=False).head(10)
        bars = axes[1, 1].bar(range(len(top_chemokines)), top_chemokines.values)
        axes[1, 1].set_xlabel('Chemokine Target')
        axes[1, 1].set_ylabel('Total Counts Across Samples')
        axes[1, 1].set_title('Top 10 Detected Chemokines')
        axes[1, 1].set_xticks(range(len(top_chemokines)))
        axes[1, 1].set_xticklabels(top_chemokines.index, rotation=45, ha='right')
        
        # Plot 6: Quality score distribution by sample
        quality_data = []
        sample_labels = []
        for sample_name in self.quality_stats.keys():
            if 'quality_distribution' in self.quality_stats[sample_name]:
                quality_data.append(self.quality_stats[sample_name]['quality_distribution'])
                sample_labels.append(sample_name)
        
        if quality_data:
            bp = axes[2, 0].boxplot(quality_data, tick_labels=sample_labels)
            axes[2, 0].set_xlabel('Sample')
            axes[2, 0].set_ylabel('Average Quality Score')
            axes[2, 0].set_title('Quality Distribution by Sample')
            axes[2, 0].tick_params(axis='x', rotation=45)
        
        # Plot 7: Detection efficiency per sample
        sample_names = list(self.detection_stats.keys())
        detection_rates = [self.detection_stats[s]['detection_rate'] for s in sample_names]
        
        bars = axes[2, 1].bar(range(len(sample_names)), detection_rates)
        axes[2, 1].set_xlabel('Sample')
        axes[2, 1].set_ylabel('Detection Rate (%)')
        axes[2, 1].set_title('Barcode Detection Rate per Sample')
        axes[2, 1].set_xticks(range(len(sample_names)))
        axes[2, 1].set_xticklabels(sample_names, rotation=45, ha='right')
        axes[2, 1].set_ylim(0, 100)
        
        # Plot 8: Summary statistics
        axes[2, 2].axis('off')
        
        total_detections = self.count_matrix.sum().sum()
        avg_detection_rate = np.mean([self.detection_stats[s]['detection_rate'] for s in self.detection_stats])
        detected_targets = (self.count_matrix.sum() > 0).sum()
        
        summary_text = f"""
Quantification Summary:

Total Samples: {len(self.count_matrix)}
Total Chemokine Targets: {len(self.count_matrix.columns)}
Total Barcode Detections: {total_detections:,}

Detection Statistics:
  Avg Detection Rate: {avg_detection_rate:.1f}%
  Targets Detected: {detected_targets}/{len(self.count_matrix.columns)}
  
Sample Balance:
  Min Counts: {self.count_matrix.sum(axis=1).min():,}
  Max Counts: {self.count_matrix.sum(axis=1).max():,}
  Mean Counts: {self.count_matrix.sum(axis=1).mean():.0f}

Top Chemokines:
{chr(10).join([f"  {target}: {count:,}" for target, count in top_chemokines.head(5).items()])}

Parameters:
  Max Mismatches: {self.max_mismatches}
  Min Quality: {self.min_quality}

Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}
        """
        
        axes[2, 2].text(0.05, 0.95, summary_text, transform=axes[2, 2].transAxes,
                        fontsize=10, verticalalignment='top', fontfamily='monospace',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.5))
        
        plt.tight_layout()
        
        # Save plot
        plot_file = os.path.join(self.output_dirs['qc_reports'], 'phase3_quantification_qc.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        logger.info(f"✓ QC plots saved to: {plot_file}")
        
        plt.close()
    
    def save_summary_report(self, quantification_results, normalized_matrices):
        """Save comprehensive quantification report"""
        logger.info("=== Step 7: Generating Summary Report ===")
        
        # Calculate summary statistics
        total_detections = self.count_matrix.sum().sum()
        avg_detection_rate = np.mean([self.detection_stats[s]['detection_rate'] for s in self.detection_stats])
        detected_targets = (self.count_matrix.sum() > 0).sum()
        
        # Prepare detailed summary
        summary_data = {
            'pipeline_info': {
                'phase': 'Phase 3 - Chemokine Barcode Quantification',
                'date': datetime.now().isoformat(),
                'parameters': {
                    'max_mismatches': self.max_mismatches,
                    'min_quality': self.min_quality
                }
            },
            'quantification_results': {
                'total_samples': len(self.count_matrix),
                'total_targets': len(self.count_matrix.columns),
                'total_detections': int(total_detections),
                'average_detection_rate': avg_detection_rate,
                'detected_targets': int(detected_targets),
                'detection_efficiency': (detected_targets / len(self.count_matrix.columns)) * 100
            },
            'sample_statistics': dict(self.detection_stats),
            'target_statistics': {
                'total_counts_per_target': self.count_matrix.sum().to_dict(),
                'detection_rate_per_target': ((self.count_matrix > 0).sum() / len(self.count_matrix) * 100).to_dict()
            }
        }
        
        # Save JSON summary
        json_file = os.path.join(self.output_dirs['qc_reports'], 'phase3_summary.json')
        with open(json_file, 'w') as f:
            json.dump(summary_data, f, indent=2, default=str)
        logger.info(f"✓ JSON summary saved to: {json_file}")
        
        # Create human-readable report
        report_file = os.path.join(self.output_dirs['qc_reports'], 'phase3_report.txt')
        with open(report_file, 'w') as f:
            f.write("CHEMOKINE RECEPTOR DEMULTIPLEXING PIPELINE\n")
            f.write("Phase 3: Chemokine Barcode Quantification Results\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Max Mismatches: {self.max_mismatches}\n")
            f.write(f"Min Quality Score: {self.min_quality}\n\n")
            
            f.write("QUANTIFICATION OVERVIEW\n")
            f.write("-" * 30 + "\n")
            f.write(f"Total Samples Analyzed: {len(self.count_matrix)}\n")
            f.write(f"Total Chemokine Targets: {len(self.count_matrix.columns)}\n")
            f.write(f"Total Barcode Detections: {total_detections:,}\n")
            f.write(f"Average Detection Rate: {avg_detection_rate:.1f}%\n")
            f.write(f"Targets Successfully Detected: {detected_targets}/{len(self.count_matrix.columns)}\n")
            f.write(f"Target Detection Efficiency: {(detected_targets / len(self.count_matrix.columns)) * 100:.1f}%\n\n")
            
            f.write("SAMPLE PERFORMANCE\n")
            f.write("-" * 30 + "\n")
            for sample_name, stats in self.detection_stats.items():
                f.write(f"{sample_name}:\n")
                f.write(f"  Total reads: {stats['total_reads']:,}\n")
                f.write(f"  Reads with barcodes: {stats['reads_with_barcodes']:,}\n")
                f.write(f"  Detection rate: {stats['detection_rate']:.1f}%\n")
                f.write(f"  Average quality: {self.quality_stats[sample_name]['average_quality']:.1f}\n\n")
            
            f.write("TARGET ABUNDANCE RANKING\n")
            f.write("-" * 30 + "\n")
            target_totals = self.count_matrix.sum().sort_values(ascending=False)
            for i, (target, count) in enumerate(target_totals.items(), 1):
                detection_rate = ((self.count_matrix[target] > 0).sum() / len(self.count_matrix)) * 100
                f.write(f"{i:2d}. {target}: {count:,} total counts ({detection_rate:.0f}% samples)\n")
            
            f.write("\nOUTPUT FILES\n")
            f.write("-" * 30 + "\n")
            f.write("Count matrices:\n")
            f.write("  results/raw_count_matrix.csv\n")
            f.write("  results/cpm_normalized_counts.csv\n")
            f.write("  results/tpm_normalized_counts.csv\n")
            f.write("  results/log2_cpm_normalized_counts.csv\n")
            f.write("  results/relative_normalized_counts.csv\n\n")
            f.write("Quality control:\n")
            f.write("  qc_reports/phase3_quantification_qc.png\n")
            f.write("  qc_reports/phase3_summary.json\n")
            f.write("  logs/phase3_quantification.log\n\n")
            
            f.write("NEXT STEPS\n")
            f.write("-" * 30 + "\n")
            f.write("1. Review count matrix and detection efficiency\n")
            f.write("2. Proceed to Phase 4: Quality Control and Validation\n")
            f.write("3. Perform statistical analysis comparing GFP+ vs GFP- conditions\n")
            f.write("4. Validate biological relevance of detected chemokines\n")
        
        logger.info(f"✓ Human-readable report saved to: {report_file}")
    
    def run_pipeline(self):
        """Execute the complete Phase 3 pipeline"""
        logger.info("Starting Phase 3: Chemokine Barcode Quantification")
        
        # Create output directories
        for dir_name in self.output_dirs.values():
            os.makedirs(dir_name, exist_ok=True)
        
        try:
            # Execute pipeline steps
            self.load_target_barcodes()
            self.find_sample_files()
            quantification_results = self.quantify_all_samples()
            count_matrix = self.create_count_matrix(quantification_results)
            normalized_matrices = self.normalize_counts()
            self.create_qc_plots(quantification_results, normalized_matrices)
            self.save_summary_report(quantification_results, normalized_matrices)
            
            # Calculate final metrics
            total_detections = count_matrix.sum().sum()
            detected_targets = (count_matrix.sum() > 0).sum()
            detection_efficiency = (detected_targets / len(count_matrix.columns)) * 100
            
            logger.info("=" * 60)
            logger.info("Phase 3 completed successfully!")
            logger.info(f"Total barcode detections: {total_detections:,}")
            logger.info(f"Target detection efficiency: {detection_efficiency:.1f}%")
            logger.info("=" * 60)
            logger.info("Ready to proceed to Phase 4: Quality Control and Validation")
            
            return True, detection_efficiency
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False, 0

def main():
    """Main execution function"""
    print("Chemokine Receptor Demultiplexing Pipeline")
    print("Phase 3: Chemokine Barcode Quantification")
    print("=" * 60)
    
    # File paths
    barcode_file = "reference_files/chemokine_barcodes_with_rc.csv"
    sample_dir = "demultiplexed_samples"
    
    # Check if files exist
    if not os.path.exists(barcode_file):
        print(f"Error: Barcode file not found - {barcode_file}")
        print("Please ensure Phase 1 has been completed.")
        sys.exit(1)
    
    if not os.path.exists(sample_dir):
        print(f"Error: Sample directory not found - {sample_dir}")
        print("Please ensure Phase 2 has been completed.")
        sys.exit(1)
    
    # Initialize and run pipeline
    quantifier = BarcodeQuantifier(
        barcode_file=barcode_file,
        sample_dir=sample_dir,
        max_mismatches=2,  # Allow 2 mismatches for 20bp barcodes
        min_quality=20     # Minimum average quality score
    )
    
    success, efficiency = quantifier.run_pipeline()
    
    if success:
        print("\nPhase 3 completed successfully.")
        print(f"Target detection efficiency: {efficiency:.1f}%")
        print("\nGenerated files:")
        print("- results/raw_count_matrix.csv")
        print("- results/*_normalized_counts.csv (multiple normalization methods)")
        print("- qc_reports/phase3_quantification_qc.png")
        print("- qc_reports/phase3_report.txt")
        print("- qc_reports/phase3_summary.json")
        print("- logs/phase3_quantification.log")
        print("\nReady to proceed to Phase 4: Quality Control and Validation")
    else:
        print("\nPhase 3 failed. Check logs for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()