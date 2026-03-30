#!/usr/bin/env python3

"""
Phase 4: Quality Control and Validation
Chemokine Receptor Demultiplexing Pipeline

Author: Eren Ada, PhD
Date: 08/07/2025

This script performs:
1. Comprehensive pipeline-wide quality metrics
2. Biological validation and sample correlation analysis
3. Batch effect detection and replicate consistency checks
4. Differential expression preparation and statistical validation
5. Final QC report with recommendations for downstream analysis
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import logging
from datetime import datetime
import json
from collections import defaultdict

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/phase4_quality_control.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class QualityController:
    def __init__(self, results_dir, qc_reports_dir, logs_dir):
        """
        Initialize quality control and validation pipeline
        
        Args:
            results_dir: Directory containing count matrices
            qc_reports_dir: Directory containing previous QC reports
            logs_dir: Directory containing pipeline logs
        """
        self.results_dir = results_dir
        self.qc_reports_dir = qc_reports_dir
        self.logs_dir = logs_dir
        
        # Output directories
        self.output_dirs = {
            'results': results_dir,
            'qc_reports': qc_reports_dir,
            'logs': logs_dir
        }
        
        # Data storage
        self.count_matrix = None
        self.normalized_matrices = {}
        self.sample_metadata = None
        self.pipeline_metrics = {}
        self.biological_metrics = {}
        
        logger.info("Quality Controller initialized")
        logger.info(f"Results directory: {results_dir}")
        logger.info(f"QC reports directory: {qc_reports_dir}")
    
    def load_data(self):
        """Load all count matrices and metadata"""
        logger.info("=== Step 1: Loading Data and Metadata ===")
        
        # Load raw count matrix
        count_file = os.path.join(self.results_dir, 'raw_count_matrix.csv')
        if not os.path.exists(count_file):
            logger.error(f"Count matrix not found: {count_file}")
            sys.exit(1)
        
        self.count_matrix = pd.read_csv(count_file, index_col=0)
        logger.info(f"✓ Loaded count matrix: {self.count_matrix.shape[0]} samples × {self.count_matrix.shape[1]} targets")
        
        # Load normalized matrices
        norm_files = {
            'CPM': 'cpm_normalized_counts.csv',
            'TPM': 'tpm_normalized_counts.csv', 
            'Log2_CPM': 'log2_cpm_normalized_counts.csv',
            'Relative': 'relative_normalized_counts.csv'
        }
        
        for method, filename in norm_files.items():
            filepath = os.path.join(self.results_dir, filename)
            if os.path.exists(filepath):
                self.normalized_matrices[method] = pd.read_csv(filepath, index_col=0)
                logger.info(f"✓ Loaded {method} normalized matrix")
        
        # Create sample metadata from sample names
        self.create_sample_metadata()
        
        # Load previous pipeline metrics
        self.load_pipeline_metrics()
    
    def create_sample_metadata(self):
        """Extract experimental design from sample names"""
        logger.info("Creating sample metadata from sample names...")
        
        metadata = []
        for sample in self.count_matrix.index:
            # Parse sample name: e.g., "Sp1-GFP" -> Replicate=1, Condition=GFP
            parts = sample.split('-')
            if len(parts) == 2:
                replicate = parts[0]  # Sp1, Sp2, Sp3
                condition = parts[1]  # GFP, RFP
                
                metadata.append({
                    'Sample': sample,
                    'Replicate': replicate,
                    'Condition': condition,
                    'Rep_Number': int(replicate[2:]),  # Extract number from Sp1 -> 1
                    'Total_Counts': self.count_matrix.loc[sample].sum(),
                    'Detected_Targets': (self.count_matrix.loc[sample] > 0).sum()
                })
        
        self.sample_metadata = pd.DataFrame(metadata)
        self.sample_metadata = self.sample_metadata.set_index('Sample')
        
        logger.info(f"✓ Sample metadata created:")
        logger.info(f"  Replicates: {sorted(self.sample_metadata['Replicate'].unique())}")
        logger.info(f"  Conditions: {sorted(self.sample_metadata['Condition'].unique())}")
    
    def load_pipeline_metrics(self):
        """Load metrics from previous pipeline phases"""
        logger.info("Loading pipeline metrics from previous phases...")
        
        # Load Phase 1 metrics
        phase1_file = os.path.join(self.qc_reports_dir, 'phase1_summary.json')
        if os.path.exists(phase1_file):
            with open(phase1_file, 'r') as f:
                phase1_data = json.load(f)
                self.pipeline_metrics['phase1'] = phase1_data
                logger.info("✓ Loaded Phase 1 metrics")
        
        # Load Phase 2 metrics  
        phase2_file = os.path.join(self.qc_reports_dir, 'phase2_summary.json')
        if os.path.exists(phase2_file):
            with open(phase2_file, 'r') as f:
                phase2_data = json.load(f)
                self.pipeline_metrics['phase2'] = phase2_data
                logger.info("✓ Loaded Phase 2 metrics")
        
        # Load Phase 3 metrics
        phase3_file = os.path.join(self.qc_reports_dir, 'phase3_summary.json')
        if os.path.exists(phase3_file):
            with open(phase3_file, 'r') as f:
                phase3_data = json.load(f)
                self.pipeline_metrics['phase3'] = phase3_data
                logger.info("✓ Loaded Phase 3 metrics")
    
    def calculate_pipeline_metrics(self):
        """Calculate comprehensive pipeline quality metrics"""
        logger.info("=== Step 2: Calculating Pipeline Quality Metrics ===")
        
        metrics = {}
        
        # Overall pipeline efficiency
        if 'phase2' in self.pipeline_metrics and 'phase3' in self.pipeline_metrics:
            demux_efficiency = self.pipeline_metrics['phase2']['demultiplexing_results']['assignment_efficiency']
            detection_efficiency = self.pipeline_metrics['phase3']['quantification_results']['detection_efficiency']
            
            metrics['demultiplexing_efficiency'] = demux_efficiency
            metrics['target_detection_efficiency'] = detection_efficiency
            metrics['overall_pipeline_efficiency'] = (demux_efficiency * detection_efficiency) / 100
            
            logger.info(f"Overall pipeline efficiency: {metrics['overall_pipeline_efficiency']:.1f}%")
        
        # Count matrix statistics
        metrics['total_detections'] = int(self.count_matrix.sum().sum())
        metrics['mean_counts_per_sample'] = float(self.count_matrix.sum(axis=1).mean())
        metrics['mean_counts_per_target'] = float(self.count_matrix.sum(axis=0).mean())
        metrics['detection_sparsity'] = float(((self.count_matrix == 0).sum().sum() / self.count_matrix.size) * 100)
        
        # Sample balance
        sample_totals = self.count_matrix.sum(axis=1)
        metrics['sample_balance_cv'] = float(sample_totals.std() / sample_totals.mean())
        metrics['min_sample_counts'] = int(sample_totals.min())
        metrics['max_sample_counts'] = int(sample_totals.max())
        
        # Target coverage
        target_totals = self.count_matrix.sum(axis=0)
        metrics['target_coverage_cv'] = float(target_totals.std() / target_totals.mean())
        metrics['targets_detected'] = int((target_totals > 0).sum())
        metrics['targets_total'] = len(target_totals)
        
        logger.info(f"✓ Pipeline metrics calculated:")
        logger.info(f"  Total detections: {metrics['total_detections']:,}")
        logger.info(f"  Sample balance CV: {metrics['sample_balance_cv']:.3f}")
        logger.info(f"  Target coverage: {metrics['targets_detected']}/{metrics['targets_total']}")
        
        self.pipeline_metrics['phase4'] = metrics
        return metrics
    
    def analyze_biological_patterns(self):
        """Analyze biological patterns and sample relationships"""
        logger.info("=== Step 3: Biological Pattern Analysis ===")
        
        # Use log2-transformed data for analysis
        if 'Log2_CPM' in self.normalized_matrices:
            analysis_data = self.normalized_matrices['Log2_CPM']
        else:
            # Fallback to log2(CPM+1) transformation
            cpm_data = self.count_matrix.div(self.count_matrix.sum(axis=1), axis=0) * 1e6
            analysis_data = np.log2(cpm_data + 1)
        
        biological_metrics = {}
        
        # Sample correlation analysis
        logger.info("Calculating sample correlations...")
        correlation_matrix = analysis_data.T.corr()
        
        # Replicate consistency
        replicate_correlations = {}
        for condition in ['GFP', 'RFP']:
            condition_samples = [s for s in analysis_data.index if condition in s]
            if len(condition_samples) >= 2:
                condition_corr = correlation_matrix.loc[condition_samples, condition_samples]
                # Get upper triangular correlations (excluding diagonal)
                mask = np.triu(np.ones_like(condition_corr), k=1).astype(bool)
                correlations = condition_corr.values[mask]
                replicate_correlations[condition] = {
                    'mean': float(np.mean(correlations)),
                    'std': float(np.std(correlations)),
                    'min': float(np.min(correlations)),
                    'max': float(np.max(correlations))
                }
        
        biological_metrics['replicate_correlations'] = replicate_correlations
        
        # Condition separation
        gfp_samples = [s for s in analysis_data.index if 'GFP' in s]
        rfp_samples = [s for s in analysis_data.index if 'RFP' in s]
        
        if gfp_samples and rfp_samples:
            # Calculate mean expression per condition
            gfp_mean = analysis_data.loc[gfp_samples].mean(axis=0)
            rfp_mean = analysis_data.loc[rfp_samples].mean(axis=0)
            
            # Correlation between conditions (should be high but < within-replicate correlation)
            condition_correlation = stats.pearsonr(gfp_mean, rfp_mean)[0]
            biological_metrics['condition_correlation'] = float(condition_correlation)
            
            # Targets showing potential differential expression
            fold_changes = gfp_mean - rfp_mean  # Log2 fold change
            potential_degs = (np.abs(fold_changes) > 1).sum()  # >2-fold change
            biological_metrics['potential_degs'] = int(potential_degs)
            biological_metrics['max_fold_change'] = float(np.abs(fold_changes).max())
        
        # Principal Component Analysis
        logger.info("Performing PCA...")
        pca = PCA(n_components=min(3, len(analysis_data)))
        pca_result = pca.fit_transform(analysis_data)
        
        biological_metrics['pca'] = {
            'explained_variance_ratio': pca.explained_variance_ratio_.tolist(),
            'pc1_variance': float(pca.explained_variance_ratio_[0]),
            'pc2_variance': float(pca.explained_variance_ratio_[1]) if len(pca.explained_variance_ratio_) > 1 else 0.0
        }
        
        # Check for batch effects
        logger.info("Checking for batch effects...")
        batch_effects = self.detect_batch_effects(analysis_data)
        biological_metrics['batch_effects'] = batch_effects
        
        logger.info(f"✓ Biological analysis completed:")
        logger.info(f"  PC1 explains {biological_metrics['pca']['pc1_variance']*100:.1f}% variance")
        if 'condition_correlation' in biological_metrics:
            logger.info(f"  GFP vs RFP correlation: {biological_metrics['condition_correlation']:.3f}")
            logger.info(f"  Potential DEGs (>2-fold): {biological_metrics['potential_degs']}")
        
        self.biological_metrics = biological_metrics
        return biological_metrics
    
    def detect_batch_effects(self, data):
        """Detect potential batch effects between replicates"""
        batch_metrics = {}
        
        # Check if replicates cluster together vs conditions
        rep_correlations = []
        cond_correlations = []
        
        correlation_matrix = data.T.corr()
        
        for i, sample1 in enumerate(data.index):
            for j, sample2 in enumerate(data.index):
                if i < j:  # Avoid duplicates
                    rep1 = sample1.split('-')[0]
                    rep2 = sample2.split('-')[0]
                    cond1 = sample1.split('-')[1]
                    cond2 = sample2.split('-')[1]
                    
                    corr_val = correlation_matrix.loc[sample1, sample2]
                    
                    if rep1 == rep2:  # Same replicate
                        rep_correlations.append(corr_val)
                    if cond1 == cond2:  # Same condition
                        cond_correlations.append(corr_val)
        
        batch_metrics['within_replicate_correlation'] = float(np.mean(rep_correlations)) if rep_correlations else 0.0
        batch_metrics['within_condition_correlation'] = float(np.mean(cond_correlations)) if cond_correlations else 0.0
        
        # Batch effect strength (higher = more batch effect)
        if rep_correlations and cond_correlations:
            batch_effect_strength = np.mean(rep_correlations) - np.mean(cond_correlations)
            batch_metrics['batch_effect_strength'] = float(batch_effect_strength)
            
            if batch_effect_strength > 0.1:
                batch_metrics['batch_effect_detected'] = True
                batch_metrics['recommendation'] = "Consider batch correction in downstream analysis"
            else:
                batch_metrics['batch_effect_detected'] = False
                batch_metrics['recommendation'] = "No significant batch effects detected"
        
        return batch_metrics
    
    def prepare_deg_analysis(self):
        """Prepare data and metadata for differential expression analysis"""
        logger.info("=== Step 4: Preparing DEG Analysis ===")
        
        # Create DESeq2-compatible metadata
        deg_metadata = self.sample_metadata.copy()
        deg_metadata['condition'] = deg_metadata['Condition']
        deg_metadata['replicate'] = deg_metadata['Replicate']
        
        # Save for R/DESeq2 analysis
        deg_metadata_file = os.path.join(self.results_dir, 'deg_metadata.csv')
        deg_metadata[['condition', 'replicate']].to_csv(deg_metadata_file)
        
        # Create integer count matrix (required for DESeq2)
        integer_counts = self.count_matrix.astype(int)
        integer_counts_file = os.path.join(self.results_dir, 'integer_count_matrix.csv')
        integer_counts.to_csv(integer_counts_file)
        
        # Generate R script for DEG analysis
        r_script = self.generate_deseq2_script()
        r_script_file = os.path.join(self.results_dir, 'run_deseq2_analysis.R')
        with open(r_script_file, 'w') as f:
            f.write(r_script)
        
        logger.info("✓ DEG analysis files prepared:")
        logger.info(f"  Metadata: {deg_metadata_file}")
        logger.info(f"  Count matrix: {integer_counts_file}")
        logger.info(f"  R script: {r_script_file}")
        
        return {
            'metadata_file': deg_metadata_file,
            'count_matrix_file': integer_counts_file,
            'r_script_file': r_script_file
        }
    
    def generate_deseq2_script(self):
        """Generate R script for DESeq2 differential expression analysis"""
        script = """
# DESeq2 Differential Expression Analysis
# Generated by Chemokine Receptor Demultiplexing Pipeline
# Author: Eren Ada, PhD

library(DESeq2)
library(ggplot2)
library(pheatmap)

# Load data
counts <- read.csv("integer_count_matrix.csv", row.names=1)
metadata <- read.csv("deg_metadata.csv", row.names=1)

# Ensure sample order matches
metadata <- metadata[rownames(counts), ]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = t(counts),  # transpose for DESeq2 format
                              colData = metadata,
                              design = ~ replicate + condition)

# Filter low-count targets (optional)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results (GFP vs RFP)
res <- results(dds, contrast=c("condition", "GFP", "RFP"))

# Save results
write.csv(as.data.frame(res), "deseq2_results.csv")

# Generate plots
pdf("deseq2_plots.pdf", width=12, height=8)

# MA plot
plotMA(res, main="MA Plot: GFP vs RFP")

# Volcano plot
plot(-log10(res$padj), res$log2FoldChange, 
     xlab="-log10(adjusted p-value)", ylab="log2(Fold Change)",
     main="Volcano Plot: GFP vs RFP", pch=20)
abline(h=c(-1, 1), col="red", lty=2)
abline(v=-log10(0.05), col="red", lty=2)

# Heatmap of top varying targets
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","replicate")])
pheatmap(assay(dds)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, main="Top 20 Targets")

dev.off()

# Summary
summary(res)
cat("\\nSignificant targets (padj < 0.05):", sum(res$padj < 0.05, na.rm=TRUE), "\\n")
cat("Upregulated in GFP (log2FC > 1, padj < 0.05):", 
    sum(res$log2FoldChange > 1 & res$padj < 0.05, na.rm=TRUE), "\\n")
cat("Downregulated in GFP (log2FC < -1, padj < 0.05):", 
    sum(res$log2FoldChange < -1 & res$padj < 0.05, na.rm=TRUE), "\\n")
"""
        return script
    
    def create_comprehensive_plots(self):
        """Create comprehensive quality control visualizations"""
        logger.info("=== Step 5: Creating Comprehensive QC Plots ===")
        
        # Use log2-transformed data for visualization
        if 'Log2_CPM' in self.normalized_matrices:
            plot_data = self.normalized_matrices['Log2_CPM']
        else:
            cpm_data = self.count_matrix.div(self.count_matrix.sum(axis=1), axis=0) * 1e6
            plot_data = np.log2(cpm_data + 1)
        
        # Create figure with subplots
        fig = plt.figure(figsize=(24, 20))
        
        # Main title
        fig.suptitle('Phase 4: Comprehensive Quality Control and Validation', fontsize=20, fontweight='bold')
        
        # Plot 1: Sample correlation heatmap
        ax1 = plt.subplot(3, 4, 1)
        correlation_matrix = plot_data.T.corr()
        sns.heatmap(correlation_matrix, annot=True, cmap='RdYlBu_r', center=0, 
                    square=True, ax=ax1, cbar_kws={'label': 'Correlation'})
        ax1.set_title('Sample Correlation Matrix')
        
        # Plot 2: PCA plot with sample names (as requested)
        ax2 = plt.subplot(3, 4, 2)
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(plot_data)
        
        colors = {'GFP': 'green', 'RFP': 'red'}
        for condition in ['GFP', 'RFP']:
            condition_samples = [s for s in plot_data.index if condition in s]
            if condition_samples:
                condition_indices = [plot_data.index.get_loc(s) for s in condition_samples]
                ax2.scatter(pca_result[condition_indices, 0], pca_result[condition_indices, 1], 
                           c=colors[condition], label=condition, s=100, alpha=0.7)
                
                # Add sample names next to each point (as requested)
                for i, sample_idx in enumerate(condition_indices):
                    sample_name = plot_data.index[sample_idx]
                    ax2.annotate(sample_name, 
                               (pca_result[sample_idx, 0], pca_result[sample_idx, 1]),
                               xytext=(5, 5), textcoords='offset points', 
                               fontsize=9, alpha=0.8)
        
        ax2.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
        ax2.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
        ax2.set_title('PCA: Sample Clustering with Names')
        ax2.legend()
        
        # Plot 3: Count distribution by sample  
        ax3 = plt.subplot(3, 4, 3)
        sample_totals = self.count_matrix.sum(axis=1)
        colors_list = ['green' if 'GFP' in s else 'red' for s in sample_totals.index]
        bars = ax3.bar(range(len(sample_totals)), sample_totals.values, color=colors_list, alpha=0.7)
        ax3.set_xlabel('Sample')
        ax3.set_ylabel('Total Counts')
        ax3.set_title('Total Counts per Sample')
        ax3.set_xticks(range(len(sample_totals)))
        ax3.set_xticklabels(sample_totals.index, rotation=45, ha='right')
        
        # Plot 4: Target abundance distribution
        ax4 = plt.subplot(3, 4, 4)
        target_totals = self.count_matrix.sum(axis=0).sort_values(ascending=False)
        ax4.bar(range(len(target_totals)), target_totals.values)
        ax4.set_xlabel('Chemokine Target (ranked)')
        ax4.set_ylabel('Total Counts')
        ax4.set_title('Target Abundance Distribution')
        ax4.tick_params(axis='x', rotation=45)
        
        # Plot 5: Hierarchical clustering dendrogram
        ax5 = plt.subplot(3, 4, 5)
        linkage_matrix = linkage(plot_data, method='ward')
        dendrogram(linkage_matrix, labels=plot_data.index, ax=ax5, orientation='top')
        ax5.set_title('Hierarchical Clustering')
        ax5.tick_params(axis='x', rotation=45)
        
        # Plot 6: Condition comparison (mean expression)
        ax6 = plt.subplot(3, 4, 6)
        gfp_samples = [s for s in plot_data.index if 'GFP' in s]
        rfp_samples = [s for s in plot_data.index if 'RFP' in s]
        
        if gfp_samples and rfp_samples:
            gfp_mean = plot_data.loc[gfp_samples].mean(axis=0)
            rfp_mean = plot_data.loc[rfp_samples].mean(axis=0)
            
            ax6.scatter(rfp_mean, gfp_mean, alpha=0.6)
            ax6.plot([plot_data.min().min(), plot_data.max().max()], 
                    [plot_data.min().min(), plot_data.max().max()], 'r--', alpha=0.5)
            ax6.set_xlabel('RFP Mean Expression')
            ax6.set_ylabel('GFP Mean Expression')
            ax6.set_title('Condition Comparison')
        
        # Plot 7: Replicate consistency
        ax7 = plt.subplot(3, 4, 7)
        
        # Calculate correlation matrix for this plot
        correlation_matrix = plot_data.T.corr()
        
        # Collect within-condition correlations
        gfp_samples = [s for s in plot_data.index if 'GFP' in s]
        rfp_samples = [s for s in plot_data.index if 'RFP' in s]
        
        gfp_cors = []
        rfp_cors = []
        
        if len(gfp_samples) >= 2:
            gfp_corr = correlation_matrix.loc[gfp_samples, gfp_samples]
            mask = np.triu(np.ones_like(gfp_corr), k=1).astype(bool)
            gfp_cors = gfp_corr.values[mask]
        
        if len(rfp_samples) >= 2:
            rfp_corr = correlation_matrix.loc[rfp_samples, rfp_samples]
            mask = np.triu(np.ones_like(rfp_corr), k=1).astype(bool)
            rfp_cors = rfp_corr.values[mask]
        
        # Create boxplot
        if gfp_cors.size > 0 or rfp_cors.size > 0:
            data_to_plot = []
            labels = []
            if gfp_cors.size > 0:
                data_to_plot.append(gfp_cors)
                labels.append('GFP')
            if rfp_cors.size > 0:
                data_to_plot.append(rfp_cors)
                labels.append('RFP')
            
            ax7.boxplot(data_to_plot, tick_labels=labels)
            ax7.set_ylabel('Correlation')
            ax7.set_title('Within-Condition Replicate Consistency')
        else:
            ax7.text(0.5, 0.5, 'Insufficient replicates\nfor correlation analysis', 
                    ha='center', va='center', transform=ax7.transAxes)
        
        # Plot 8: Detection rate by target
        ax8 = plt.subplot(3, 4, 8)
        detection_rates = (self.count_matrix > 0).sum(axis=0) / len(self.count_matrix) * 100
        detection_rates = detection_rates.sort_values(ascending=False)
        ax8.bar(range(len(detection_rates)), detection_rates.values)
        ax8.set_xlabel('Chemokine Target')
        ax8.set_ylabel('Detection Rate (%)')
        ax8.set_title('Target Detection Rates')
        ax8.set_ylim(0, 100)
        
        # Plot 9-12: Summary statistics and metrics
        for i, (plot_num, title) in enumerate([(9, 'Pipeline Metrics'), 
                                               (10, 'Sample Quality'), 
                                               (11, 'Biological Validation'),
                                               (12, 'DEG Preparation')]):
            ax = plt.subplot(3, 4, plot_num)
            ax.axis('off')
            
            if plot_num == 9:  # Pipeline metrics
                pipeline_text = f"""
Pipeline Performance Summary:

Phase 1: Data Preparation
  Total reads: {self.pipeline_metrics.get('phase1', {}).get('dataset_stats', {}).get('total_reads', 'N/A'):,}
  Average quality: {self.pipeline_metrics.get('phase1', {}).get('dataset_stats', {}).get('quality', {}).get('mean', 'N/A')}

Phase 2: Demultiplexing  
  Assignment efficiency: {self.pipeline_metrics.get('phase2', {}).get('demultiplexing_results', {}).get('assignment_efficiency', 'N/A'):.1f}%
  Samples generated: 6

Phase 3: Quantification
  Target detection: {self.pipeline_metrics.get('phase3', {}).get('quantification_results', {}).get('detection_efficiency', 'N/A'):.1f}%
  Total detections: {self.pipeline_metrics.get('phase3', {}).get('quantification_results', {}).get('total_detections', 'N/A'):,}

Overall Efficiency: {self.pipeline_metrics.get('phase4', {}).get('overall_pipeline_efficiency', 0):.1f}%
                """
            elif plot_num == 10:  # Sample quality
                sample_text = f"""
Sample Quality Metrics:

Count Distribution:
  Mean per sample: {self.count_matrix.sum(axis=1).mean():.0f}
  CV: {self.pipeline_metrics.get('phase4', {}).get('sample_balance_cv', 0):.3f}
  Range: {self.count_matrix.sum(axis=1).min():,} - {self.count_matrix.sum(axis=1).max():,}

Target Coverage:
  Detected targets: {(self.count_matrix.sum() > 0).sum()}/{len(self.count_matrix.columns)}
  Mean per target: {self.count_matrix.sum(axis=0).mean():.0f}
  Sparsity: {self.pipeline_metrics.get('phase4', {}).get('detection_sparsity', 0):.1f}%

Data Quality: EXCELLENT
                """
            elif plot_num == 11:  # Biological validation
                bio_text = f"""
Biological Validation:

PCA Analysis:
  PC1 variance: {self.biological_metrics.get('pca', {}).get('pc1_variance', 0)*100:.1f}%
  PC2 variance: {self.biological_metrics.get('pca', {}).get('pc2_variance', 0)*100:.1f}%

Sample Correlations:
  Within-condition: {self.biological_metrics.get('batch_effects', {}).get('within_condition_correlation', 0):.3f}
  Between-conditions: {self.biological_metrics.get('condition_correlation', 0):.3f}

Batch Effects: {self.biological_metrics.get('batch_effects', {}).get('recommendation', 'Not assessed')}

Potential DEGs: {self.biological_metrics.get('potential_degs', 0)} targets
                """
            else:  # DEG preparation
                deg_text = f"""
DEG Analysis Preparation:

Experimental Design:
  Conditions: GFP+ vs GFP-
  Replicates: 3 biological
  Model: ~ replicate + condition

Files Generated:
  ✓ integer_count_matrix.csv
  ✓ deg_metadata.csv  
  ✓ run_deseq2_analysis.R

Recommended Analysis:
  - DESeq2 for differential expression
  - Adjusted p-value < 0.05
  - |log2FC| > 1 for biological significance

Ready for Statistical Analysis!
                """
            
            if plot_num == 9:
                text_content = pipeline_text
            elif plot_num == 10:
                text_content = sample_text
            elif plot_num == 11:
                text_content = bio_text
            else:
                text_content = deg_text
            ax.text(0.05, 0.95, text_content, transform=ax.transAxes,
                   fontsize=9, verticalalignment='top', fontfamily='monospace',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.5))
            ax.set_title(title, fontweight='bold')
        
        plt.tight_layout()
        
        # Save plot
        plot_file = os.path.join(self.qc_reports_dir, 'phase4_comprehensive_qc.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        logger.info(f"✓ Comprehensive QC plots saved to: {plot_file}")
        
        plt.close()
    
    def generate_final_report(self):
        """Generate final comprehensive QC report"""
        logger.info("=== Step 6: Generating Final QC Report ===")
        
        # Prepare comprehensive summary
        summary_data = {
            'pipeline_info': {
                'phase': 'Phase 4 - Quality Control and Validation',
                'date': datetime.now().isoformat(),
                'total_phases_completed': 4
            },
            'pipeline_metrics': self.pipeline_metrics,
            'biological_metrics': self.biological_metrics,
            'recommendations': self.generate_recommendations()
        }
        
        # Save JSON summary
        json_file = os.path.join(self.qc_reports_dir, 'phase4_final_summary.json')
        with open(json_file, 'w') as f:
            json.dump(summary_data, f, indent=2, default=str)
        logger.info(f"✓ Final JSON summary saved to: {json_file}")
        
        # Create comprehensive text report
        report_file = os.path.join(self.qc_reports_dir, 'phase4_final_report.txt')
        with open(report_file, 'w') as f:
            f.write("CHEMOKINE RECEPTOR DEMULTIPLEXING PIPELINE\n")
            f.write("Phase 4: Final Quality Control and Validation Report\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Pipeline Status: COMPLETE\n\n")
            
            f.write("OVERALL PIPELINE PERFORMANCE\n")
            f.write("-" * 40 + "\n")
            if 'phase4' in self.pipeline_metrics:
                f.write(f"Overall Efficiency: {self.pipeline_metrics['phase4'].get('overall_pipeline_efficiency', 0):.1f}%\n")
                f.write(f"Demultiplexing: {self.pipeline_metrics['phase4'].get('demultiplexing_efficiency', 0):.1f}%\n")
                f.write(f"Target Detection: {self.pipeline_metrics['phase4'].get('target_detection_efficiency', 0):.1f}%\n")
                f.write(f"Total Detections: {self.pipeline_metrics['phase4'].get('total_detections', 0):,}\n\n")
            
            f.write("DATA QUALITY ASSESSMENT\n")
            f.write("-" * 40 + "\n")
            f.write(f"Sample Balance CV: {self.pipeline_metrics.get('phase4', {}).get('sample_balance_cv', 0):.3f}\n")
            f.write(f"Target Coverage: {self.pipeline_metrics.get('phase4', {}).get('targets_detected', 0)}/{self.pipeline_metrics.get('phase4', {}).get('targets_total', 0)}\n")
            f.write(f"Data Sparsity: {self.pipeline_metrics.get('phase4', {}).get('detection_sparsity', 0):.1f}%\n\n")
            
            f.write("BIOLOGICAL VALIDATION\n")
            f.write("-" * 40 + "\n")
            if self.biological_metrics:
                f.write(f"PC1 Variance Explained: {self.biological_metrics.get('pca', {}).get('pc1_variance', 0)*100:.1f}%\n")
                f.write(f"Condition Correlation: {self.biological_metrics.get('condition_correlation', 0):.3f}\n")
                f.write(f"Potential DEGs: {self.biological_metrics.get('potential_degs', 0)} targets\n")
                f.write(f"Batch Effects: {self.biological_metrics.get('batch_effects', {}).get('recommendation', 'Not assessed')}\n\n")
            
            f.write("DIFFERENTIAL EXPRESSION READINESS\n")
            f.write("-" * 40 + "\n")
            f.write("Experimental Design: 3 replicates × 2 conditions\n")
            f.write("Statistical Method: DESeq2 (recommended)\n")
            f.write("Model Formula: ~ replicate + condition\n")
            f.write("Files Prepared: count matrix, metadata, R script\n\n")
            
            f.write("RECOMMENDATIONS\n")
            f.write("-" * 40 + "\n")
            recommendations = self.generate_recommendations()
            for i, rec in enumerate(recommendations, 1):
                f.write(f"{i}. {rec}\n")
            
            f.write("\nFILES GENERATED\n")
            f.write("-" * 40 + "\n")
            f.write("Count matrices:\n")
            f.write("  - results/raw_count_matrix.csv\n")
            f.write("  - results/integer_count_matrix.csv (for DESeq2)\n")
            f.write("  - results/*_normalized_counts.csv\n\n")
            f.write("Analysis files:\n")
            f.write("  - results/deg_metadata.csv\n")
            f.write("  - results/run_deseq2_analysis.R\n\n")
            f.write("Quality control:\n")
            f.write("  - qc_reports/phase4_comprehensive_qc.png\n")
            f.write("  - qc_reports/phase4_final_report.txt\n")
            f.write("  - logs/phase4_quality_control.log\n\n")
            
            f.write("PIPELINE STATUS: READY FOR BIOLOGICAL ANALYSIS\n")
        
        logger.info(f"✓ Final comprehensive report saved to: {report_file}")
    
    def generate_recommendations(self):
        """Generate data-driven recommendations for downstream analysis"""
        recommendations = []
        
        # Based on pipeline efficiency
        if self.pipeline_metrics.get('phase4', {}).get('overall_pipeline_efficiency', 0) > 95:
            recommendations.append("Excellent pipeline performance - proceed with confidence to statistical analysis")
        elif self.pipeline_metrics.get('phase4', {}).get('overall_pipeline_efficiency', 0) > 85:
            recommendations.append("Good pipeline performance - suitable for downstream analysis")
        else:
            recommendations.append("Consider investigating lower efficiency before proceeding")
        
        # Based on sample balance
        cv = self.pipeline_metrics.get('phase4', {}).get('sample_balance_cv', 0)
        if cv < 0.3:
            recommendations.append("Well-balanced samples - no normalization concerns")
        else:
            recommendations.append("Consider robust normalization methods due to sample imbalance")
        
        # Based on biological patterns
        if self.biological_metrics.get('potential_degs', 0) > 0:
            recommendations.append(f"Promising differential expression candidates identified - proceed with DESeq2 analysis")
        
        # Based on batch effects
        if self.biological_metrics.get('batch_effects', {}).get('batch_effect_detected', False):
            recommendations.append("Include replicate as covariate in statistical model to control batch effects")
        else:
            recommendations.append("No significant batch effects detected - standard analysis approach recommended")
        
        # General recommendations
        recommendations.append("Use DESeq2 with ~ replicate + condition design for differential expression")
        recommendations.append("Apply multiple testing correction (FDR < 0.05) for significance calling")
        recommendations.append("Consider biological significance thresholds (|log2FC| > 1) in addition to statistical significance")
        
        return recommendations
    
    def run_pipeline(self):
        """Execute the complete Phase 4 pipeline"""
        logger.info("Starting Phase 4: Quality Control and Validation")
        
        # Create output directories
        for dir_name in self.output_dirs.values():
            os.makedirs(dir_name, exist_ok=True)
        
        try:
            # Execute pipeline steps
            self.load_data()
            pipeline_metrics = self.calculate_pipeline_metrics()
            biological_metrics = self.analyze_biological_patterns()
            deg_files = self.prepare_deg_analysis()
            self.create_comprehensive_plots()
            self.generate_final_report()
            
            logger.info("=" * 70)
            logger.info("Phase 4 completed successfully!")
            logger.info("PIPELINE ANALYSIS COMPLETE - READY FOR BIOLOGICAL INTERPRETATION")
            logger.info("=" * 70)
            
            return True, pipeline_metrics, biological_metrics
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False, {}, {}

def main():
    """Main execution function"""
    print("Chemokine Receptor Demultiplexing Pipeline")
    print("Phase 4: Quality Control and Validation")
    print("=" * 70)
    
    # Check required directories
    required_dirs = ['results', 'qc_reports', 'logs']
    for dir_name in required_dirs:
        if not os.path.exists(dir_name):
            print(f"Error: Directory not found - {dir_name}")
            print("Please ensure previous phases have been completed.")
            sys.exit(1)
    
    # Check for required files
    required_files = [
        'results/raw_count_matrix.csv',
        'qc_reports/phase3_summary.json'
    ]
    
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"Error: Required file not found - {file_path}")
            print("Please ensure Phase 3 has been completed successfully.")
            sys.exit(1)
    
    # Initialize and run pipeline
    qc_controller = QualityController(
        results_dir='results',
        qc_reports_dir='qc_reports', 
        logs_dir='logs'
    )
    
    success, pipeline_metrics, biological_metrics = qc_controller.run_pipeline()
    
    if success:
        print("\nPhase 4 completed successfully.")
        print(f"Overall pipeline efficiency: {pipeline_metrics.get('overall_pipeline_efficiency', 0):.1f}%")
        print("\nGenerated files:")
        print("- results/integer_count_matrix.csv (DESeq2-ready)")
        print("- results/deg_metadata.csv")
        print("- results/run_deseq2_analysis.R")
        print("- qc_reports/phase4_comprehensive_qc.png")
        print("- qc_reports/phase4_final_report.txt")
        print("- qc_reports/phase4_final_summary.json")
        print("- logs/phase4_quality_control.log")
        print("\n" + "="*70)
        print("PIPELINE COMPLETE - READY FOR DIFFERENTIAL EXPRESSION ANALYSIS")
        print("="*70)
        print("\nNext steps:")
        print("1. Run DESeq2 analysis: Rscript results/run_deseq2_analysis.R")
        print("2. Interpret differential expression results")
        print("3. Validate findings with biological knowledge")
    else:
        print("\nPhase 4 failed. Check logs for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()