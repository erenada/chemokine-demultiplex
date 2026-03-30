#!/usr/bin/env python3

"""
Adapted Chemokine Receptor Pipeline
Author: Eren Ada, PhD (Adapted by AI Assistant)
Date: 02/10/2026

This script adapts the original 4-phase pipeline for the new dataset structure:
1.  Reference Prep: Adapts 'Gene' column to 'Chemokine' and generates RCs.
2.  Quantification: Processes R1 FASTQ files directly (no demultiplexing needed).
3.  QC: Generates reports based on new metadata (Source, Condition, Replicate).
"""

import os
import sys
import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
import edlib
import logging
from datetime import datetime
import json
from collections import defaultdict
from tqdm import tqdm
import glob
# from scipy import stats
# from scipy.cluster.hierarchy import linkage, dendrogram
# from sklearn.decomposition import PCA

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# --- Configuration ---
BARCODE_INPUT_FILE = "/Users/eren/Desktop/HMS/EricLoo/project_miseq/barcodes/barcode_list.csv"
SEQ_DATA_DIR = "/Users/eren/Desktop/HMS/EricLoo/project_miseq/seq_data"
OUTPUT_BASE_DIR = "/Users/eren/Desktop/HMS/EricLoo/project_miseq"

# --- Step 1: Reference Preparation ---

class ReferencePrep:
    def __init__(self):
        self.output_dir = os.path.join(OUTPUT_BASE_DIR, 'reference_files')
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Setup file logging after output dir is created
        log_dir = os.path.join(OUTPUT_BASE_DIR, 'logs')
        os.makedirs(log_dir, exist_ok=True)
        file_handler = logging.FileHandler(os.path.join(log_dir, 'adapted_pipeline.log'))
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
        
    def run(self):
        logger.info("=== Step 1: Reference Preparation ===")
        
        if not os.path.exists(BARCODE_INPUT_FILE):
            logger.error(f"Barcode file not found: {BARCODE_INPUT_FILE}")
            sys.exit(1)
            
        try:
            # Load original file
            df = pd.read_csv(BARCODE_INPUT_FILE)
            logger.info(f"Loaded {len(df)} barcodes from {BARCODE_INPUT_FILE}")
            
            # Rename 'Gene' to 'Chemokine' if needed
            if 'Gene' in df.columns:
                df = df.rename(columns={'Gene': 'Chemokine'})
                logger.info("Renamed 'Gene' column to 'Chemokine'")
            
            if 'Chemokine' not in df.columns or 'Barcode' not in df.columns:
                logger.error("Required columns 'Chemokine' (or 'Gene') and 'Barcode' not found")
                sys.exit(1)
                
            # Generate Reverse Complements
            df['Barcode_rc'] = df['Barcode'].apply(lambda x: str(Seq(x).reverse_complement()))
            
            # Save
            output_file = os.path.join(self.output_dir, 'chemokine_barcodes_with_rc.csv')
            df.to_csv(output_file, index=False)
            logger.info(f"✓ Saved reference file with RCs to: {output_file}")
            return output_file
            
        except Exception as e:
            logger.error(f"Reference prep failed: {e}")
            sys.exit(1)

# --- Step 2: Quantification ---

class BarcodeQuantifier:
    def __init__(self, barcode_file, max_mismatches=2, min_quality=20):
        self.barcode_file = barcode_file
        self.max_mismatches = max_mismatches
        self.min_quality = min_quality
        self.output_dirs = {
            'results': os.path.join(OUTPUT_BASE_DIR, 'results'),
            'qc_reports': os.path.join(OUTPUT_BASE_DIR, 'qc_reports')
        }
        for d in self.output_dirs.values():
            os.makedirs(d, exist_ok=True)
            
        self.barcode_patterns = {}
        self.count_matrix = None
        self.detection_stats = defaultdict(dict)
        
    def load_barcodes(self):
        df = pd.read_csv(self.barcode_file)
        for _, row in df.iterrows():
            self.barcode_patterns[row['Chemokine']] = {
                'forward': row['Barcode'],
                'reverse': row['Barcode_rc']
            }
        logger.info(f"Loaded {len(self.barcode_patterns)} target patterns")

    def get_sample_name(self, filename):
        # Extract "Invitro-GPCR-Neg-1" from "LIB067987_GEN00315253_Invitro-GPCR-Neg-1_S1_L001_R1_001.fastq"
        # Strategy: Split by underscore, take the 3rd element (index 2)
        try:
            parts = os.path.basename(filename).split('_')
            if len(parts) >= 3:
                return parts[2]
            return os.path.basename(filename).split('.')[0]
        except:
            return os.path.basename(filename).split('.')[0]

    def process_file(self, filepath):
        sample_name = self.get_sample_name(filepath)
        logger.info(f"Processing sample: {sample_name} (File: {os.path.basename(filepath)})")
        
        counts = defaultdict(int)
        
        # Count total reads first for progress bar (optional, can be skipped for speed)
        # total_reads = sum(1 for _ in SeqIO.parse(filepath, "fastq"))
        # Using a simpler approach for progress to avoid double reading
        
        processed = 0
        detected = 0
        
        # Use SeqIO instead of pyfastx
        with open(filepath, "r") as handle:
            for read in tqdm(SeqIO.parse(handle, "fastq"), desc=sample_name, leave=False):
                processed += 1
                seq = str(read.seq)
                
                # Simple quality check (skip if avg quality < threshold)
                # Biopython stores quality as integers in letter_annotations
                qual = read.letter_annotations["phred_quality"]
                if np.mean(qual) < self.min_quality:
                    continue
                    
                # Find matches
                best_match = None
                best_dist = self.max_mismatches + 1
                
                # Optimization: Could use Aho-Corasick or similar for speed, 
                # but sticking to original fuzzy logic for consistency
                
                found_in_read = False
                for target, patterns in self.barcode_patterns.items():
                    # Check forward
                    res_f = edlib.align(patterns['forward'], seq, mode="HW", task="locations", k=self.max_mismatches)
                    if res_f['editDistance'] != -1:
                        if res_f['editDistance'] < best_dist:
                            best_match = target
                            best_dist = res_f['editDistance']
                            found_in_read = True
                    
                    # Check reverse (if not found better already)
                    if not found_in_read:
                        res_r = edlib.align(patterns['reverse'], seq, mode="HW", task="locations", k=self.max_mismatches)
                        if res_r['editDistance'] != -1:
                            if res_r['editDistance'] < best_dist:
                                best_match = target
                                best_dist = res_r['editDistance']
                                found_in_read = True
                
                if best_match:
                    counts[best_match] += 1
                    detected += 1
                
        self.detection_stats[sample_name] = {
            'total_reads': processed,
            'detected_reads': detected,
            'rate': (detected/processed*100) if processed > 0 else 0
        }
        return sample_name, counts

    def run(self, read_type="R1"):
        logger.info(f"=== Step 2: Quantification ({read_type}) ===")
        self.load_barcodes()
        
        # Find files based on read type
        pattern = os.path.join(SEQ_DATA_DIR, f"*_{read_type}_*.fastq")
        files = sorted(glob.glob(pattern))
        
        if not files:
            logger.error(f"No {read_type} FASTQ files found in {SEQ_DATA_DIR}")
            return None
            
        logger.info(f"Found {len(files)} {read_type} files to process")
        
        all_counts = {}
        for f in files:
            name, counts = self.process_file(f)
            all_counts[name] = counts
            
        # Create Matrix
        self.count_matrix = pd.DataFrame(all_counts).fillna(0).astype(int)
        # Sort index (targets) and columns (samples)
        self.count_matrix = self.count_matrix.sort_index().sort_index(axis=1)
        
        # Save Matrix to specific subfolder
        results_subdir = os.path.join(self.output_dirs['results'], read_type)
        os.makedirs(results_subdir, exist_ok=True)
        
        out_file = os.path.join(results_subdir, 'raw_count_matrix.csv')
        self.count_matrix.to_csv(out_file)
        logger.info(f"✓ Saved {read_type} count matrix to: {out_file}")
        
        return self.count_matrix

# --- Step 3: Quality Control ---

class QualityControl:
    def __init__(self, count_matrix, read_type="R1"):
        self.count_matrix = count_matrix
        self.read_type = read_type
        self.output_dir = os.path.join(OUTPUT_BASE_DIR, 'qc_reports', read_type)
        os.makedirs(self.output_dir, exist_ok=True)
        
    def parse_metadata(self):
        # Parse: Invitro-GPCR-Neg-1
        # Source: Invitro
        # Condition: Neg
        # Replicate: 1
        meta = []
        for sample in self.count_matrix.columns:
            parts = sample.split('-')
            # Expected format: Source-GPCR-Condition-Replicate
            # e.g. Invitro-GPCR-Neg-1
            # e.g. MC38-GPCR-POS-1
            # e.g. Spleen-GPCR-Neg-1
            
            if len(parts) >= 4:
                source = parts[0]
                condition = parts[2]
                replicate = parts[3]
                group = f"{source}-{condition}"
            else:
                # Fallback
                source = "Unknown"
                condition = "Unknown"
                replicate = "Unknown"
                group = "Unknown"
                
            meta.append({
                'Sample': sample,
                'Source': source,
                'Condition': condition,
                'Replicate': replicate,
                'Group': group,
                'Total_Counts': self.count_matrix[sample].sum()
            })
        return pd.DataFrame(meta).set_index('Sample')

    def run(self):
        logger.info("=== Step 3: Quality Control ===")
        metadata = self.parse_metadata()
        
        # 1. PCA
        try:
            from sklearn.decomposition import PCA
            import matplotlib.pyplot as plt
            import seaborn as sns
            self.plot_pca(metadata)
        except ImportError:
            logger.warning("Skipping PCA plot: sklearn/matplotlib/seaborn not installed or crashing")
        except Exception as e:
            logger.warning(f"Skipping PCA plot: {e}")
        
        # 2. Correlation Heatmap
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            self.plot_correlation()
        except ImportError:
            logger.warning("Skipping Correlation plot: matplotlib/seaborn not installed or crashing")
        except Exception as e:
            logger.warning(f"Skipping Correlation plot: {e}")
        
        # 3. Summary Stats
        self.save_summary(metadata)
        
    def plot_pca(self, metadata):
        # Log2(CPM+1)
        cpm = self.count_matrix.div(self.count_matrix.sum(), axis=1) * 1e6
        log_cpm = np.log2(cpm + 1)
        
        from sklearn.decomposition import PCA
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        pca = PCA(n_components=2)
        coords = pca.fit_transform(log_cpm.T)
        
        plt.figure(figsize=(10, 8))
        
        # Map groups to colors
        groups = metadata['Group'].unique()
        colors = sns.color_palette("husl", len(groups))
        group_map = dict(zip(groups, colors))
        
        for i, sample in enumerate(log_cpm.columns):
            grp = metadata.loc[sample, 'Group']
            plt.scatter(coords[i,0], coords[i,1], color=group_map[grp], label=grp if grp not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.text(coords[i,0]+0.1, coords[i,1]+0.1, sample, fontsize=8)
            
        plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        plt.title("PCA of Samples")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'pca_plot.png'))
        plt.close()
        logger.info("✓ Saved PCA plot")

    def plot_correlation(self):
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        cpm = self.count_matrix.div(self.count_matrix.sum(), axis=1) * 1e6
        log_cpm = np.log2(cpm + 1)
        corr = log_cpm.corr()
        
        plt.figure(figsize=(12, 10))
        sns.heatmap(corr, annot=False, cmap='RdYlBu_r')
        plt.title("Sample Correlation Matrix")
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'correlation_plot.png'))
        plt.close()
        logger.info("✓ Saved correlation plot")

    def save_summary(self, metadata):
        summary_file = os.path.join(self.output_dir, 'run_summary.txt')
        with open(summary_file, 'w') as f:
            f.write("Adapted Pipeline Summary\n")
            f.write("========================\n\n")
            f.write(f"Date: {datetime.now()}\n")
            f.write(f"Total Samples: {len(self.count_matrix.columns)}\n")
            f.write(f"Total Targets: {len(self.count_matrix)}\n\n")
            
            f.write("Sample Statistics:\n")
            f.write(metadata[['Total_Counts', 'Source', 'Condition']].to_string())
            
        logger.info(f"✓ Saved summary text to: {summary_file}")

# --- Main Execution ---

def main():
    print("Starting Adapted Pipeline...")
    
    # 1. Prep
    prep = ReferencePrep()
    ref_file = prep.run()
    
    # 2. Quantify & QC for R1
    print("\n--- Processing R1 Reads ---")
    quant = BarcodeQuantifier(ref_file)
    matrix_r1 = quant.run(read_type="R1")
    
    if matrix_r1 is not None:
        qc_r1 = QualityControl(matrix_r1, read_type="R1")
        qc_r1.run()
    
    # 3. Quantify & QC for R2
    print("\n--- Processing R2 Reads ---")
    # Re-initialize to clear stats
    quant = BarcodeQuantifier(ref_file) 
    matrix_r2 = quant.run(read_type="R2")
    
    if matrix_r2 is not None:
        qc_r2 = QualityControl(matrix_r2, read_type="R2")
        qc_r2.run()
    
    print("\nPipeline Completed Successfully!")
    print(f"Results in: {os.path.join(OUTPUT_BASE_DIR, 'results')}")
    print(f"Reports in: {os.path.join(OUTPUT_BASE_DIR, 'qc_reports')}")

if __name__ == "__main__":
    main()
