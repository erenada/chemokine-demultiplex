#!/usr/bin/env python3

"""
Cross-Check R1 vs R2 Count Matrices
Author: Eren Ada, PhD (Adapted by AI Assistant)
Date: 02/10/2026

This script compares the count matrices from R1 and R2 reads to validate the pipeline results.
It calculates the correlation between the two matrices and generates a scatter plot.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# Configuration
BASE_DIR = "/Users/eren/Desktop/HMS/EricLoo/project_miseq"
R1_FILE = os.path.join(BASE_DIR, "results", "R1", "raw_count_matrix.csv")
R2_FILE = os.path.join(BASE_DIR, "results", "R2", "raw_count_matrix.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "qc_reports", "comparison")

def main():
    print("Starting R1 vs R2 Cross-Check...")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1. Load Matrices
    if not os.path.exists(R1_FILE) or not os.path.exists(R2_FILE):
        print("Error: Input files not found.")
        return

    df_r1 = pd.read_csv(R1_FILE, index_col=0)
    df_r2 = pd.read_csv(R2_FILE, index_col=0)

    print(f"Loaded R1 Matrix: {df_r1.shape}")
    print(f"Loaded R2 Matrix: {df_r2.shape}")

    # 2. Ensure alignment
    # Keep only shared indices (targets) and columns (samples)
    shared_index = df_r1.index.intersection(df_r2.index)
    shared_cols = df_r1.columns.intersection(df_r2.columns)

    df_r1 = df_r1.loc[shared_index, shared_cols]
    df_r2 = df_r2.loc[shared_index, shared_cols]

    # 3. Flatten for correlation
    r1_values = df_r1.values.flatten()
    r2_values = df_r2.values.flatten()

    # 4. Calculate Statistics
    correlation, p_value = pearsonr(r1_values, r2_values)
    total_r1 = r1_values.sum()
    total_r2 = r2_values.sum()
    diff_percent = abs(total_r1 - total_r2) / ((total_r1 + total_r2) / 2) * 100

    print("\n--- Validation Statistics ---")
    print(f"Pearson Correlation (R): {correlation:.4f}")
    print(f"Total Counts R1: {total_r1:,}")
    print(f"Total Counts R2: {total_r2:,}")
    print(f"Difference: {diff_percent:.2f}%")

    # 5. Generate Plot
    plt.figure(figsize=(10, 8))
    plt.scatter(r1_values, r2_values, alpha=0.5, s=10)
    
    # Add diagonal line
    max_val = max(r1_values.max(), r2_values.max())
    plt.plot([0, max_val], [0, max_val], 'r--', alpha=0.7, label='Perfect Match')
    
    plt.xlabel('R1 Counts')
    plt.ylabel('R2 Counts')
    plt.title(f'R1 vs R2 Count Comparison\nR = {correlation:.4f}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Use log scale if range is huge
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1, max_val * 1.5)
    plt.ylim(1, max_val * 1.5)

    plot_file = os.path.join(OUTPUT_DIR, "r1_r2_comparison.png")
    plt.savefig(plot_file)
    print(f"\nSaved comparison plot to: {plot_file}")

    # 6. Save Difference Matrix
    # (R1 - R2) to see specific discrepancies
    diff_matrix = df_r1 - df_r2
    diff_file = os.path.join(OUTPUT_DIR, "r1_r2_difference_matrix.csv")
    diff_matrix.to_csv(diff_file)
    print(f"Saved difference matrix to: {diff_file}")

if __name__ == "__main__":
    main()
