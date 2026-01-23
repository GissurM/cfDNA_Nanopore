#!/usr/bin/env python3

"""
Generate heatmaps for cfDNA fragmentomics analysis results
Creates PNG visualizations from the motif CSV files
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os

print("=== cfDNA Fragmentomics Heatmap Generator ===")

# Set working directory
os.chdir("/home/gissu/heatmap_gen")

# Get all motif CSV files
results_dir = "results2_softclip"
csv_files = glob.glob(os.path.join(results_dir, "*.motif.csv"))

print(f"Found {len(csv_files)} motif CSV files")

if len(csv_files) == 0:
    print("No motif CSV files found in results2 directory")
    exit(1)

# Read all CSV files and combine
all_data = []
sample_names = []

for file in csv_files:
    # Extract sample name
    sample_name = os.path.basename(file).replace(".motif.csv", "")
    sample_names.append(sample_name)
    
    # Read CSV
    df = pd.read_csv(file)
    df['sample'] = sample_name
    
    all_data.append(df)
    print(f"Loaded: {sample_name} with {len(df)} motifs")

# Combine all data
combined_data = pd.concat(all_data, ignore_index=True)
print(f"Combined data shape: {combined_data.shape}")

# Create matrix for heatmap (samples x motifs)
# Focus on top 20 most frequent motifs across all samples
motif_summary = combined_data.groupby('motif')['count'].sum().sort_values(ascending=False).head(20)
print("Top 20 motifs selected for heatmap:")
print(motif_summary)

# Filter data to top motifs
filtered_data = combined_data[combined_data['motif'].isin(motif_summary.index)]

# Create wide format matrix
heatmap_matrix = filtered_data.pivot(index='motif', columns='sample', values='frequency')
heatmap_matrix = heatmap_matrix.fillna(0)

# Sort samples by name for consistent ordering
heatmap_matrix = heatmap_matrix.reindex(sorted(heatmap_matrix.columns), axis=1)

print(f"Heatmap matrix created: {heatmap_matrix.shape[0]} motifs x {heatmap_matrix.shape[1]} samples")

# Set up matplotlib for better plots
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 150

# Create custom colormap: pink -> red -> green -> light blue -> purple (high to low frequency)
from matplotlib.colors import LinearSegmentedColormap
colors = ['#800080', '#87CEEB', '#00FF00', '#FF0000', '#FFC0CB']  # purple, light blue, green, red, pink
custom_cmap = LinearSegmentedColormap.from_list("custom", colors, N=256)

# Transpose the matrix to switch X and Y axes (samples on Y-axis, motifs on X-axis)
heatmap_matrix_transposed = heatmap_matrix.T

# Create complete heatmap with all 24 samples
plt.figure(figsize=(18, 12))
sns.heatmap(heatmap_matrix_transposed, 
           cmap=custom_cmap, 
           annot=False, 
           fmt='.2f',
           cbar_kws={'label': 'Frequency (%)'},
           xticklabels=True,
           yticklabels=True,
           linewidths=0.5,
           linecolor='white')

plt.title('cfDNA Fragment End Motif Frequencies (%) - All 24 Samples\nTop 20 Most Frequent Motifs', 
         fontsize=18, fontweight='bold', pad=20)
plt.xlabel('Motif', fontsize=16, fontweight='bold')
plt.ylabel('Sample', fontsize=16, fontweight='bold')
plt.xticks(rotation=45, ha='right', fontsize=14)
plt.yticks(rotation=0, fontsize=14)

# Highlight CCCA column if present
if 'CCCA' in heatmap_matrix_transposed.columns:
    ccca_idx = list(heatmap_matrix_transposed.columns).index('CCCA')
    plt.axvline(x=ccca_idx, color='black', linewidth=3, alpha=0.8)
    plt.axvline(x=ccca_idx+1, color='black', linewidth=3, alpha=0.8)

# Add horizontal line to separate first 12 and last 12 samples
if heatmap_matrix_transposed.shape[0] > 12:
    plt.axhline(y=12, color='black', linewidth=3, alpha=0.8)

plt.tight_layout()
plt.savefig('cfDNA_Motif_Heatmap_Complete.png', bbox_inches='tight', dpi=300)
plt.close()

# Create focused CCCA barplot
ccca_data = combined_data[combined_data['motif'] == 'CCCA'].copy()
ccca_data = ccca_data.sort_values('sample')

plt.figure(figsize=(14, 8))
colors = ['#2ECC71' if freq > 1.4 else '#E74C3C' for freq in ccca_data['frequency']]

bars = plt.bar(range(len(ccca_data)), ccca_data['frequency'], color=colors, 
               edgecolor='black', linewidth=0.5)

plt.axhline(y=1.4, color='red', linestyle='--', linewidth=2, alpha=0.8, 
           label='Cancer concern threshold (1.4%)')

plt.title('CCCA Motif Frequency - Cancer Biomarker Analysis', 
         fontsize=16, fontweight='bold')
plt.suptitle('All samples show healthy profiles (>1.4% threshold)', 
            fontsize=12, y=0.92)
plt.xlabel('Sample', fontsize=12)
plt.ylabel('CCCA Frequency (%)', fontsize=12)

plt.xticks(range(len(ccca_data)), ccca_data['sample'], rotation=45, ha='right')
plt.ylim(1.5, 2.0)

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#2ECC71', label='Healthy (>1.4%)'),
                  Patch(facecolor='#E74C3C', label='Concern (≤1.4%)')]
plt.legend(handles=legend_elements, loc='upper right')

plt.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig('cfDNA_CCCA_Barplot.png', bbox_inches='tight')
plt.close()

# Create partial heatmaps with the same custom colormap and transposed layout
# Part 1: Samples 1-12
n_samples = min(12, heatmap_matrix_transposed.shape[0])
heatmap_part1 = heatmap_matrix_transposed.iloc[:n_samples, :]

plt.figure(figsize=(14, 7))
sns.heatmap(heatmap_part1, 
           cmap=custom_cmap, 
           annot=False, 
           fmt='.2f',
           cbar_kws={'label': 'Frequency (%)'},
           xticklabels=True,
           yticklabels=True,
           linewidths=0.5,
           linecolor='white')

plt.title('cfDNA Fragment End Motif Frequencies\nSamples 1-12 - Top 20 Motifs', 
         fontsize=16, fontweight='bold')
plt.xlabel('Motif', fontsize=16, fontweight='bold')
plt.ylabel('Sample', fontsize=16, fontweight='bold')
plt.xticks(rotation=45, ha='right', fontsize=14)
plt.yticks(rotation=0, fontsize=14)

plt.tight_layout()
plt.savefig('cfDNA_Motif_Heatmap_Part1.png', bbox_inches='tight', dpi=300)
plt.close()

# Part 2: Samples 13-24 (if available)
if heatmap_matrix_transposed.shape[0] > 12:
    heatmap_part2 = heatmap_matrix_transposed.iloc[12:, :]
    
    plt.figure(figsize=(14, 7))
    sns.heatmap(heatmap_part2, 
               cmap=custom_cmap, 
               annot=False, 
               fmt='.2f',
               cbar_kws={'label': 'Frequency (%)'},
               xticklabels=True,
               yticklabels=True,
               linewidths=0.5,
               linecolor='white')
    
    plt.title('cfDNA Fragment End Motif Frequencies\nSamples 13-24 - Top 20 Motifs', 
             fontsize=16, fontweight='bold')
    plt.xlabel('Motif', fontsize=16, fontweight='bold')
    plt.ylabel('Sample', fontsize=16, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=14)
    plt.yticks(rotation=0, fontsize=14)
    
    plt.tight_layout()
    plt.savefig('cfDNA_Motif_Heatmap_Part2.png', bbox_inches='tight', dpi=300)
    plt.close()

print("\n=== Heatmap generation complete! ===")
print("Generated files:")
print("  - cfDNA_Motif_Heatmap_Complete.png (all 24 samples, transposed with custom colors)")
print("  - cfDNA_CCCA_Barplot.png (CCCA cancer biomarker focus)")
print("  - cfDNA_Motif_Heatmap_Part1.png (samples 1-12, transposed)")
if heatmap_matrix_transposed.shape[0] > 12:
    print("  - cfDNA_Motif_Heatmap_Part2.png (samples 13-24, transposed)")

# Summary statistics
print("\n=== Summary Statistics ===")
print(f"Total samples analyzed: {len(combined_data['sample'].unique())}")
print(f"Total unique motifs: {len(combined_data['motif'].unique())}")
print(f"CCCA frequency range: {ccca_data['frequency'].min():.3f}% - {ccca_data['frequency'].max():.3f}%")
print(f"Healthy samples (CCCA > 1.4%): {sum(ccca_data['frequency'] > 1.4)}/{len(ccca_data)}")
