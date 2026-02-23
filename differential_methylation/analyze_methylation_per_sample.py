#!/usr/bin/env python3
"""
Promoter Methylation Analysis - Per-Sample with Permutation Testing
Improved statistical approach with sample-level analysis and permutation tests
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict
import subprocess
from typing import Dict, List, Tuple
import warnings
import requests
from scipy.stats import chi2_contingency
warnings.filterwarnings('ignore')

# Configuration
GENES_OF_INTEREST = {
    'ACTL6A': None, 'APC': None, 'ATR': None, 'BRCA1': None, 'BRCA2': None,
    'CCND2': None, 'DKK3': None, 'ERCC4': None, 'ERCC5': None, 'ESR1': None,
    'EYA4': None, 'FANCA': None, 'GAPDH': None, 'GSTP1': None, 'HIN1': None,
    'ITIH5': None, 'PARP4': None, 'PIK3CA': None, 'PI3K': None, 'PUM1': None,
    'RAD50': None, 'RARB': None, 'RFC1': None, 'RFC3': None, 'RRM2B': None,
    'RPPH1': None, 'SFRP1': None, 'SFRP2': None, 'SFRP5': None, 'SHLD1': None,
    'SUMO1': None, 'TFPT': None, 'TWIST1': None, 'UBE2W': None, 'UBR5': None, 'WIF1': None
}

MIN_CPGS = 3  # Minimum valid CpG sites (n_valid_cov) per promoter per sample for inclusion (configurable)
MIN_SAMPLES = 6  # Minimum samples per group to perform analysis
PERMUTATIONS = 10000  # Number of permutations for p-value calculation
# Note: MIN_CPGS is based on column 10 (n_valid_cov) from bedMethyl files - total CpG sites in promoter region

class GeneLocationFetcher:
    """Fetch gene locations from Ensembl"""
    
    def __init__(self, assembly='GRCh38'):
        self.assembly = assembly
        self.ensembl_url = "https://rest.ensembl.org"
        self.gene_cache = {}
        
    def get_gene_info(self, gene_name: str) -> Dict:
        """Get gene coordinates from Ensembl"""
        if gene_name in self.gene_cache:
            return self.gene_cache[gene_name]
        
        try:
            url = f"{self.ensembl_url}/lookup/symbol/homo_sapiens/{gene_name}"
            headers = {"Content-Type": "application/json"}
            response = requests.get(url, headers=headers, timeout=10)
            
            if response.ok:
                data = response.json()
                gene_info = {
                    'gene_name': gene_name,
                    'chromosome': data.get('seq_region_name'),
                    'start': data.get('start'),
                    'end': data.get('end'),
                    'strand': data.get('strand'),
                    'ensembl_id': data.get('id')
                }
                # Define promoter as TSS ± 2kb
                if gene_info['strand'] == 1:
                    promoter_start = max(1, gene_info['start'] - 2000)
                    promoter_end = gene_info['start'] + 500
                else:
                    promoter_start = max(1, gene_info['end'] - 500)
                    promoter_end = gene_info['end'] + 2000
                
                gene_info['promoter_start'] = promoter_start
                gene_info['promoter_end'] = promoter_end
                
                self.gene_cache[gene_name] = gene_info
                return gene_info
        except Exception as e:
            print(f"Warning: Could not fetch {gene_name} from Ensembl: {e}")
        
        return None


class PermutationTester:
    """Permutation-based statistical testing"""
    
    @staticmethod
    def extract_sample_methylation(bed_file: str, chrom: str, start: int, end: int) -> Tuple[int, int]:
        """
        Extract methylation metrics for a region.
        Returns: (n_mod, n_valid_cov) where:
          - n_mod: number of methylated CpGs (column 12)
          - n_valid_cov: total valid CpG sites (column 10) - used for coverage weighting
        """
        try:
            cmd = f"grep '^{chrom}\t' {bed_file} | awk '$2 >= {start} && $3 <= {end} {{m+=$12; t+=$10}} END {{print m\"\t\"t}}'"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
            
            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split()
                m = int(parts[0]) if parts and parts[0] else 0
                t = int(parts[1]) if len(parts) > 1 and parts[1] else 0
                return m, t
        except Exception as e:
            pass
        return 0, 0
    
    @staticmethod
    def weighted_mean(values, weights):
        weights = np.array(weights)
        values = np.array(values)
        if np.sum(weights) == 0:
            return np.nan
        return np.sum(values * weights) / np.sum(weights)

    @staticmethod
    def weighted_std(values, weights):
        values = np.array(values)
        weights = np.array(weights)
        if len(values) < 2:
            return np.nan
        w_sum = np.sum(weights)
        if w_sum == 0:
            return np.nan
        mean = np.sum(weights * values) / w_sum

        # Unbiased weighted sample variance
        numerator = np.sum(weights * (values - mean)**2)
        denominator = w_sum - (np.sum(weights**2) / w_sum)
        if denominator <= 0:
            return np.nan
        variance = numerator / denominator
        return np.sqrt(variance)    
    
    @staticmethod
    def permutation_test(ctrl_values: List[float], brca2_values: List[float], 
                        ctrl_weights: List[float], brca2_weights: List[float],
                        n_perms: int = 10000) -> float:
        """
        Two-tailed permutation test with coverage-weighted means
        Samples with higher coverage (more CpGs) have greater influence on the mean.
        
        Args:
            ctrl_values: methylation percentages for control samples
            brca2_values: methylation percentages for BRCA2 samples
            ctrl_weights: n_valid_cov (total CpGs) for each control sample
            brca2_weights: n_valid_cov (total CpGs) for each BRCA2 sample
            n_perms: number of permutations
        
        Returns p-value based on permutation distribution
        """
        # Calculate observed difference in weighted means
        ctrl_wmean = PermutationTester.weighted_mean(ctrl_values, ctrl_weights)
        brca2_wmean = PermutationTester.weighted_mean(brca2_values, brca2_weights)
        observed_diff = abs(brca2_wmean - ctrl_wmean)
        
        # Combine all values and weights
        all_values = np.array(ctrl_values + brca2_values)
        all_weights = np.array(ctrl_weights + brca2_weights)
        group_labels = np.array([0] * len(ctrl_values) + [1] * len(brca2_values))  # 0=ctrl, 1=BRCA2
        n_ctrl = len(ctrl_values)
        
        # Permutation test
        count = 0
        for _ in range(n_perms):
            # Shuffle group labels while keeping values and weights paired
            shuffled_labels = np.random.permutation(group_labels)
            
            # Separate by permuted labels
            perm_ctrl_vals = all_values[shuffled_labels == 0]
            perm_ctrl_wts = all_weights[shuffled_labels == 0]
            perm_brca2_vals = all_values[shuffled_labels == 1]
            perm_brca2_wts = all_weights[shuffled_labels == 1]
            
            # Calculate weighted means
            perm_ctrl_wmean = PermutationTester.weighted_mean(perm_ctrl_vals, perm_ctrl_wts)
            perm_brca2_wmean = PermutationTester.weighted_mean(perm_brca2_vals, perm_brca2_wts)
            perm_diff = abs(perm_brca2_wmean - perm_ctrl_wmean)
            
            if perm_diff >= observed_diff:
                count += 1
        
        p_value = (count + 1) / (n_perms + 1)
        return p_value
    
    @staticmethod
    def benjamini_hochberg(p_values: List[float], alpha: float = 0.05) -> Tuple[List[float], List[bool]]:
        """BH FDR correction"""
        n = len(p_values)
        sorted_idx = np.argsort(p_values)
        sorted_p = np.array(p_values)[sorted_idx]
        
        adjusted_p = np.ones(n)
        for i, idx in enumerate(sorted_idx):
            adjusted_p[idx] = sorted_p[i] * n / (i + 1)
        
        adjusted_p = np.minimum.accumulate(adjusted_p[::-1])[::-1]
        adjusted_p = np.minimum(adjusted_p, 1.0)
        
        return list(adjusted_p), list(adjusted_p <= alpha)


class PerSampleMethylationAnalyzer:
    """Analyze methylation per-sample with permutation testing"""
    
    def __init__(self, bed_dir: str, output_dir: str = './methylation_results',
                 min_cpgs: int = MIN_CPGS, min_samples: int = MIN_SAMPLES, 
                 n_perms: int = PERMUTATIONS):
        self.bed_dir = Path(bed_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        self.min_cpgs = min_cpgs
        self.min_samples = min_samples
        self.n_perms = n_perms
        
        self.ctrl_bed_files = []
        self.brca2_bed_files = []
        self.gene_fetcher = GeneLocationFetcher()
        self.perm_tester = PermutationTester()
        self.analysis_results = defaultdict(dict)
        
    def find_bed_files(self):
        """Find and segregate bed files by group"""
        bed_files = sorted(self.bed_dir.glob('*-barcode*.bed'))
        
        if not bed_files:
            raise FileNotFoundError(f"No bed files found in {self.bed_dir}")
        
        print(f"Found {len(bed_files)} bed files")
        
        for bed_file in bed_files:
            filename = bed_file.stem
            match = re.match(r'(ctrl|BRCA2)-barcode(\d+)', filename, re.IGNORECASE)
            if match:
                group = match.group(1).lower()
                barcode = match.group(2)
                
                if group == 'ctrl':
                    self.ctrl_bed_files.append((barcode, str(bed_file)))
                else:
                    self.brca2_bed_files.append((barcode, str(bed_file)))
        
        print(f"Loaded {len(self.ctrl_bed_files)} ctrl samples and {len(self.brca2_bed_files)} BRCA2 samples")
        
        if len(self.ctrl_bed_files) < self.min_samples or len(self.brca2_bed_files) < self.min_samples:
            print(f"WARNING: Less than {self.min_samples} samples in at least one group!")
    
    def analyze_gene_promoters(self):
        """Per-sample analysis with permutation testing"""
        print("\n" + "="*70)
        print(f"Per-Sample Methylation Analysis (Min CpGs: {self.min_cpgs}, Min Samples: {self.min_samples})")
        print("="*70)
        
        gene_list = list(GENES_OF_INTEREST.keys())
        p_values = []
        genes_analyzed = []
        
        for i, gene_name in enumerate(gene_list, 1):
            print(f"[{i}/{len(gene_list)}] {gene_name}...", end=' ')
            
            gene_info = self.gene_fetcher.get_gene_info(gene_name)
            if not gene_info:
                print("SKIPPED (no coordinates)")
                continue
            
            chrom = f"chr{gene_info['chromosome']}"
            start = gene_info['promoter_start']
            end = gene_info['promoter_end']
            
            # Per-sample methylation percentages with coverage weights
            ctrl_meth_pcts = []
            ctrl_weights = []  # n_valid_cov for each sample
            brca2_meth_pcts = []
            brca2_weights = []  # n_valid_cov for each sample
            
            # Extract for control samples
            for barcode, bed_file in self.ctrl_bed_files:
                meth, total = self.perm_tester.extract_sample_methylation(bed_file, chrom, start, end)
                if total >= self.min_cpgs:  # Filter by minimum valid CpGs in region
                    meth_pct = (meth / total * 100)
                    ctrl_meth_pcts.append(meth_pct)
                    ctrl_weights.append(total)  # Weight by coverage (n_valid_cov)
            
            # Extract for BRCA2 samples
            for barcode, bed_file in self.brca2_bed_files:
                meth, total = self.perm_tester.extract_sample_methylation(bed_file, chrom, start, end)
                if total >= self.min_cpgs:  # Filter by minimum valid CpGs in region
                    meth_pct = (meth / total * 100)
                    brca2_meth_pcts.append(meth_pct)
                    brca2_weights.append(total)  # Weight by coverage (n_valid_cov)
            
            # Check if enough samples passed filter
            if len(ctrl_meth_pcts) < self.min_samples or len(brca2_meth_pcts) < self.min_samples:
                print(f"SKIPPED (insufficient samples: ctrl={len(ctrl_meth_pcts)}, brca2={len(brca2_meth_pcts)})")
                continue
            
            # Calculate statistics (weighted by coverage)
            ctrl_mean = self.perm_tester.weighted_mean(ctrl_meth_pcts, ctrl_weights)
            ctrl_sd = self.perm_tester.weighted_std(ctrl_meth_pcts, ctrl_weights)
            brca2_mean = self.perm_tester.weighted_mean(brca2_meth_pcts, brca2_weights)
            brca2_sd = self.perm_tester.weighted_std(brca2_meth_pcts, brca2_weights)
            
            # Permutation test (with coverage weighting)
            p_value = self.perm_tester.permutation_test(ctrl_meth_pcts, brca2_meth_pcts, 
                                                        ctrl_weights, brca2_weights, self.n_perms)
            p_values.append(p_value)
            genes_analyzed.append(gene_name)
            
            # Determine direction
            diff = brca2_mean - ctrl_mean
            direction = "↑" if diff > 0 else "↓"
            
            self.analysis_results[gene_name] = {
                'chromosome': chrom,
                'promoter_start': start,
                'promoter_end': end,
                'ctrl_mean': ctrl_mean,
                'ctrl_sd': ctrl_sd,
                'ctrl_n_samples': len(ctrl_meth_pcts),
                'brca2_mean': brca2_mean,
                'brca2_sd': brca2_sd,
                'brca2_n_samples': len(brca2_meth_pcts),
                'difference': diff,
                'p_value': p_value,
                'adjusted_p_value': 1.0,
                'is_significant': False,
                'status': 'NOT SIGNIFICANT',
                'ensembl_id': gene_info.get('ensembl_id', 'N/A'),
                'ctrl_values': ctrl_meth_pcts,
                'brca2_values': brca2_meth_pcts
            }
            
            print(f"{direction} p={p_value:.4f} | Ctrl: {ctrl_mean:.1f}±{ctrl_sd:.1f}% | BRCA2: {brca2_mean:.1f}±{brca2_sd:.1f}%")
        
        # Apply FDR correction
        if p_values:
            print("\n" + "="*70)
            print(f"Benjamini-Hochberg FDR Correction (alpha=0.05)")
            print("="*70)
            
            adjusted_p, is_sig = self.perm_tester.benjamini_hochberg(p_values)
            
            for gene_name, adj_p, sig in zip(genes_analyzed, adjusted_p, is_sig):
                data = self.analysis_results[gene_name]
                data['adjusted_p_value'] = adj_p
                data['is_significant'] = sig
                
                if sig:
                    if data['difference'] > 0:
                        data['status'] = "HYPERMETHYLATED *"
                    else:
                        data['status'] = "HYPOMETHYLATED *"
                else:
                    data['status'] = "NOT SIGNIFICANT"
                
                if sig:
                    print(f"{gene_name}: p={data['p_value']:.6f}, FDR-adj={adj_p:.6f} - {data['status']}")
    
    def generate_text_report(self):
        """Generate detailed text report"""
        report_path = self.output_dir / 'methylation_analysis_report.txt'
        
        with open(report_path, 'w') as f:
            f.write("="*80 + "\n")
            f.write("PROMOTER METHYLATION ANALYSIS - PER-SAMPLE WITH PERMUTATION TESTING\n")
            f.write("="*80 + "\n\n")
            
            f.write("METHODOLOGY:\n")
            f.write("-"*80 + "\n")
            f.write(f"- Sample analysis: Per-sample methylation % calculated individually\n")
            f.write(f"- Minimum valid CpG sites per sample: {self.min_cpgs} (based on n_valid_cov from bedMethyl files)\n")
            f.write(f"- Coverage weighting: Mean calculations weighted by n_valid_cov (higher coverage = higher weight)\n")
            f.write(f"- Minimum samples per group: {self.min_samples}\n")
            f.write(f"- Statistical test: Two-tailed permutation test with coverage weighting ({self.n_perms:,} permutations)\n")
            f.write(f"- Multiple testing correction: Benjamini-Hochberg FDR (alpha=0.05)\n")
            f.write(f"- Results marked with * have FDR-adjusted p-value < 0.05\n\n")
            
            f.write(f"SAMPLES:\n")
            f.write("-"*80 + "\n")
            f.write(f"Control (n={len(self.ctrl_bed_files)}): {', '.join([b for b,_ in self.ctrl_bed_files])}\n")
            f.write(f"BRCA2 (n={len(self.brca2_bed_files)}): {', '.join([b for b,_ in self.brca2_bed_files])}\n\n")
            
            # Summary
            status_counts = defaultdict(int)
            for gene_data in self.analysis_results.values():
                status_counts[gene_data['status']] += 1
            
            f.write("SUMMARY\n")
            f.write("-"*80 + "\n")
            f.write(f"Genes analyzed: {len(self.analysis_results)}\n")
            for status in sorted(status_counts.keys()):
                f.write(f"  {status}: {status_counts[status]}\n")
            f.write("\n")
            
            # Detailed results
            f.write("DETAILED RESULTS\n")
            f.write("-"*80 + "\n\n")
            
            sorted_genes = sorted(self.analysis_results.keys(), 
                                 key=lambda g: self.analysis_results[g]['adjusted_p_value'])
            
            for gene_name in sorted_genes:
                data = self.analysis_results[gene_name]
                f.write(f"Gene: {gene_name}\n")
                f.write(f"  Ensembl ID: {data['ensembl_id']}\n")
                f.write(f"  Location: {data['chromosome']}:{data['promoter_start']}-{data['promoter_end']}\n")
                f.write(f"  Control: {data['ctrl_mean']:.2f}% ± {data['ctrl_sd']:.2f}% (n={data['ctrl_n_samples']} samples, coverage-weighted mean)\n")
                f.write(f"  BRCA2:   {data['brca2_mean']:.2f}% ± {data['brca2_sd']:.2f}% (n={data['brca2_n_samples']} samples, coverage-weighted mean)\n")
                f.write(f"  Difference: {data['difference']:.2f}%\n")
                f.write(f"  P-value (permutation): {data['p_value']:.6f}\n")
                f.write(f"  P-value (FDR-adjusted): {data['adjusted_p_value']:.6f}\n")
                f.write(f"  Status: {data['status']}\n")
                f.write("\n")
        
        print(f"\nReport saved to {report_path}")
        return report_path
    
    def generate_histogram(self):
        """Generate histogram for significant genes"""
        sig_genes = sorted([g for g in self.analysis_results.keys() 
                           if self.analysis_results[g]['is_significant']])
        
        if not sig_genes:
            print("No statistically significant genes found")
            return None
        
        fig, ax = plt.subplots(figsize=(14, 8))
        
        x = np.arange(len(sig_genes))
        width = 0.35
        
        ctrl_means = [self.analysis_results[g]['ctrl_mean'] for g in sig_genes]
        ctrl_sds = [self.analysis_results[g]['ctrl_sd'] for g in sig_genes]
        brca2_means = [self.analysis_results[g]['brca2_mean'] for g in sig_genes]
        brca2_sds = [self.analysis_results[g]['brca2_sd'] for g in sig_genes]
        
        bars1 = ax.bar(x - width/2, ctrl_means, width, yerr=ctrl_sds, label='Control', 
                      alpha=0.8, color='#2E86AB', capsize=5)
        bars2 = ax.bar(x + width/2, brca2_means, width, yerr=brca2_sds, label='BRCA2', 
                      alpha=0.8, color='#A23B72', capsize=5)
        
        ax.set_xlabel('Gene', fontsize=12, fontweight='bold')
        ax.set_ylabel('Methylation Percentage (%)', fontsize=12, fontweight='bold')
        ax.set_title(f'Significant Promoter Methylation Differences (FDR < 0.05)\nMean ± SD (n={len(self.ctrl_bed_files)} ctrl, {len(self.brca2_bed_files)} BRCA2)', 
                     fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(sig_genes, rotation=45, ha='right')
        ax.legend(fontsize=11)
        ax.grid(axis='y', alpha=0.3)
        ax.set_ylim(0, 100)
        
        plt.tight_layout()
        histogram_path = self.output_dir / 'methylation_histogram_significant.png'
        plt.savefig(histogram_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Histogram saved to {histogram_path}")
        return histogram_path
    
    def generate_summary_table(self):
        """Generate summary table"""
        data = []
        for gene_name in sorted(self.analysis_results.keys()):
            d = self.analysis_results[gene_name]
            data.append({
                'Gene': gene_name,
                'Chromosome': d['chromosome'],
                'Ctrl Mean (%)': f"{d['ctrl_mean']:.2f}",
                'Ctrl SD': f"{d['ctrl_sd']:.2f}",
                'Ctrl N': d['ctrl_n_samples'],
                'BRCA2 Mean (%)': f"{d['brca2_mean']:.2f}",
                'BRCA2 SD': f"{d['brca2_sd']:.2f}",
                'BRCA2 N': d['brca2_n_samples'],
                'Difference (%)': f"{d['difference']:.2f}",
                'P-value': f"{d['p_value']:.6f}",
                'FDR P-value': f"{d['adjusted_p_value']:.6f}",
                'Significant': 'Yes' if d['is_significant'] else 'No',
                'Status': d['status']
            })
        
        df = pd.DataFrame(data)
        
        csv_path = self.output_dir / 'methylation_summary_table.csv'
        df.to_csv(csv_path, index=False)
        
        html_path = self.output_dir / 'methylation_summary_table.html'
        html_content = df.to_html(index=False, escape=False)
        html_content = html_content.replace('HYPERMETHYLATED *', '<span style="color: red; font-weight: bold;">HYPERMETHYLATED *</span>')
        html_content = html_content.replace('HYPOMETHYLATED *', '<span style="color: blue; font-weight: bold;">HYPOMETHYLATED *</span>')
        html_content = html_content.replace('NOT SIGNIFICANT', '<span style="color: gray;">NOT SIGNIFICANT</span>')
        
        with open(html_path, 'w') as f:
            f.write(f"<html><head><style>table {{border-collapse: collapse}} th, td {{border: 1px solid #999; padding: 6px}}</style></head>" \
                    f"<body><h2>Per-Sample Methylation Analysis with Permutation Testing</h2>" \
                    f"<p>Min CpGs: {self.min_cpgs} | Min Samples: {self.min_samples} | Permutations: {self.n_perms:,}</p>" \
                    f"{html_content}</body></html>")
        
        print(f"Summary table saved to {csv_path}")
        return csv_path
    
    def run_analysis(self):
        """Run complete analysis"""
        print("Starting per-sample analysis with permutation testing...")
        np.random.seed(42)
        self.find_bed_files()
        self.analyze_gene_promoters()
        self.generate_text_report()
        self.generate_histogram()
        self.generate_summary_table()
        print("\n" + "="*70)
        print("Analysis Complete!")
        print(f"Results saved to {self.output_dir}")
        print("="*70)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Per-Sample Methylation Analysis with Permutation Testing')
    parser.add_argument('bed_dir', help='Directory containing bedMethyl files')
    parser.add_argument('--output-dir', default='./methylation_results', help='Output directory')
    parser.add_argument('--min-cpgs', type=int, default=MIN_CPGS, 
                       help=f'Minimum CpGs per sample to include (default: {MIN_CPGS})')
    parser.add_argument('--min-samples', type=int, default=MIN_SAMPLES,
                       help=f'Minimum samples per group to analyze (default: {MIN_SAMPLES})')
    parser.add_argument('--permutations', type=int, default=PERMUTATIONS,
                       help=f'Number of permutations (default: {PERMUTATIONS})')
    
    args = parser.parse_args()
    
    if not Path(args.bed_dir).exists():
        print(f"Error: Directory {args.bed_dir} does not exist")
        return
    
    analyzer = PerSampleMethylationAnalyzer(args.bed_dir, args.output_dir, 
                                           args.min_cpgs, args.min_samples, 
                                           args.permutations)
    analyzer.run_analysis()


if __name__ == '__main__':
    main()
