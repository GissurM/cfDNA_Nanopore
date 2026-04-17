"""
IchorCNA Analysis Pipeline for Nanopore BAM Files (WSL Version)
Based on published study parameters for low-coverage Nanopore cfDNA analysis

This pipeline implements the exact filtering and analysis steps described in:
- Samtools filtering: unmapped reads, secondary/supplementary reads, MAPQ < 20, reads > 700bp
- IchorCNA analysis with specified parameters
- Automated processing of multiple BAM files

Requirements:
- samtools (version 1.9 or later)
- R with ichorCNA package installed
- HMMcopy suite for read counting
"""

import os
import subprocess
import pandas as pd
import glob
import sys
import shutil
import shlex
import html
from pathlib import Path
import logging
import argparse
from datetime import datetime

class IchorCNAPipelineWSL:
    def __init__(self, input_dir, output_dir, reference_genome=None, normal_restarts='c(0.5,0.8)'):
        """
        Initialize IchorCNA pipeline for WSL environment.
        
        Args:
            input_dir: Directory containing BAM files (can be Windows or WSL path)
            output_dir: Directory for output files (WSL path)
            reference_genome: Path to reference genome (optional)
            normal_restarts: ichorCNA normal restart vector, e.g. 'c(0.5,0.8)'
        """
        # Convert Windows paths to WSL paths if needed
        self.input_dir = self.convert_to_wsl_path(input_dir)
        self.output_dir = Path(output_dir)
        self.reference_genome = reference_genome
        
        # Create output directories
        self.filtered_dir = self.output_dir / "filtered_bams"
        self.readcount_dir = self.output_dir / "readcounts"
        self.ichorcna_dir = self.output_dir / "ichorCNA_results"
        self.logs_dir = self.output_dir / "logs"
        
        # Create directories using WSL
        for dir_path in [self.filtered_dir, self.readcount_dir, self.ichorcna_dir, self.logs_dir]:
            self.create_wsl_directory(str(dir_path))
        
        # Setup logging
        log_file = self.logs_dir / f'ichorcna_pipeline_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(str(log_file)),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
        # Prefer less-sensitive normal restarts to reduce minor false-positive
        # TF calls in low-coverage cfDNA samples.
        if not isinstance(normal_restarts, str) or not normal_restarts.strip().startswith('c('):
            raise ValueError("normal_restarts must be an R vector string like 'c(0.5,0.8)'.")

        self.filter_params = {
            'min_mapq': 20,
            'max_read_length': 700,
        }
        self.readcount_params = {
            'window': 1000000,
            'quality': 20,
        }

        self.ichorcna_params = {
            'ploidy': 'c(2)',
            'normal': normal_restarts,
            'maxCN': '7',
            'includeHOMD': 'FALSE',
            'estimateNormal': 'TRUE',
            'estimatePloidy': 'TRUE',
            'estimateScPrevalence': 'FALSE',
            'altFracThreshold': '0.05',
            'rmCentromereFlankLength': '1000000',
            'minMapScore': '0.95',
            'minSegmentBins': '100',
            'maxFracCNASubclone': '0.50',
            'maxFracGenomeSubclone': '0.30',
            'txnE': '0.99999999',
            'txnStrength': '100000000',
            'genomeBuild': 'hg38',
        }
    
    def convert_to_wsl_path(self, windows_path):
        """Convert Windows path to WSL path."""
        if isinstance(windows_path, str) and windows_path.startswith('C:'):
            # Convert C:\path to /mnt/c/path
            wsl_path = windows_path.replace('C:', '/mnt/c').replace('\\', '/')
            return Path(wsl_path)
        return Path(windows_path)
    
    def create_wsl_directory(self, path):
        """Create directory natively in WSL."""
        cmd = f'mkdir -p {path}'
        subprocess.run(cmd, shell=True, check=True)
    
    def run_wsl_command(self, command, capture_output=True):
        """Run command natively in WSL, handling encoding errors."""
        if capture_output:
            return subprocess.run(command, shell=True, capture_output=True, text=True, encoding='utf-8', errors='replace')
        else:
            return subprocess.run(command, shell=True)
    
    def check_dependencies(self):
        """Check if required tools are available in WSL."""
        self.logger.info("Checking dependencies in WSL...")
        
        dependencies = [
            ('samtools --version', 'samtools'),
            ('R --version', 'R')
        ]
        
        missing = []
        
        for cmd, name in dependencies:
            result = self.run_wsl_command(cmd)
            if result.returncode == 0:
                self.logger.info(f"✓ {name} is available in WSL")
            else:
                missing.append(name)
                self.logger.error(f"✗ {name} not found in WSL")
        
        if missing:
            self.logger.error(f"Missing dependencies: {missing}")
            return False
        
        # Check R packages
        r_check_cmd = '''
        if (!require("ichorCNA", quietly = TRUE)) {
            cat("ichorCNA package not found\\n")
            quit(status = 1)
        }
        if (!require("HMMcopy", quietly = TRUE)) {
            cat("HMMcopy package not found\\n") 
            quit(status = 1)
        }
        cat("R packages OK\\n")
        '''
        
        result = self.run_wsl_command(f'R --slave --no-restore --no-save -e \'{r_check_cmd}\'')
        if result.returncode == 0:
            self.logger.info("✓ Required R packages are available")
            return True
        else:
            self.logger.warning("Some R packages may be missing. Will attempt to install them.")
            return self.install_r_packages()
    
    def install_r_packages(self):
        """Install required R packages."""
        self.logger.info("Installing required R packages...")
        
        install_script = '''
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos="https://cran.r-project.org")
        
        BiocManager::install(c("HMMcopy"))
        
        if (!requireNamespace("devtools", quietly = TRUE))
            install.packages("devtools", repos="https://cran.r-project.org")
        
        devtools::install_github("broadinstitute/ichorCNA")
        
        cat("Installation completed\\n")
        '''
        
        result = self.run_wsl_command(f'R --slave --no-restore --no-save -e \'{install_script}\'')
        
        if result.returncode == 0:
            self.logger.info("✓ R packages installed successfully")
            return True
        else:
            self.logger.error(f"Failed to install R packages: {result.stderr}")
            return False
    
    def find_bam_files(self):
        """Find all BAM files in input directory."""
        # Use WSL to find BAM files
        find_cmd = f'find {self.input_dir} -name "*.bam"'
        result = self.run_wsl_command(find_cmd)
        
        if result.returncode == 0:
            bam_files = [Path(line.strip()) for line in result.stdout.strip().split('\n') if line.strip()]
            self.logger.info(f"Found {len(bam_files)} BAM files in {self.input_dir}")
            return bam_files
        else:
            self.logger.error("Failed to find BAM files")
            return []
    
    def filter_nanopore_bam(self, input_bam):
        """
        Filter BAM file according to study parameters using WSL samtools.
        """
        sample_name = input_bam.stem
        filtered_bam = self.filtered_dir / f"{sample_name}_filtered.bam"
        min_mapq = self.filter_params['min_mapq']
        max_read_length = self.filter_params['max_read_length']
        
        self.logger.info(f"Filtering {input_bam.name}...")
        self.logger.info(f"Using filter thresholds: MAPQ>={min_mapq}, read_length<={max_read_length}")
        
        # Combined samtools filter command for WSL
        filter_cmd = f'''
        samtools view -h -b -F 0x4 -F 0x100 -F 0x800 -F 0x400 -q {min_mapq} {input_bam} | \
        samtools view -h - | \\
        awk 'BEGIN{{OFS="\\t"}} /^@/ || length($10) <= {max_read_length}' | \
        samtools view -b - > {filtered_bam} && \\
        samtools index {filtered_bam}
        '''
        
        result = self.run_wsl_command(filter_cmd)
        
        if result.returncode == 0:
            self.logger.info(f"✓ Filtered BAM saved: {filtered_bam}")
            return filtered_bam
        else:
            self.logger.error(f"Error filtering {input_bam.name}: {result.stderr}")
            return None
    
    def generate_read_counts(self, filtered_bam):
        """Generate read counts using HMMcopy readCounter."""
        sample_name = filtered_bam.stem.replace('_filtered', '')
        readcount_file = self.readcount_dir / f"{sample_name}.wig"
        window_size = self.readcount_params['window']
        readcount_quality = self.readcount_params['quality']
        
        self.logger.info(f"Generating read counts for {sample_name}...")
        self.logger.info(f"Using readCounter settings: window={window_size}, quality>={readcount_quality}")
        
        # Check if readCounter is available, if not use alternative method
        readcounter_check = self.run_wsl_command('which readCounter')
        
        if readcounter_check.returncode == 0:
            # Use HMMcopy readCounter
            readcounter_cmd = f'''
            readCounter --window {window_size} --quality {readcount_quality} \
            --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \\
            {filtered_bam} > {readcount_file}
            '''
        else:
            # Alternative: create read counts using samtools and custom script
            self.logger.warning("readCounter not found, using alternative method")
            readcounter_cmd = f'''
            samtools view -b {filtered_bam} | \\
            samtools bedcov <(awk 'BEGIN{{for(i=1;i<=22;i++) print "chr"i"\\t0\\t1000000"; print "chrX\\t0\\t1000000"; print "chrY\\t0\\t1000000"}}') - > {readcount_file}
            '''
        
        result = self.run_wsl_command(readcounter_cmd)
        
        if result.returncode == 0:
            self.logger.info(f"✓ Read counts generated: {readcount_file}")
            return readcount_file
        else:
            self.logger.error(f"Error generating read counts: {result.stderr}")
            return None
    
    def run_ichorcna(self, readcount_file):
        """Run IchorCNA analysis with study parameters."""
        sample_name = readcount_file.stem
        self.logger.info(f"Running IchorCNA for {sample_name}...")

        # Clear stale sample outputs so old n0.99/n0.95 files are not mixed into new runs.
        sample_out_root = self.ichorcna_dir / sample_name
        if sample_out_root.exists():
            shutil.rmtree(sample_out_root)
        sample_out_root.mkdir(parents=True, exist_ok=True)

        out_dir = f"{sample_out_root}/"
        gc_wig = "/home/gissu/ichorCNA/extdata/gc_hg38_1000kb.wig"
        map_wig = "/home/gissu/ichorCNA/extdata/map_hg38_1000kb.wig"
        self.logger.info(f"Using ichorCNA normal restart(s): {self.ichorcna_params['normal']}")

        # Use the downloaded runichorCNA.R script with GC and map wig files
        cmd = (
            f'Rscript /home/gissu/ichorCNA/runichorCNA.R '
            f'--WIG {readcount_file} '
            f'--id {sample_name} '
            f'--ploidy "{self.ichorcna_params["ploidy"]}" '
            f'--normal "{self.ichorcna_params["normal"]}" '
            f'--maxCN {self.ichorcna_params["maxCN"]} '
            f'--includeHOMD {self.ichorcna_params["includeHOMD"]} '
            f'--estimateNormal {self.ichorcna_params["estimateNormal"]} '
            f'--estimatePloidy {self.ichorcna_params["estimatePloidy"]} '
            f'--estimateScPrevalence {self.ichorcna_params["estimateScPrevalence"]} '
            f'--altFracThreshold {self.ichorcna_params["altFracThreshold"]} '
            f'--rmCentromereFlankLength {self.ichorcna_params["rmCentromereFlankLength"]} '
            f'--minMapScore {self.ichorcna_params["minMapScore"]} '
            f'--minSegmentBins {self.ichorcna_params["minSegmentBins"]} '
            f'--maxFracCNASubclone {self.ichorcna_params["maxFracCNASubclone"]} '
            f'--maxFracGenomeSubclone {self.ichorcna_params["maxFracGenomeSubclone"]} '
            f'--txnE {self.ichorcna_params["txnE"]} '
            f'--txnStrength {self.ichorcna_params["txnStrength"]} '
            f'--genomeBuild {self.ichorcna_params["genomeBuild"]} '
            f'--gcWig {gc_wig} '
            f'--mapWig {map_wig} '
            f'--outDir {out_dir}'
        )
        result = self.run_wsl_command(cmd)
        if result.returncode == 0:
            self.logger.info(f"✓ IchorCNA completed for {sample_name}")
            return True
        else:
            self.logger.error(f"Error running IchorCNA: STDERR: {result.stderr}\nSTDOUT: {result.stdout}")
            return False

    def apply_cna_instability_filter(self, sample_name):
        """
        Post-processing stopgap: if fewer than 15% of the genome is covered by
        copy-number alterations the tumor-fraction estimate is considered unstable
        and is set to 0.

        Reference: used in published ichorCNA analyses of similar low-pass cfDNA
        data — "If the percentage of genome covered by CN alterations was less than
        15%, then the tumor fraction was determined to be unstable and set to 0."

        Returns:
            (cna_fraction, adjusted_tf) where adjusted_tf is 0 if the threshold
            was triggered, or None if the original TF is retained / filter skipped.
        """
        CNA_THRESHOLD = 0.15
        sample_dir = self.ichorcna_dir / sample_name

        # ichorCNA writes *.cna.seg; fall back to any *.seg if not found
        seg_files = list(sample_dir.glob("*.cna.seg"))
        if not seg_files:
            seg_files = list(sample_dir.glob("*.seg"))
        if not seg_files:
            self.logger.warning(
                f"{sample_name}: No seg file found — skipping CNA instability filter"
            )
            return None, None

        seg_file = seg_files[0]
        try:
            seg_df = pd.read_csv(str(seg_file), sep='\t')
        except Exception as exc:
            self.logger.error(f"{sample_name}: Could not read {seg_file.name}: {exc}")
            return None, None

        # Map column names case-insensitively, including sample-prefixed columns
        # like "<sample>.event" and "<sample>.copy.number" from .cna.seg outputs.
        lower_to_orig = {c.lower(): c for c in seg_df.columns}

        def first_existing(names):
            for name in names:
                if name in lower_to_orig:
                    return lower_to_orig[name]
            return None

        def first_suffix(suffixes):
            for col in seg_df.columns:
                low = col.lower()
                if any(low.endswith(sfx) for sfx in suffixes):
                    return col
            return None

        start_col = first_existing(['loc.start', 'start', 'chromstart']) or first_suffix(['.start', '_start'])
        end_col = first_existing(['loc.end', 'end', 'chromend']) or first_suffix(['.end', '_end'])

        if not start_col or not end_col:
            self.logger.warning(
                f"{sample_name}: Cannot identify genomic coordinate columns in "
                f"{seg_file.name} — skipping CNA instability filter"
            )
            return None, None

        seg_df['_len'] = (
            pd.to_numeric(seg_df[end_col], errors='coerce')
            - pd.to_numeric(seg_df[start_col], errors='coerce')
            + 1
        ).clip(lower=0)
        seg_df = seg_df[seg_df['_len'] > 0].copy()
        total_length = seg_df['_len'].sum()
        if total_length == 0:
            self.logger.warning(f"{sample_name}: Total segment length is 0 — skipping filter")
            return None, None

        call_col = (
            first_existing(['call', 'event', 'corrected_call'])
            or first_suffix(['.event', '.call', '.corrected_call'])
        )
        cn_col = (
            first_existing(['copy.number', 'copynumber', 'copy_number', 'corrected_copy_number'])
            or first_suffix(['.copy.number', '.copynumber', '.copy_number', '.corrected_copy_number'])
        )

        if call_col:
            call_values = seg_df[call_col].astype(str).str.upper()
            altered_len = seg_df.loc[
                call_values != 'NEUT', '_len'
            ].sum()
        elif cn_col:
            cn_values = pd.to_numeric(seg_df[cn_col], errors='coerce')
            altered_len = seg_df.loc[cn_values != 2, '_len'].sum()
        else:
            self.logger.warning(
                f"{sample_name}: No call/copy-number column recognized in {seg_file.name} "
                f"(columns: {', '.join(seg_df.columns)}) — skipping CNA instability filter"
            )
            return None, None

        cna_fraction = altered_len / total_length
        self.logger.info(
            f"{sample_name}: {cna_fraction * 100:.1f}% of genome covered by CN alterations"
        )

        if cna_fraction < CNA_THRESHOLD:
            self.logger.warning(
                f"{sample_name}: CNA fraction ({cna_fraction * 100:.1f}%) is below the "
                f"{CNA_THRESHOLD * 100:.0f}% threshold — tumor fraction is unstable, set to 0"
            )
            self._set_tumor_fraction_zero(sample_name)
            return cna_fraction, 0

        self.logger.info(
            f"{sample_name}: CNA fraction ({cna_fraction * 100:.1f}%) passes threshold — TF retained"
        )
        return cna_fraction, None

    def _set_tumor_fraction_zero(self, sample_name):
        """Rewrite the Tumor Fraction line in the ichorCNA params.txt to 0."""
        sample_dir = self.ichorcna_dir / sample_name
        params_files = list(sample_dir.glob("*.params.txt"))
        if not params_files:
            self.logger.warning(
                f"{sample_name}: No params.txt found — cannot override tumor fraction in file"
            )
            return
        params_file = params_files[0]
        try:
            lines = params_file.read_text().splitlines()
            updated = []
            for line in lines:
                if line.startswith("Tumor Fraction:"):
                    updated.append(
                        "Tumor Fraction: 0  "
                        "# overridden: <15% of genome covered by CN alterations (unstable estimate)"
                    )
                else:
                    updated.append(line)
            params_file.write_text('\n'.join(updated) + '\n')
            self.logger.info(f"  Updated {params_file.name}: Tumor Fraction set to 0")
        except Exception as exc:
            self.logger.error(f"Could not update params.txt for {sample_name}: {exc}")
    def _regenerate_plot_with_zero_tf(self, sample_name, cna_fraction=None, force_tf_zero=False):
        """
        Re-render a genome-wide PDF with explicit title details.

        Uses the same best-solution selection logic as runichorCNA.R so the
        first-line title matches the standard genome-wide plot, then appends a
        second line with:
          Tumor fraction: X, Ploidy: X, CNA fraction: X

        If force_tf_zero=True, displayed tumor fraction is overridden to 0.
        Saves a new plot named
        <sample>_genomeWide_TF0_override.pdf alongside the other solution PDFs.
        """
        sample_dir = self.ichorcna_dir / sample_name
        rdata_file = sample_dir / f"{sample_name}.RData"
        if not rdata_file.exists():
            self.logger.warning(
                f"{sample_name}: No .RData file found — cannot regenerate plot with TF=0"
            )
            return

        plot_dir = sample_dir / sample_name
        plot_dir.mkdir(parents=True, exist_ok=True)
        out_plot = plot_dir / f"{sample_name}_genomeWide_TF0_override"
        cna_fraction_label = "N/A" if cna_fraction is None else f"{cna_fraction:.4f}"
        force_tf_zero_r = "TRUE" if force_tf_zero else "FALSE"

        r_script = f"""\
library(ichorCNA)
library(HMMcopy)
library(GenomicRanges)
library(GenomeInfoDb)
options(bitmapType='cairo')

# Restore the complete workspace saved by runichorCNA.R
load("{rdata_file}")

mainTitle <- NULL

# Recreate runichorCNA.R best-solution selection for title/metadata consistency
if (exists("results") && exists("loglik")) {{
    ll <- loglik
    if (is.data.frame(ll) && "init" %in% colnames(ll) && "loglik" %in% colnames(ll)) {{
        ll <- ll[!is.na(ll$init), , drop = FALSE]
        if (nrow(ll) > 0) {{
            ord <- order(as.numeric(ll$loglik), decreasing = TRUE)
            best_idx <- ord[1]
            if (length(results) >= best_idx) {{
                hmmResults.cor <- results[[best_idx]]
            }}
            if (exists("mainName") && length(mainName) >= best_idx) {{
                mainTitle <- as.character(mainName[best_idx])
            }}
        }}
    }}
}}

if (!is.null(mainTitle) && !is.na(mainTitle) && mainTitle != "") {{
    mainTitle <- sub(", n: [^,]+", "", mainTitle)
}}

iter <- hmmResults.cor$results$iter
if (is.null(mainTitle) || is.na(mainTitle) || mainTitle == "") {{
    p_est <- as.numeric(hmmResults.cor$results$phi[1, iter])
    ll_val <- as.numeric(hmmResults.cor$results$loglik[iter])
    mainTitle <- sprintf(
        "%s, p: %.4g, log likelihood: %.4g",
        "{sample_name}", p_est, ll_val
    )
}}

if ({force_tf_zero_r}) {{
    hmmResults.cor$results$n[1, iter] <- 1.0
}}

cnaFractionLabel <- "{cna_fraction_label}"
cnaLine <- sprintf("CNA fraction: %s", cnaFractionLabel)

pdf(file = paste0("{out_plot}", ".pdf"), width = 14, height = 6)

plotGWSolution(
  hmmResults.cor, s = 1,
  outPlotFile = "{out_plot}",
  plotFileType = "pdf",
  logR.column = "logR",
  call.column = "Corrected_Call",
  plotYLim = plotYLim,
  estimateScPrevalence = estimateScPrevalence,
  seqinfo = seqinfo,
    turnDevOn = FALSE,
    turnDevOff = FALSE,
    main = mainTitle
)
mtext(line = -2.2, cnaLine, cex = 1.2)
dev.off()
cat("Plot saved\\n")
"""
        tmp_script = sample_dir / f"_replot_{sample_name}.R"
        try:
            tmp_script.write_text(r_script)
            result = self.run_wsl_command(f"Rscript {tmp_script}")
            if result.returncode == 0:
                self.logger.info(
                    f"  Regenerated TF override plot with title details: {out_plot}.pdf"
                )
            else:
                self.logger.warning(
                    f"  Plot regeneration failed for {sample_name}: "
                    f"{result.stderr[:300]}"
                )
        except Exception as exc:
            self.logger.error(f"Could not regenerate plot for {sample_name}: {exc}")
        finally:
            if tmp_script.exists():
                tmp_script.unlink()

    def process_sample(self, input_bam):
        """Process a single BAM file through the complete pipeline."""
        self.logger.info(f"Processing {input_bam.name}...")
        
        # Step 1: Filter BAM
        filtered_bam = self.filter_nanopore_bam(input_bam)
        if not filtered_bam:
            return None
        
        # Step 2: Generate read counts
        readcount_file = self.generate_read_counts(filtered_bam)
        if not readcount_file:
            return None
        
        # Step 3: Run IchorCNA
        ichorcna_success = self.run_ichorcna(readcount_file)
        if not ichorcna_success:
            return None

        # Step 4: CNA instability stopgap — set TF=0 if <15% genome is altered
        cna_fraction, adjusted_tf = self.apply_cna_instability_filter(input_bam.stem)

        # Always generate a TF override plot per sample for cross-sample viewing.
        self._regenerate_plot_with_zero_tf(
            input_bam.stem,
            cna_fraction=cna_fraction,
            force_tf_zero=(adjusted_tf == 0),
        )

        return {
            'sample': input_bam.stem,
            'filtered_bam': filtered_bam,
            'readcount_file': readcount_file,
            'ichorcna_success': ichorcna_success,
            'cna_fraction': cna_fraction,
            'adjusted_tf': adjusted_tf,
        }
    
    def run_pipeline(self):
        """Run the complete IchorCNA pipeline on all BAM files."""
        self.logger.info("Starting IchorCNA pipeline in WSL...")
        
        # Check dependencies
        if not self.check_dependencies():
            self.logger.error("Dependency check failed. Please install required tools.")
            return False
        
        # Find BAM files
        bam_files = self.find_bam_files()
        if not bam_files:
            self.logger.error("No BAM files found in input directory")
            return False
        
        # Process each BAM file
        results = []
        successful = 0
        
        for bam_file in bam_files:
            try:
                result = self.process_sample(bam_file)
                if result:
                    results.append(result)
                    successful += 1
                    self.logger.info(f"✓ Successfully processed {bam_file.name}")
                else:
                    self.logger.error(f"✗ Failed to process {bam_file.name}")
            except Exception as e:
                self.logger.error(f"Error processing {bam_file.name}: {e}")
        
        # Create summary report
        self.create_summary_report(results)

        # Create a single cross-sample view of TF0 override plots.
        self.create_tf0_override_collective_file(results)

        cna_passed = sum(
            1 for r in results
            if r.get('cna_fraction') is not None and r['cna_fraction'] >= 0.15
        )
        self.logger.info(f"Pipeline completed: {successful}/{len(bam_files)} samples processed successfully")
        self.logger.info(f"CNA stability filter: {cna_passed}/{successful} samples passed (≥15% genome altered — TF considered reliable)")
        return successful > 0
    
    def create_summary_report(self, results):
        """Create a summary report of the pipeline results."""
        summary_file = self.output_dir / "pipeline_summary.txt"
        
        summary_content = f"""IchorCNA Pipeline Summary
==================================================

Input directory: {self.input_dir}
Output directory: {self.output_dir}
Processed samples: {len(results)}
    Normal restart vector used: {self.ichorcna_params['normal']}

IchorCNA Parameters Used:
"""
        for key, value in self.ichorcna_params.items():
            summary_content += f"  {key}: {value}\n"

        summary_content += "\nPreprocessing Parameters:\n"
        summary_content += f"  min_mapq: {self.filter_params['min_mapq']}\n"
        summary_content += f"  max_read_length: {self.filter_params['max_read_length']}\n"
        summary_content += f"  readcount_window: {self.readcount_params['window']}\n"
        summary_content += f"  readcount_quality: {self.readcount_params['quality']}\n"
        
        summary_content += "\nProcessed Samples:\n"
        for result in results:
            cna_pct = (
                f"{result['cna_fraction'] * 100:.1f}%"
                if result.get('cna_fraction') is not None
                else "N/A"
            )
            tf_note = (
                "TF set to 0 (CNA fraction <15% — unstable)"
                if result.get('adjusted_tf') == 0
                else ("CNA filter skipped" if result.get('cna_fraction') is None else "TF retained")
            )
            summary_content += f"  - {result['sample']}  |  Genome altered: {cna_pct}  |  {tf_note}\n"
        
        summary_content += f"""
Output Files:
  - Filtered BAMs: {self.filtered_dir}
  - Read counts: {self.readcount_dir}
  - IchorCNA results: {self.ichorcna_dir}
    - TF0 override gallery (all samples): {self.ichorcna_dir / 'genomeWide_TF0_override_all_samples.html'}
    - TF0 override merged PDF (all samples, if pdfunite available): {self.ichorcna_dir / 'genomeWide_TF0_override_all_samples.pdf'}
  - Logs: {self.logs_dir}
"""
        
        # Write using WSL
        write_cmd = f'cat > {summary_file} << \'EOF\'\n{summary_content}\nEOF'
        self.run_wsl_command(write_cmd)
        
        self.logger.info(f"Summary report saved: {summary_file}")

    def _find_tf0_override_plots(self, sample_names):
        """Locate per-sample TF0 override plot files in a stable order."""
        candidates = []

        for sample_name in sample_names:
            sample_plot_dir = self.ichorcna_dir / sample_name / sample_name
            matches = sorted(sample_plot_dir.glob(f"{sample_name}_genomeWide_TF0_override.pdf"))
            if not matches:
                # Fallback for minor naming differences or additional suffixes.
                matches = sorted(sample_plot_dir.glob(f"{sample_name}_*TF0_override*.pdf"))
            if matches:
                candidates.append(matches[0])

        return candidates

    def _write_tf0_override_gallery_html(self, tf0_plots):
        """Write one HTML file that embeds all TF0 override PDFs for quick browsing."""
        gallery_file = self.ichorcna_dir / "genomeWide_TF0_override_all_samples.html"

        sections = []
        for plot_file in tf0_plots:
            sample_name = plot_file.parent.name
            rel_plot_path = os.path.relpath(plot_file, self.ichorcna_dir)
            section = f"""
<section class=\"card\">
  <h2>{html.escape(sample_name)}</h2>
  <object data=\"{html.escape(rel_plot_path)}\" type=\"application/pdf\" width=\"100%\" height=\"650\">
    <a href=\"{html.escape(rel_plot_path)}\">Open {html.escape(sample_name)} TF0 override PDF</a>
  </object>
</section>
"""
            sections.append(section)

        html_content = f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <title>genomeWide TF0 Override Plots</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; background: #f8fafc; color: #0f172a; }}
    h1 {{ margin-bottom: 20px; }}
    .grid {{ display: grid; gap: 20px; }}
    .card {{ background: #ffffff; border: 1px solid #dbe3ef; border-radius: 10px; padding: 12px; box-shadow: 0 1px 2px rgba(0,0,0,0.06); }}
    h2 {{ margin: 0 0 10px; font-size: 18px; }}
  </style>
</head>
<body>
  <h1>genomeWide TF0 Override Plots ({len(tf0_plots)} samples)</h1>
  <div class=\"grid\">
    {''.join(sections)}
  </div>
</body>
</html>
"""

        gallery_file.write_text(html_content)
        self.logger.info(f"✓ TF0 override gallery written: {gallery_file}")
        return gallery_file

    def _merge_tf0_override_pdfs(self, tf0_plots):
        """Attempt to merge TF0 override PDFs into one multi-page PDF using pdfunite."""
        merged_pdf = self.ichorcna_dir / "genomeWide_TF0_override_all_samples.pdf"

        if len(tf0_plots) == 1:
            shutil.copy2(tf0_plots[0], merged_pdf)
            self.logger.info(f"✓ Single TF0 override plot copied to: {merged_pdf}")
            return merged_pdf

        pdfunite_check = self.run_wsl_command("which pdfunite")
        if pdfunite_check.returncode != 0:
            self.logger.warning(
                "pdfunite is not available; skipping merged PDF generation (HTML gallery still created)."
            )
            return None

        quoted_inputs = ' '.join(shlex.quote(str(p)) for p in tf0_plots)
        cmd = f"pdfunite {quoted_inputs} {shlex.quote(str(merged_pdf))}"
        result = self.run_wsl_command(cmd)
        if result.returncode == 0:
            self.logger.info(f"✓ Merged TF0 override PDF created: {merged_pdf}")
            return merged_pdf

        self.logger.warning(
            f"Could not merge TF0 override PDFs with pdfunite: {result.stderr[:300]}"
        )
        return None

    def create_tf0_override_collective_file(self, results):
        """
        Create collective cross-sample TF0 override output files for easy review.

        Always writes an HTML gallery and also writes a merged multi-page PDF when
        pdfunite is available.
        """
        sample_names = [r['sample'] for r in results]
        tf0_plots = self._find_tf0_override_plots(sample_names)

        if not tf0_plots:
            self.logger.warning(
                "No genomeWide_TF0_override plots found; collective file was not created."
            )
            return

        self._write_tf0_override_gallery_html(tf0_plots)
        self._merge_tf0_override_pdfs(tf0_plots)

def main():
    """Main function for WSL pipeline."""
    import argparse
    
    parser = argparse.ArgumentParser(description='IchorCNA pipeline for Nanopore BAM files (WSL)')
    parser.add_argument('--input-dir', required=True, help='Directory containing BAM files')
    parser.add_argument('--output-dir', required=True, help='Output directory (WSL path)')
    parser.add_argument('--reference', help='Reference genome file (optional)')
    parser.add_argument(
        '--normal-restarts',
        default='c(0.5,0.8)',
        help="ichorCNA normal restart vector, e.g. 'c(0.5,0.8)' or 'c(0.8)'",
    )
    
    args = parser.parse_args()
    
    # Initialize and run pipeline
    pipeline = IchorCNAPipelineWSL(
        args.input_dir,
        args.output_dir,
        args.reference,
        normal_restarts=args.normal_restarts,
    )
    success = pipeline.run_pipeline()
    
    if success:
        print("Pipeline completed successfully!")
    else:
        print("Pipeline failed. Check logs for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()
