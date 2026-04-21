This directory creates a heatmap depicting the most common 4-mer motifs on cfDNA 5' ends. The script requires several dependencies, mostly R packages and uses a modified
version of codes from the fragmentomics_GenomBiol (https://github.com/Puputnik/Fragmentomics_GenomBiol) package which you want to take a look at.

In order to simplify I made some changes to their scripts and methods, for your convenience I will include those changes in this directory. The codes in this directory do not rely on the fragmentomics_GenomBiol package being downloaded to work. 
Still credit for the method goes to them, make sure to acredit them in all intended publications featuring data generated with this. 

In order to run this pipeline you will have to change the hardcoded paths in the run_fragmentomics_analysis.sh or run_fragmentomics_analysis_softclip.sh file to match your environment. These scripts work as wrappers and will run the whole analysis when properly set up. The run_fragmentoimcs_analysis.sh wrapper is recommended in most cases, use run_fragmentomics_analysis_softclip.sh if you get suspicious results and want to see if that version leads to improvement, if the softclip version gives bettere results you can likely conclude there is some sort of trimming issue in your bam files that should be looked into. For best performanced download the "4-mer_heatmap" directory as a whole to some location and run the script from there. The run_fragmentomic_analysis.sh requires the rest of the directory as it is a wrapper that runs the nanopore_fragmentomics.py script to create .motif and .stats files from your bam files and then runs count_motif.R to creat .motif.csv files for analysis. You can then run the generate_heatmaps.py, 4-mer_PERMANOVA.py and 4-mer_PCA_analysis.py scripts on the directory containing the .motif.csv files for statistical analysis and visualization of the results. When properly set up simply run "bash run_fragmentomic_analysis.sh" on a linux shell to 

The code expects a bam file input, a matching .fa index to the bam files and a chromosome list so make sure to have those on hand. The scripts in this directory have been tested with the following dependencies

Python 3.9+
  - pandas (v2.2.3)
  - numpy (v2.0.2)
  - pysam (v0.23.0)
R (v4.3.3)
  - data.table (v1.18.2.1)
Used for statistics scripts only
  - scipy (v1.13.1)
  - sklearn (v1.6.1)
