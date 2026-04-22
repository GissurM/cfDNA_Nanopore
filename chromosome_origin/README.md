# Chromosomal origin of cfDNA

This is a custom analysis method that uses the header information given to each individual cfDNA thread as part of the dorado alignment process. Because of this the code requires genome aligned .bam files, you should know that your alignment was successful if the file includes both .bam and .bam.bai files. The bam.bai files may not be necessary in this analysis however it is better to have them in the directory with the .bam files than to not have them as samtools will sometimes rely on them. The code might work on animal genomes as well since it's designed to be broadly applicable but I have no current data on this.
To further explain each thread is given a flag by the likes of chr 1, chr  etc. there are also some abnormal headers, which i assume dorado likes using because there are a lot of them, so they will be placed in a seperate column in the analysis. In all likelihood this will be the tallest column.  

Chromosome_Histogram.py extracts the chromosome info and creates a histogram. To run simply run the command python Chromosome_histogram.py --counts-dir "/path/to/your/bam-file-dir". 

For a more detailed output run Chromosome_counts.r set the bam_dir and bam_files functions in the script to match your files than run the script with R. This script outputs a .csv for each sample which can be used with the chr_statistics.py script to get visual and statistics output. Set up the chr_statistics.py script to work in your environment and then run it using python. None of these scripts require high processing power but can take a while to run. 

## Dependencies
R version 4.5.0 (This version is tested others may work just as well)
R libraries:
   - Rsamtools (v2.18.0)
   - ggplot2 (v4.0.2)
Python (v3.9+)
   - numpy (v2.0.2)
   - pandas (v2.2.3)
   - scipy (v1.13.1)
   - matplotlib (v3.9.1)
