# cfDNA thread size analysis 
This method is fairly simple similar to the chromsome counts code as it uses a flag which dorado automatically inserts onto the header of each individual cfDNA thread. 
The flag in question is called "qwidth" and is something i believe to be a universal flag across NGS sequencing and not unique to Nanopore. The "qwidth" flag is representitive of the length of the sequenced read. 
In a normal case this would not be considered the optimal flag for cfDNA size analysis from sequenced data due to it only representing one strand of the cfDNA thread, however both Nanopore cfDNA protocols
found on the Nanoporetech.com website as "Ligation sequencing V14 — Human cfDNA singleplex (SQK-LSK114)" and "Ligation sequencing V14 — Human cfDNA multiplex (SQK-NBD114.24)" include an end prep step which repairs any jagged edges. 
If this step is not performed or otherwise went awry the qwidth flag cannot be considered a true represntetive of the cfDNA thread size. The reason why the qwidth flag needs to be used is that
the isize flag which represents the full length of both strands of a double stranded cfDNA thread does not appear in Nanopore data due to it's unique nature as a single stranded DNA thread sequencing method. 
If performing this analysis on Illumina data or other dsDNA sequencing methods the "isize" flag should be utilized instead. 

Another thing to watch out for is that if the --no-trim flag is used as parts of both dorado basecall and dorado demux commands you will end up with bam files showing cfDNA threads with 
barcodes/native adapters on both ends this will cause the lengths to be incorrectly represented. To get around it run the dorado basecall command with the --no-trim flag as the barcode data is still necessary for demux and then run 
dorado demux without the --no-trim flag, this will ensure your final demultiplexed results will be in the correct sizes. This will also be necessary for other methods detailed in this repository and will generally cause a cleaner dataset for you to use.

The scripts themselves are fairly simple to run first thread_size_analysis.r is used to create individual histograms for each bam file containing information on cfDNA thread sizes. To run this script change  the bam_dir parameter to a path containing bam files and then run the script using R. For a more detailed output run extract_detailed_cfDNA_sizes.py you can do so via a command line "python extract_detailed_cfDNA_sizes.py --bam-dir /your/bam/dir --output-dir /your/output/dir" there are additional parameters you can change but running it like this will be successful alternitively it's possible to change the default directories in the script if preferred. Once this script has finished running you should now have size bins for each of the mains fragment size groups, mono, di, tri and HMW. Run the cfDNA_size_KS_and_MWU.py and cfDNA_size_PCA.py scripts on the output from the previous script to receive interpretable statistics and graphs.

## Dependencies
- R studios - version 4.3.3 (this version is tested others may work as well) 
- R packages:
    - Rsamtools (v2.18.0)
    - ggplot (v4.0.2)
    - dplyr  (v1.2.0)
- Python - version 3.9+
  - pandas (v2.2.3)
  - numpy (v2.0.2)
  - scipy (v1.13.1)
  - sckit-learn (v1.6.1)
  - pysam (v0.23.0)
  - seaborn (v0.13.2)
