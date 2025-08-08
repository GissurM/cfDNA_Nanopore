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
Alternitavely it should be possible to run basecall, demux and alignment in one command line and skip the --no-trim command entirely however I could not get it to work and eventually gave up. If you can get it to work that's great! 

The code itself is fairly simple and does not rely on too many packages. It will output a individual barplot for each bam file in a directory. Like all codes i created for this analysis it was designed to run on a directory of files. 
If you only have one file simply set the path to the directory that file is contained in and it should work as intended. The code also outputs a comprehensive .csv file which compiles numerical data from all files into a row. 
Specifically it outputs mean, median and sd both with only fragments sized 100-1000 bp and with all fragments in the file, than it outputs how many threads in each file are representitve of the known thread sizes mono-, di- and trinucleosomal as well as high molecular weight DNA. 
It also outputs the ratio of the thread sizes in the file. This .csv file is created for easy access to most factors necessary for statistical analysis. Keep in mind that when working with Nanopore the absolute amount of cfDNA threads is often valueless, input amoutn of DNA is highly controlled by the nature of the analysis and is heavily normalized in ideal conditions, 
therefore the most valuable data to be gathered here is the average size of cfDNA fragments in bp and ratios of cfDNA sizes. The shape of the graph can also be valuable, if you see a sharp mononucleosomal peak and then two more softer mounds for di- and trinuclosomal peaks. The dinucleosomal peak should be slightly higher and the tri nucleosomal is often more spread out. 
Sometimes a fourth peak is visible possibly representing a tetranucleosomal peak which is a known and interesting phenomenon. 

## Dependencies
- R studios - version 4.5.0 (this version is tested others may work as well) 
- R packages:
    - Rsamtools
    - ggplot
    - dplyr 
