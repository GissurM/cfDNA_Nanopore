# Dorado Basecalling

In this branch i will cover basecalling using Dorado. In essance this process is pretty simple. The Dorado application is easily accessible through open source means:
https://github.com/nanoporetech/dorado and is the Oxford Nanopore Technologies official basecaller as of now. For those who will perform the Nanopore analysis from scratch
the Nanopore basecaller and data acquisition unit themselves are usually fully able to basecall by themselves. It will be necessary to run the demux command on barcoded multiplexed
samples. Dorado should also be used to align data to a genome of your choice for the following pipelines I would recommend either the Telomere to Telomere (T2T) consortium 
genome which has less gaps around telomeric regions making it better for analysis of telomeric analysis, as an aside this genome might be higher quality overall but as of
now hg38 is more tested and likely more compatible with a lot of tools I would recommend using it for analyses like those performed with MethAtlas which is as of now incompatible with T2T. 
hg19 is also compatible with such tools but have not used it personally.
