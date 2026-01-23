This file creates a heatmap depicting the most common 4-mer motifs on cfDNA 5' ends. The script requires several dependencies, mostly R packages and uses a modified
version of codes from the fragmentomics_GenomBiol (https://github.com/Puputnik/Fragmentomics_GenomBiol) package which you want to take a look at.

In order to simplify I made some changes to their scripts and methods, for your convenience I will include those changes in this directory. The codes in this directory do not rely on the fragmentomics_GenomBiol package being downloaded to work. 
Still credit for the method goes to them, make sure to acredit them in all intended publications featuring data generated with this. 

The code expects a bam file input, a matching .fa index to the bam files and a chromosome list so make sure to have those on hand.
