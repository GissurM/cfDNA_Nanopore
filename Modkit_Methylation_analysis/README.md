# Modkit methylation analysis

Modkit is an extremely useful tool for methylation analysis that i have yet to fully delve into personally. 
For now we will utilize it's function to creath bed or bedgraph files from aligned bam files from Dorado. These bed files will then be used further downstream in this analysis. 
Modkit is mainly able to differentiate between 3 methylation types 5mC (normal methylated Cytosines), 5HmC (hydroxymethylated cytosines) and 6mA (methylated Adenosines). 
For current downstream analysis in this project only 5mC is necessary. Like Dorado Modkit is also a tool made by Oxford Nanopore technologies 
and is simple and intuitive to use once you understand it. In order to achieve our purposes we need only use the "modkit pileup" command with 
some additional parameters which i will detail in the wrapper script. Other commands are available and potentially useful, this page will be updated if i find a good use for them.

Modkit is available here as an easy to download and set up executable: https://github.com/nanoporetech/modkit
