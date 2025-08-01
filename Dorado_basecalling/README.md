# Using Dorado to basecall pod5/fast5 files

Most pipelines made here were create based on multiplexed samples. Singleplex samples will work as well or better than the multiplex variants. 
If multiplexing remember to run the Dorado demux command and specify the kit you used. For alignment i would recommend using either the telomere to telomere consortium (T2T)
or the hg38 genome. T2T has less holes and is more robust overall making it better for telomeric analysis, hg38 however has more citations backing it and is more comfortable to use with some tools used for this method.
I will also post a sbatch wrapper here which details the manner in which i did the run. It's a 3 stage operation that basecalls, demultiplexes and aligns in one run.
It's also possible to finish each step separately. This pipeline is recommended to run on a cluster as it is really intensive and is used to process really big data. Can also be performed using the data acquisition unit but it won't be able to rebasecall in case you want to align the raw data to a different genome or for some other reason. 

## KEEP IN MIND ONLY THE "Dorado basecall" COMMAND NEEDS TO RUN IN A GPU IF PERFORMING SEPARATELY RUN THE OTHER COMMANDS IN A CPU ENVIRONMENT UNLESS YOU HAVE A POWERFUL GPU READILY AVAILABLE.

A downloadable Dorado executable can be found here: https://github.com/nanoporetech/dorado

A downloadable version of all genomes mentioned above can be found here: https://hgdownload.soe.ucsc.edu/downloads.html
(When downloading ref. genomes from USCS keep in mind you will be receiveing tar.gz files which are best to unseal using the gunzip command in a Linux environment. If on windows you can use Ubuntu/WSL. Once unzipped the fasta files in the directory are readily accessible in all systems.)
