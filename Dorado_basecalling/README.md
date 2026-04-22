# Using Dorado to basecall pod5/fast5 files

Most pipelines made here were create based on multiplexed samples. Singleplex samples will work as well or better than the multiplex variants. 
If multiplexing remember to run the Dorado demux command and specify the kit you used. For alignment i would recommend using either the telomere to telomere consortium (T2T)
or the hg38 genome. T2T has less holes and is more robust overall making it better for telomeric analysis, hg38 however has more citations backing it and is more comfortable to use with some tools used for this method.
I will also post a sbatch wrapper here which details the manner in which i did the run. It's a 3 stage operation that basecalls, demultiplexes and aligns in one run.
It's also possible to finish each step separately. This pipeline is recommended to run on a cluster as it is really intensive and is used to process really big data. Can also be performed using the PromethION data acquisition unit through the MinKNOW software but it will be inconvenient to re-basecall in case you want to align the raw data to a different reference fasta. Once the Basecall_wrapper.sh has been properly set up it should be possible to run it via "sbatch Basecall_wrapper.sh" if you are using an sbatch compatible system. Otherwise "bash Basecall_wrapper.sh" should work fine but again most machines cannot handle directly running dorado due to high processing costs so don't run this via WSL shell unless you are confident in your machines capabilities. 

The "dorado basecall" command requires a powerful GPU to run effectively as it uses machine learning based methods. The other functions "dorado demux" and "dorado aligner" run fine on a CPU. Given the size of the files being processed by this command and the complexity of the operation I would always recommend running this in a cluster.

A downloadable Dorado executable can be found here: https://github.com/nanoporetech/dorado

A downloadable version of all genomes mentioned above can be found here: https://hgdownload.soe.ucsc.edu/downloads.html
(When downloading ref. genomes from USCS keep in mind you will be receiveing tar.gz files which are best to unseal using the gunzip command in a Linux environment. If on windows you can use Ubuntu/WSL. Once unzipped the fasta files in the directory are readily accessible in all systems.)

Dependencies:
Dorado (v1.3.1) other versions will likely work but the wrapper may require slight alterations
