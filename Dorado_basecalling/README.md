# Using Dorado to basecall pod5/fast5 files

Most pipelines made here were create based on multiplexed samples. Singleplex samples will work as well or better than the multiplex variants. 
If multiplexing remember to run the Dorado demux command and specify the kit you used. For alignment i would recommend using either the telomere to telomere consortium (T2T)
or the hg38 genome. T2T has less holes and is more robust overall making it better for telomeric analysis, hg38 however has more citations backing it and is more comfortable to use with some tools used for this method.
I will also post a sbatch wrapper here which details the manner in which i did the run. It's a 3 stage operation that basecalls, demultiplexes and aligns in one run.
It's also possible to finish each step separately. This pipeline is recommended to run on a cluste as it is really intensive and is used to process really big data.

## KEEP IN MIND ONLY THE "Dorado basecall" COMMAND NEEDS TO RUN IN A GPU IF PERFORMING SEPARATELY RUN THE OTHER COMMANDS IN A CPU ENVIRONMENT UNLESS YOU HAVE A POWERFUL GPU READILY AVAILABLE.
