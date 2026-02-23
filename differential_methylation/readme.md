This code is intended for differential methylation analysis on Nanopore methylation data processed from POD5 or fast5 files through a pipeline of Dorado for basecalling and modkit for methylation value extraction. 

The pipeline is long but in essence should be easy to use. Currently it's adjusted to BRCA2 samples, it will be necessary to change the input functions to match your input bed file locations and grouping strategy. The analysis is built to compare 2 groups in this case BRCA2 and control it expects you to provide data for both groups there is no inbuilt baseline in this code so you must have a valid and matching ctrl dataset to compare with your actual focus group. This will be the most complicated part of the pipeline to adapt to your needs

The script analyses promoter regions, this can be changed according to your needs. It uses an ensembl API to extract promoter region information, at the top of the script there is a configuration function here you can edit your genese of interest simply list genes whose promoters you would like to analyze as so:
GENES_OF_INTEREST = {
    'ACTL6A': None, 'APC': None, 'ATR': None, 'BRCA1': None, 'BRCA2': None,
    'CCND2': None
}
Also as part of the config there are 3 parameters you can set MIN_CPGS = X, MIN_SAMPLES = X and PERMUTATIONS = XXXXXX. MIN_CPGS should be set according to the depth of coverage in your data the higher you can make this number the better, it is based on the Nvalid_cov value in your modkit pileup processed bed files, this essentially shows amount of CpGs in the promoter region per sample Ex. if a promoter region has 20 CpGs and each CpG is covered by 10 threads in a sample file that sample would be able to pass a MIN_CPGS value set to 200 or lower. Realistically when running this type of analysis on cfDNA that is Nanopore sequenced your min_cpg will have to be much lower if a large amount of your samples consistently pass with MIN_CPGS set to 10 that is good and set to 20 would be great. MIN_SAMPLES is the minimum number of samples per group that need to pass the min_cpgs check for a gene in order to include that gene in analysis. Prevents inaccurate analysis I recommend setting this value to at least half of the total sample amount per group in your dataset. PERMUTATIONS refers to how many permutations will be performed for P-value calculation. Higher value will lead to higher quality of analysis but will also increase run time, permutations = 10000 is recommended but setting it lower is acceptable. 

It is also possible to control parameters via command line a full possible command line is: 
**"python analyze_methylation_per_sample.py /path/to/your/bed_dir --output-dir /path/to/your/out_dir --min-cpgs x --min-samples XX --permutations XXXXX"**
A bed_dir path is necessary needs to be a single path that leads to all bed files of both groups. Grouping strategy must be defined in script, --output dir is not needed but is recommended, othewise output will go to a default directory of ./methylation_results. 

Dependencies for this script are: pands, numpy, matplotlib and requests. These are all python packages install with pip or bioconda if needed.

Currently code is not fully tested am uploading it here simply for my own convenience


