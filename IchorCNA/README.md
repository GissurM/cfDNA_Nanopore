This pipeline is made to run ichorCNA on low coverage Nanopore samples. To run this download the entire IchorCNA directory from this repository. Then make alterations to the INPUT_DIR, OUTPUT_DIR and NORMAL_RESTARTS functions to suit your needs. INPUT_DIR should be changed to the path containing your bam files, OUTPUT_DIR should be changed to the path you want to deposit the output files to and NORMAL_RESTARTS is fine as is but you cann add normal restart values from 0.0 to 1.0 if you would like. Adding a NORMAL_RESTART value above 0.8 will likely lead to unnecessarily high sensitivity leading to false positive tumor fraction detection if working with low coverage Nanopore data so woud not recommend it in that case. 

Once run_ichorCNA.py has been set up simply run it with "python run_ichorCNA.py" this script serves as a wrapper and will then proceed to run ichorCNA_pipeline.py and ichorCNA.r, the former ensures that the ichorCNA.r script will run with the correct settings and the latter is a script from the official ichorCNA repository: https://github.com/broadinstitute/ichorCNA. In this repository one minor change was added so it would fallback to a local reference fasta, but all credit for that script goes toward the original creators. 

Dependencies
Python (v3.0+)
  -pandas (v2.2.3)
R (v4.3.3)
  -ichorCNA (v0.3.2)
  -GenomeInfoDb (v1.38.8)
  -optparse (v1.7.5)
  -HMMcopy (v1.44.0)
Samtools (v1.19.2)
