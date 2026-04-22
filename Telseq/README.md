# Telseq post-processing scripts

These scripts are used to process results from the Telseq tool which can be found here: https://github.com/zd1/telseq. These scripts are highly personalized to function on my specific database and utilize a lot of specific information that I had access to for the samples I was working with, feel free to use the scripts as reference but these scripts are unlikely to work on an alternate dataset.

The pipeline for Telseq resutls analysis is as follows. First run telseq_parser.py on your telseq output this will transform it into a more maneagable and readable format.
After this run telseq_metadata_intergration.py on the telseq_parser.py output to inject metadata information into the dataset for statistical analysis. 
Lastly run telseq_visualization.py in the metadata containing .csv file created using the precvious script. This should generate an output containing several plots and statistics that can be used to interpret the telseq results. For these scripts either paths can be hardcoded into the script by editing the "default" function under each parser or by calling out each path in the command line such as by running "python telseq_parser.py --input-file /path/to/your file --output-csv /intended/output/file/location.csv etc." The telseq_parser.py should work on most datasets, the telseq_metadata_intergration.py script requiers a metadata.csv with information on, age, gender and BMI for your samples so will not be compatible without that information the script will not work without alteration. The telseq_visualization.py script requries output form telseq_metadata_intergration.py so will not work without the full pipeline. Should still be possible to alter it to work with datasets that don't include the metadata.

Dependencies
TelSeq (v0.0.1)
GCC (≥v4.8)
bamtools (v2.5.3)
python (v3.9+)
  -pandas (v2.2.3)
  -numpy (v2.0.2)
  -scipy (v1.13.1)
  -sckit-learn (v1.6.1)
