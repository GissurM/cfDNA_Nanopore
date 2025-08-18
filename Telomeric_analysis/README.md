# Telomeric analysis from cfDNA
## Disclaimer: While I believe the logic behind this method is mostly sound it seems to be inaccurate as of now. Analyze results with a grain of salt. There is a lack of pure telomeric cfDNA threads. This is likely due to there being no unique binding spots for these fragments on the genome. Currently looking for a workaround would rather use aligned data if possible.

This method is created by myself and gives a result in the form of total telomeric DNA % of cfDNA in the blood stream. It is very hard to fact check whether the results of this particular method are correct. 
The current other methods used to detect similar metrics qPCR/dPCR are not directly translateable as they give metrics such as a T/S ratio which is not directly translatable to a total %. 
In order to fact check how correct this is i have been comparing results with a previously performed qPCR analysis performed on the samples and checking for similar ratios of results. 
Currently I need more power for a full confirmation. 

To find this percentage i created a rule based algorhithm which defines several functions that are representitive of telomeric regions. These rules are set up in a way to attempt to balance
sensitivity and specificity. The major rules are that there needs to be a motif matching the forward or reverse motif of telomeric repeat sequences. The first appearance of this motif needs to be within 15 bp of either end of the cfDNA thread, 
the motif than has to repeat at least 4 times with a certain density requirement. I also account for alternative Telomeric motifs found here: https://pmc.ncbi.nlm.nih.gov/articles/PMC8256856. The motifs also have to be a certain percentage of the thread
so larger threads need longer motif repeats. There is a 8% error rate allowance in the code to account for the known basecalling error rat of Nanopore sequencing. 

I then created a test function with several positive and negative controls synthesized with an AI. The controls are made to be roughly as long as the average cfDNA thread and both positive and negative controls include both surefire cases which should be easy to discover as well as complete edge cases which the algorhithm should struggle with. 
I ran the test function a bunch of times and tweaked the algorithm until the test run with about 90% accuracy with a relative balance between sensitivity and specificity. 
It is likely possible to make the algorhithm more accurate but it will take some work. 3 major options come to mind first is to add additional rules I currently have no real idea what rules those could be,
second it should be possible to add additional context such as accounting for subtelomeric regions (hard to do as they do not contain motifs) or maybe find some methylation pattern unique to telomeres.
Thirdly it would be possible to train an AI model either KNN or something like a random forest model to get really good at detecting these telomeric threads, this would both be really time consuming and would require a lot of data but would definitly be worth doing.
When running this code i get consistent results at around 0.0010-0.0030% this is a very low percentage but it is known that telomeric DNA is rather underrepresented in blood under normal conditions. 
The ratios however are seemingly consistent with the theories that older people have shorter telomeres and that BRCA2-999del5 carriers have shorter telomeres on average currently I don't have enough power to state this as fact.

I would recommend running this code in a cluster as it is fairly resource intensive due to having to read throuh a lot of bases in a Nanopore bam file which is a lot of data. 
It's not too difficult to set up and i will also post a wrapper which should allow you to set up a proper virtual environment within the HPC and than clean up afterwards. 

## Dependencies
- Python 3.9 (tested)
- Samtools (Needs to be available on system path, if possible running through module library is  recommended)
- Python packages:
    - Pysam
