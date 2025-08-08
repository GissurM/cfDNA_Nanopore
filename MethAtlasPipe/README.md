# Meth atlas pipeline
This pipeline relies on a tool known as MethAtlas to perform methylome deconvolution on Nanopore data the tool can be found here: https://github.com/nloyfer/meth_atlas. 
### Disclaimer: I personally had to make slight changes to the code of the tool to make it work in my environment. In particular I had to update the names of some python packages to use the newer version of such a package.
If you do not know how to do this it should be very easy to use a AI coding assistant to do it for you simply specify your current python version and show it the code for the tool. I will not post the edited version of the code i made 
due to distribution guidelines detailed on the github site of it's original creators, the tool is the intellectual property of Yissum Research Development Company Of The Hebrew University of Jerusalem Ltd. and should be treated as such in all cases. 
The code is however free to use and you can feel free to make any edits or changes without facing legal complications.

This is also not the first time Nanopore data has been prepped to be used for MethAtlas and an older pipeline i took a lot of inspiration from can be found here: https://github.com/methylgrammarlab/cfdna-ont/tree/main/deconvolution_code. 
This Github page was released as part of the following publication i would recommend reading as it has been hugely helpful with my data analysis: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02710-1. 
However the only thing you will require from the github page posted above is the manifest found under "Infinium_HumanMethylation450k_manifests" on the page. 
There are two manifests one for hg38 and one for hg19 I personally used hg38 as that was the genome I had aligned my data to, should you insist upon using T2T or other reference fastas I would recommend creating your own manifest following the same format. 
Where my pipeline however seperates from the older pipeline is in terms of what tools I will utilize for further analysis. 

The pipeline mentioned above utilizes guppy for basecalling and than megalodon for methylation analysis. These tools are currently considered inferior to Dorado and Modkit since Dorado has been touted as the gold standard methylation calling tool by Oxford Nanopore Technologies since 2024. 
As to how to run the pipelines simply basecall your files and then transform them to a BED format as described in the Dorado and Modkit chapters in this repository. Then run your Modkit created BED files through the python code I created which will be posted in this repository. 
This should convert the methylation values in your BED files to the beta value format given by the 450K illumina array. Remember to set the manifest variable in the python code to the path of the manifest you downloaded from: https://github.com/methylgrammarlab/cfdna-ont/tree/main/deconvolution_code. 
After doing this you should be able to run the MethAtlas tool using the following command line: "python deconvolve.py -a reference_atlas.csv C:\path\to\your\array-transformed-nanopore-data.csv". I would recommend running this within the meth_atlas directory you should have git cloned from: https://github.com/nloyfer/meth_atlas 

This is the stage where it gets a little bit annoying to set up your virtual environment for further analysis. First of all you are going to want to set up Python on your path. On windows this can be done easily by finding the newest Python package in the windows store, python can also be downloaded from python.org.
If you get the option make sure to check the box to install pip along with python. Once downloaded check if you have both installed via python --version and pip --version. Now you can pip install pandas and numpy via the command line "pip install pandas numpy" these will be necessary for this analysis. 
You should than install git so you can git clone MethAtlas. You can install git from: https://git-scm.com/downloads, you should make sure to set it into your path it should be possible to do via the installer so just read through it and make sure not to miss checking any boxes that mention path. 
Then you can run the git command directly in cmd or powershell and run "git clone https://github.com/nloyfer/meth_atlas" this will set up a directory named meth_atlas on your computer. Git can install stuff to unpredictable places but it is likely in your user home directory or program files. This directory will include everything you need to run MethAtlas and i would always recommend running the 
"python deconvolve.py -a reference_atlas.csv C:\path\to\your\array-transformed-nanopore-data.csv" command line from within this directory as it is more convenient. This will also ensure your deconvolution results will be in this directory making them easy tofind. 
This entire pipeline aside from Dorado and modkit to start does not require high processing and can be ran locally on most decent hardware setups. When pip installing I would recommend doing so in a python environment best to create via conda. 

## Dependencies
- Python (most versions over 3.0 should be compatible, some things might need to be adapted based on your python version.)
- Deconvolve.py (this is a program that should be included in the meth_atlas directory you git cloned)
- "Infinium_HumanMethylation450k_manifests" (either hg38, hg19 or custom)
- Python packages:
    - pandas
    - numpy

