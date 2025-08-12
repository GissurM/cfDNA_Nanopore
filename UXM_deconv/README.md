# UXM_deconv
This pipeline is still relatively unfit for use. 

Uses the UXM and wgbs tools which are both necessary for this analysis, these tools are also attributed to Netanel Loyfer sem as MethAtlas. As this analysis relies in very small part on any work I have done the only additional content in this directory
will be a wrapper script I made to be ran on modbam files made by Nanopore (modbam is a term for bam files with epigenetic modification data, therefore assume all bam files from Nanopore are modbam unless you have removed modification data in one way or another).
wgbstools and UXM_deconv can be found through the following links: https://github.com/nloyfer/wgbs_tools, https://github.com/nloyfer/UXM_deconv. 

Using UXM_deconv is fairly straightforward as it is simply a .py python script with some additional commands available. However setting up wgbs_tools can be a problem. 
Wgbs_tools has the following dependencies (also listed at the link): 
- Python 3+
- Python libraries:
    - Pandas version 1.0 or newer
    - numpy
    - scipy
 - samtools
 - tabix / bgzip
 - bedtools (bedtools is listed as optional on the tools official website but you will need it for this analysis.)

You then need to set up the code onto your path using the setup.py code that you should get in the file you git clone from the website.
When you run the setup.py code using Python make sure to do so in a virtual environment with all the listed dependencies present otherwise it will not work, running it on linux or Mac is preferred as there are no easy ways to set up Samtools / bedtools on windows. 

## If the setup.py code gives error messages or complaints DO NOT make changes to the setup.py code or any other parts of wgbs_tools as it is unlikely to fix the issue. Rather uninstall all traces of wgbs_tools and then try again in the correct environment! Am speaking from personal experience here.

If the above statement sounds harsh that's because i altered the setup.py code and got the wgbstools command onto path but kept running into issues with individual parts of the wgbstools package and it took hours to sort through it all. I then gave up and deleted everything only for it all to work perfectly after running setup.py, if on windows make sure to use WSL using MSYS2 will only cause problems.
Anyways once you have gotten the wgbs_tools command onto your path the hard part is over. Then you can use bam2pat to transform bam files to pat.gz files. 

Currently the wrapper I use gives very suspect results that are unlikely to be biologically relevant. I am working on making it work better but am unsure if it will work. I don't expect many people to see this text considering this is a private repository for now but you know jsut talking to myself.
