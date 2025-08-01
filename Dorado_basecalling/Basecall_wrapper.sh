#!/bin/bash

#SBATCH --job-name=base
#SBATCH --partition=gpu-1xA100  # request node from a specific partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32      # 48 cores per node (96 in total)
#SBATCH --output=messages.align1-12.out.txt   
#SBATCH --error=messages.align1-12.err.txt 
 

#
# 1. Create a temporary directory with a unique identifier associated with your jobid
# 
scratchlocation=/scratch/users
if [ ! -d $scratchlocation/$USER ]; then
mkdir -p $scratchlocation/$USER
fi
tdir=$(mktemp -d $scratchlocation/$USER/$SLURM_JOB_ID-XXXX)
echo "the scratch dir is:"
echo $tdir

mkdir $tdir/myAoutputfiles1-12
mkdir $tdir/myAoutputfiles1-12/outputAdemux1-12
mkdir $tdir/outputdemux1-12
 
# Go to the temporary directory
cd $tdir
# Make sure to copy the files you need for the job
time cp -rf /path/to/your/dorado/executable $tdir/. # Not needed if running a module

time cp -rf /path/to/your/fast5/files $tdir/.

time cp -rf /path/to/your/reference/fasta $tdir/.

# Execute Dorado commands
time dorado-0.9.1-linux-x64/bin/dorado basecaller sup@v5.0.0 --kit-name your-kit-name name-of-pod5-directory --modified-bases 5mCG_5hmCG -r --no-trim -v --device cuda:all > $tdir/myAoutputfiles1-12/FC1-12Goutput.bam # When specifying kit name make sure you use "-" instead of "_" if your kit name has a _ in it
 
time dorado-0.9.1-linux-x64/bin/dorado demux --kit-name your-kit-name $tdir/myAoutputfiles1-12/FC1-12Goutput.bam --output-dir $tdir/outputdemux1-12 -v # Don't run --no-trim barcodes unnecessary after this point and get in the way of true size analysis of cfDNA
 
time dorado-0.9.1-linux-x64/bin/dorado aligner chm13v2.0_unmasked.fa $tdir/outputdemux1-12 -r --output-dir $tdir/myAoutputfiles1-12/outputAdemux1-12 -v
# After the job is completed make sure to copy the output to your submit directory.
time cp -r  $tdir/myAoutputfiles1-12 /hpcdata/Mimir/gmi6/STRBasecall_output1-12/.

 
# If the program produces many output files you can add a separate line for each file.
# Please try to only copy the files that you need.
 
#
# 5. clean scratch
#
echo ""
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date
