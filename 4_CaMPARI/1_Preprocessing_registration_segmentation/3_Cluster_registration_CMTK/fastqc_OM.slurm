#!/bin/bash

#SBATCH --job-name=arrayJob                   #This is the name of your job
#SBATCH --array=1-12
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --cpus-per-task=16                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=4G              #This is the memory reserved per core.


# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH -o hostname_%A_%a.out # Standard out goes to this file
#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname


#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


#load your required modules below
#################################


#export your required environment variables below
#################################################


#add your command lines below
module load Perl/5.28.0-foss-2018b
/scicore/home/schiera/GROUP/ImageRegistration/cmtk-3.2.3/build/bin/munger -b /scicore/home/schiera/GROUP/ImageRegistration/cmtk-3.2.3/build/bin/ -T 32 -awr 010203 -X 52 -G 80 -R 3 -A '--accuracy 0.8' -W '--accuracy 1.6' -s /scicore/home/schiera/maysel0000/registration/20230815_CaMPARI_draQ52500_ON/ref_20230620CamPARIE302_32_01.nrrd  -d .$SLURM_ARRAY_TASK_ID -v images/$SLURM_ARRAY_TASK_ID/

#module load fastqc/0.11.5-fasrc01

#fastqc --casava -o fastqc_reports A01_R1.pair.fastq.gz

~
~
~
~
~
~
~
~
~
"fastqc_OM.slurm" 40L, 1481C
