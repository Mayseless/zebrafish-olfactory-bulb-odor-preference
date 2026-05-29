#!/bin/bash

#SBATCH --job-name=arr_analyze_FLASHseq                   #This is the name of your job
#SBATCH --cpus-per-task=2                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16G              #This is the memory reserved per core.
#Total memory reserved: 32GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=./log_files/myoutput%j.out     #These are the STDOUT and STDERR files
#SBATCH --error=./log_files/myerror%j.err

#You selected an array of jobs from 1 to 3075
#SBATCH --array=1-3075

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


# Set the base directory
p_dir="/scicore/home/schiera/maysel0000/20240412_CampariFS_raw"
samples_file="${p_dir}/samples_to_process.txt"

# Ensure necessary directories are created
mkdir -p ${p_dir}/merged
mkdir -p ${p_dir}/processed_data/STAR


# Read the directory and identifier for this job using the SLURM_ARRAY_TASK_ID
sample_info=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $samples_file)
sample_dir=$(echo $sample_info | cut -d':' -f1)
full_identifier=$(echo $sample_info | cut -d':' -f2)

# Skip processing for the first line if it doesn't contain a full identifier
if [ -z "$full_identifier" ]; then
    echo "Skipping invalid entry: $sample_info"
    exit 0
fi

echo "Processing directory: $sample_dir"
echo "Full identifier: $full_identifier"

R1="${p_dir}/merged/${full_identifier}_merged_R1_001.fastq.gz"

# Check if the merged file exists and merge if it does not
if [ ! -f "$R1" ]; then
    echo "Merging FastQ files for $full_identifier"
    cat ${sample_dir}/*_R1_*.fastq.gz > $R1
fi

module load STAR/2.7.1a-foss-2018b
STAR_REF="${p_dir}/STAR_indexed_genome103/"
STAR --runThreadN 4 --genomeDir $STAR_REF --readFilesIn $R1 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${p_dir}/processed_data/STAR/${full_identifier}_

