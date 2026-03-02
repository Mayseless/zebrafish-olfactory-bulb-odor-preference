#!/bin/bash

#SBATCH --job-name=OM_20230131_CaMPARI_Cad_FS #This is the name of your job

#SBATCH --cpus-per-task=4 #This is the number of cores reserved, 20 is max for one node.

#SBATCH --mem-per-cpu=16G #This is the memory reserved per core. 256G is max total memory for one node.

#SBATCH --time=24:00:00 #This is the time that your task will run

#SBATCH --qos=1day #You will run in this queue

#SBATCH --output=myrun.o%j #These are the STDOUT and STDERR files

#SBATCH --error=myrun.e%j

#SBATCH --mail-type=END,FAIL,TIME_LIMIT

#SBATCH --mail-user=oded.mayseless@unibas.ch #You will be notified via email when your task ends or fails

 

module load Subread/1.6.4-foss-2018b

 

GTF_path="/scicore/home/schiera/maysel0000/20240412_CampariFS_raw/Danio_rerio.GRCz11.103.gtf"

out_dir="/scicore/home/schiera/maysel0000/20240412_CampariFS_raw/processed_data/"

out_file="${out_dir}featureCounts/count_gene.txt"

mkdir -p "${out_dir}featureCounts/"

## since gene_name id not always in the same position (it can behind gene_id, or behind transcript_id, or behind exon_number) in column 9 of the ## gtf file. so, the featureCounts can't get the gene_name information. We use gene_id instead since it's position is the always the same in column 9

featureCounts -T 4 -B --primary -t gene -g gene_id -a ${GTF_path} -o ${out_file} ${out_dir}STAR/*.bam
