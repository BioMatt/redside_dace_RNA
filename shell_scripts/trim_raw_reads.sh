#!/bin/bash
#SBATCH --time=07:00:00
#SBATCH --account=def-kmj477
#SBATCH --array=0-29
#SBATCH --mem=96000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL


module load nixpkgs/16.09
module load trimmomatic/0.36

dir=/home/biomatt/scratch/dace/raw_reads/

out_dir=/home/biomatt/scratch/dace/trimmed_reads


forward_files=(${dir}*_R1.fastq.gz)
reverse_files=(${dir}*_R2.fastq.gz)

forward_file=${forward_files[${SLURM_ARRAY_TASK_ID}]}
reverse_file=${reverse_files[${SLURM_ARRAY_TASK_ID}]}

echo ${forward_file}
echo ${reverse_file}

name=`basename ${forward_file} _R1.fastq.gz`

echo ${name}

sample=$(echo ${name} | cut -d "." -f5)

echo ${sample}

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE \
${forward_file} ${reverse_file} \
${out_dir}/paired/${sample}_R1_paired.fq.gz ${out_dir}/unpaired/${sample}_R1_unpaired.fq.gz \
${out_dir}/paired/${sample}_R2_paired.fq.gz ${out_dir}/unpaired/${sample}_R2_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:36