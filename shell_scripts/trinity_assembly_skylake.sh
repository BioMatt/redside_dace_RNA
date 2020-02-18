#!/bin/bash
#SBATCH --time=20-00:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=0
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --constraint=skylake
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --propagate=STACK

module load nixpkgs/16.09  gcc/7.3.0  openmpi/3.1.2 salmon/0.14.2-1
module load bowtie2/2.3.5.1 trinity/2.9.0 samtools/1.9
module load scipy-stack
pip list | grep -a numpy
module load java/1.8.0_192

Trinity --seqType fq --max_memory 50G --CPU 2 --samples_file /home/biomatt/scratch/dace/trinity_sample_list.txt \
--output /home/biomatt/scratch/dace/trinity_outputs_skylake3