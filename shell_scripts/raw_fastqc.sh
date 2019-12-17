#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=10000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-60

module load nixpkgs/16.09
module load fastqc/0.11.8

raw_read_list=raw_dace_reads.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${raw_read_list}" 
read=$($string) 

fastqc $read --outdir /home/biomatt/scratch/dace/raw_fastqc_output