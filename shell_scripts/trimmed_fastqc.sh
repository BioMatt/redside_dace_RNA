#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=10000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-60

module load nixpkgs/16.09
module load fastqc/0.11.8

trimmed_read_list=trimmed_read_list.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${trimmed_read_list}" 
read=$($string) 

fastqc $read --outdir /home/biomatt/scratch/dace/trimmed_fastqc_output