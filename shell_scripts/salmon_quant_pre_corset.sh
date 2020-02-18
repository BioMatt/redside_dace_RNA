#!/bin/bash
#SBATCH --time=18:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=0
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --constraint=skylake
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-30

module load nixpkgs/16.09 gcc/7.3.0 openmpi/3.1.2
module load salmon/1.1.0

list=array_paired_trimmed_reads.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${list}" 
str=$($string) 

var=$(echo $str | awk -F"\t" '{print $1, $2}') 
set -- $var 

c1=$1 
c2=$2

echo "$c1" 
echo "$c2"

salmon quant -i /home/biomatt/scratch/dace/redside_dace_precorset_index -l IU --dumpEq \
-1 ${c1} \
-2 ${c2} \
--validateMappings \
-o /home/biomatt/scratch/dace/salm_outputs_precorset/salm_output_${c1}