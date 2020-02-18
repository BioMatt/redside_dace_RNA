#!/bin/bash
#SBATCH --time=18:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=20000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL

module load nixpkgs/16.09
module load gcc/7.3.0
module load openmpi/3.1.2
module load salmon/1.1.0

salmon index -t /home/biomatt/scratch/dace/trinity_outputs_skylake3/redside_dace_transcriptome_Trinity.fasta \
-i /home/biomatt/scratch/dace/redside_dace_precorset_index2
