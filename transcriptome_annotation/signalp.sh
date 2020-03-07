#!/bin/bash
#SBATCH --account=def-kmj477
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mem=50000M

module load nixpkgs/16.09 signalp-max/4.1f

perl $EBROOTSIGNALP/signalp -f short -n /home/biomatt/scratch/dace/signalp.out /home/biomatt/scratch/dace/transdecoder_out/longest_orfs.pep