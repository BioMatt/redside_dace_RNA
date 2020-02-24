#!/bin/bash
#SBATCH --account=def-kmj477
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --time=120:00:00
#SBATCH --mem=200000M

module load nixpkgs/16.09 intel/2018.3
module load corset/1.07

corset -D 99999999999 -g 1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3 \
-n CTmax.10,CTmax.1,CTmax.2,CTmax.3,CTmax.4,CTmax.5,CTmax.6,CTmax.7,CTmax.8,CTmax.9,\
Handle.10,Handle.1,Handle.2,Handle.3,Handle.4,Handle.5,Handle.6,Handle.7,Handle.8,Handle.9,\
Wild.10,Wild.1,Wild.2,Wild.3,Wild.4,Wild.5,Wild.6,Wild.7,Wild.8,Wild.9 \
-i salmon_eq_classes /home/biomatt/scratch/dace/salm_outputs_precorset/*/aux_info/eq_classes.txt \
-p /home/biomatt/scratch/dace/dace_corset2/dace \
-l 5 \
-x 1000
