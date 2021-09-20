#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J BKB
# Queue (Partition):
#SBATCH --partition=medium
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#
#SBATCH --mail-type=all
#SBATCH --mail-user=akkoush@fhi-berlin.mpg.de
#
# Wall clock limit:
#SBATCH --time=00:50:00
module load anaconda/3/5.1
conda create -n py37 python=3.7 numpy scipy mkl
conda activate py37


vibes run phonopy | tee log.phonopy
vibes output phonopy phonopy/trajectory.son --full

