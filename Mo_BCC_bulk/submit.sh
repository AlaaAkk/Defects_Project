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
#SBATCH --time=02:30:00

module load anaconda/3/5.1
source activate vibes


vibes run phonopy | tee log.phonopy
vibes output phonopy phonopy/trajectory.son --full

vibes utils fc remap
vibes utils fc frequencies
