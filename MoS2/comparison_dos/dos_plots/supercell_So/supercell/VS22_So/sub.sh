#!/bin/bash -l

# Standard output and error:
#SBATCH -o ./out.%j
#SBATCH -e ./err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J Alaa
# Queue (Partition):
#SBATCH --partition=general
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#
#SBATCH --mail-type=all
#SBATCH --mail-user=akkosh@fhi-berlin.mpg.de
#
# Wall clock limit:
#SBATCH --time=00:40:00

#Run the program:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64_lin/
module load anaconda/3/5.1
module load qe/6.4
python3 unfold.py --post
