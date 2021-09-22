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
#SBATCH --ntasks-per-node=40
#
#SBATCH --mail-type=all
#SBATCH --mail-user=akkosh@fhi-berlin.mpg.de
#
# Wall clock limit:
#SBATCH --time=00:30:00

#Run the program:
module load qe/6.4

export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INTEL_HOME/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel6

srun  dos.x < dos.in > dos.out
