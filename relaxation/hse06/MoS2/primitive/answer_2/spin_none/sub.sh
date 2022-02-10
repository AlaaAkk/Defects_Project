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
#SBATCH --mail-type=none
#SBATCH --mail-user=<userid>@rzg.mpg.de
#
# Wall clock limit:
#SBATCH --time=01:40:00

# Run the program:
module load intel/19.0.5
module load mkl/2019.5
module load impi/2019.5
module load anaconda/3/2020.02
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INTEL_HOME/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel6


AIMSPATH=/u/alaa/band/FHIaims/build
EXE=aims.x
export aimsbin=${AIMSPATH}/${EXE}
srun $aimsbin > out_2



