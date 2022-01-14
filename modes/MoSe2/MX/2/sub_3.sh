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
#SBATCH --time=24:00:00

# Run the program:
module load intel/19.1.3
module load mkl/2020.4
module load impi/2019.9
module load anaconda/3/2020.02
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INTEL_HOME/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel6


AIMSPATH=/u/alaa/band/FHIaims/build
EXE=aims_latest.x
export aimsbin=${AIMSPATH}/${EXE}
python get_vibrations_3.py  WSe2 1  -d 0.0025 -r $aimsbin -s srun   > get_output



