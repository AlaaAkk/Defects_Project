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
#SBATCH --time=00:40:00

# Run the program:

export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INTEL_HOME/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel6

date=$(date +"%s")
echo "the date is: $date"


AIMSPATH=/u/alaa/do/FHIaims/build
EXE=aims.x


export aimsbin=${AIMSPATH}/${EXE} 

srun $aimsbin > out_$date
