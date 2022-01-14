#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J test_slurm
# Queue (Partition):
#SBATCH --partition=general
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#
#SBATCH --mail-type=none
#SBATCH --mail-user=<userid>@rzg.mpg.de
#
# Wall clock limit:
#SBATCH --time=14:00:00

module purge
module load mkl
module load intel
module load anaconda/2/2019.03
module load impi
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INTEL_HOME/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel6

date=$(date +"%s")
echo "the date is: $date"


AIMSPATH=/u/alaa/do/FHIaims/build
EXE=aims.x


export aimsbin=${AIMSPATH}/${EXE}

# Run the program:
python get_vibrations.py  MoS2 -r $aimsbin -s srun  -d 0.002  1  -w  > output

