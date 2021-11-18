#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --job-name=job2
#SBATCH --time=00:22:00
#SBATCH --mail-type=none
# Queue (Partition):
#SBATCH --partition=express
### 
### for OpenMP:
###SBATCH --cpus-per-task=10

module load intel
module load mkl
module load impi

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INTEL_HOME/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel6
cd $PWD

EXE=/u/alaa/band/FHIaims/build/aims.190906.scalapack.mpi.x 
JNAME=control

srun -n 126 $EXE < $JNAME'.in' > $JNAME'.out'

exit
