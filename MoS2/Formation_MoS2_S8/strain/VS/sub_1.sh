#!/bin/bash -l
#Standard output and error:
#SBATCH -o ./out.%j
#SBATCH -e ./err.%j
# Initial working directory:
#SBATCH -D ./
# Job name
#SBATCH -J Alaa
#
#SBATCH --nodes=4            # Request 1 (or more) node(s)
#SBATCH --constraint="gpu"    #    providing GPUs.
#SBATCH --ntasks-per-node=72  # Launch 72 tasks per node
#SBATCH --gres=gpu:a100:4     # Request all 4 GPUs of each node
#SBATCH --nvmps               # Launch NVIDIA MPS to enable concurrent access to the GPUs from multiple processes efficiently
#
#SBATCH --mail-type=all
#SBATCH --mail-user=akkoush@fhi-berlin.mpg.de
#SBATCH --time=24:00:00
export OMP_NUM_THREADS=1
module purge
module load mkl 
module load intel/19.1.2 
module load impi/2019.8
module list  
module load anaconda

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INTEL_HOME/compilers_and_libraries/linux/lib/intel64:$MKL_HOME/lib/intel64/:$HOME/.local/lib
#FHI
AIMSPATH=/u/alaa/FHIaims/build
EXE=aims.x


export aimsbin=${AIMSPATH}/${EXE} 
python get_vibrations_1.py  WSe2 1  -d 0.0025 -r $aimsbin -s srun   > get_output
