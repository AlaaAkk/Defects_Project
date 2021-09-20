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

from phonopy import load

from vibes.brillouin import get_bands_and_labels
from vibes.structure.convert import to_Atoms

phonon = load("phonopy.yaml")

primitive = to_Atoms(phonon.primitive)

# this is the default collection of bandpaths, can be anything else
paths = primitive.cell.get_bravais_lattice().special_path.split(",")
bands, labels = get_bands_and_labels(primitive, paths, latex=False)

# compute the new bandstructure and save to file
phonon.run_band_structure(bands, labels=labels)
phonon.write_yaml_band_structure(filename="custom_band.yaml")


