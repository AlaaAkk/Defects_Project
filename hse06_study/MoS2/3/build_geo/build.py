from ase import Atoms
import numpy
import sys
from ase.build import make_supercell
from ase.build import mx2
from ase.io import read , write
import numpy as np
from numpy.linalg import inv
from ase.build import fcc111, root_surface

geo = read("mos2.in", 0, "aims")
geo.center(vacuum=50,axis=(2))
slab = make_supercell(geo,numpy.diag([5,5,1]))
#slab.center(vacuum=50,axis=(2))
slab.write('geometry_mos2.in',format='aims')
