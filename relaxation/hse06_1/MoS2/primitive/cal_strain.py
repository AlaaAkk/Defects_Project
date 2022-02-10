from ase import Atoms
import numpy
import sys
from ase.build import make_supercell
from ase.build import mx2
from ase.io import read , write
import numpy as np
from numpy.linalg import inv
from ase.build import fcc111, root_surface

latvec = []
F = []
F1 = []

geo =fcc111('Au', (1, 1,4), a=4.14509402)
#geo = read("geometry.in", 0, "aims")
slab = make_supercell(geo,numpy.diag([12,12,4]))
slab.center(vacuum=50,axis=(2))
#slab.center()
slab.write('geometry_supercell.in',format='aims')
for line in open("geometry_supercell.in"):
    line = line.split("#")[0]
    words = line.split()
    if len(words) == 0:
        continue
    if words[0] == "lattice_vector":
        if len(words) != 4:
            raise Exception("geometry.in: Syntax error in line '"+line+"'")
        latvec += [ np.array(list(map(float,words[1:4]))) ]

if len(latvec) != 3:
    raise Exception("geometry.in: Must contain exactly 3 lattice vectors")

F = np.asarray(latvec)

print("Lattice vectors of film supercell :")
for i in range(3):
    print(F[i,:])
print()

F1=F[0:2,:]


latvec = []
S = []
S1 = []
geo = read("mos2.in", 0, "aims")
geo.center(vacuum=50,axis=(2))
slab = make_supercell(geo,numpy.diag([11,11,1]))
#slab.center(vacuum=50,axis=(2))
slab.write('geometry_mos2.in',format='aims')
for line in open("geometry_mos2.in"):
    line = line.split("#")[0]
    words = line.split()
    if len(words) == 0:
        continue
    if words[0] == "lattice_vector":
        if len(words) != 4:
            raise Exception("geometry.in: Syntax error in line '"+line+"'")
        latvec += [ np.array(list(map(float,words[1:4]))) ]

if len(latvec) != 3:
    raise Exception("geometry.in: Must contain exactly 3 lattice vectors")

S = np.asarray(latvec)

print("Lattice vectors of Substratei Supercell:")
for i in range(3):
    print(S[i,:])
print()
S1=S[0:2,:]
#n_array=np.matmul(F,inv(F))
print('F', F1)
print('S', S1)
n_array=np.matmul(S1, np.linalg.pinv(F1))
det = np.linalg.det(n_array)
print('det of T is', det)
