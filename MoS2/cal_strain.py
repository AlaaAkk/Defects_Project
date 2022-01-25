from ase import Atoms
import numpy
import sys
from ase.build import make_supercell
from ase.build import mx2
from ase.io import read , write
import numpy as np
from numpy.linalg import inv
from ase.build import fcc111, root_surface, add_adsorbate
import math
from sklearn.metrics import mean_absolute_percentage_error

latvec = []
F = []
F1 = []

geo =fcc111('Au', (1, 1,4), a=4.13002462)
geo.center(axis=(2))
#geo = read("geometry.in", 0, "aims")
slab = make_supercell(geo,numpy.diag([9,9,4]))
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
geo.center(axis=(2))
slab2 = make_supercell(geo,numpy.diag([8,8,1]))
#slab.write('geometry_MoS2.in',format='aims')
slab2.center(vacuum=50,axis=(2))
slab2.translate(10)
slab2.write('geometry_mos2.in',format='aims')
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

print("Lattice vectors of Substrate Supercell:")
for i in range(3):
    print(S[i,:])
print()
S1=S[0:2,:]
#n_array=np.matmul(F,inv(F))
F1 = np.array(F1)
#F1 = F1.reshape(2,3)
S1 = np.array(S1)
#S1 = S1.reshape(2,3)
print("Lattice vectors of Substrate Supercell without z:")
for i in range(2):
    print(S1[i,:])
print()
print("Lattice vectors of Film Supercell without z:")
#print(F11)
for i in range(2):
    print(F1[i,:])
print()

F_inv=np.linalg.pinv(F1)
#for i in range(3):
#    print(F_inv[i,:])
#print()
n_array=np.matmul(S1, np.linalg.pinv(F1))
numpy.around(n_array,10,n_array)
print(n_array)
det = np.linalg.det(n_array)
print('det of T is', det)
print('% of strain ', (1-det)*100)
#def mape(actual, pred):
#    actual, pred = np.array(actual), np.array(pred)
#    return np.mean(np.abs((actual - pred) / (1.0*(actual))))
#result=mape(S1[:,0:2],F1[:,0:2])
#print(result)
#
#
#r2=mean_absolute_percentage_error(S1[:,0:2],F1[:,0:2])
#print(r2)


#total=add_adsorbate(slab, slab2, 2.69, 'ontop',offset=(1.2,3.2))
#ase.io.write('geometry_Au.in', total, format='aims')
