from ase import Atoms
from ase.io import read, write
import sys
import numpy as np

__all__ = ['main']

### conversions ###
amu2me=1822.8885
ang2bohr=1.8897261
evtohartree=0.036749322
invcmtoau=4.55633e-06
kbttohartree=3.1668116e-06
atm2ap=3.4439668e-09

def inertia(system):
    ''' Calculation of moment of intertia of a molecule. ASE returns  in amu*angstrom**2. '''
    Is=system.get_moments_of_inertia()
    Is=Is*amu2me*ang2bohr*ang2bohr
    return Is

def translationat(system, kt):
    '''ase returns masses in amu'''
    mass=system.get_masses()
    totalmass=np.sum(mass)*amu2me
    p0=1*atm2ap # 1 atm is the reference pressure
    natoms=len(system)
    trans=-0.5*kt*np.log((totalmass/(2.*np.pi))**3.*(kt/p0)**5)
    return trans/natoms

def rotationsat(system, kt):
    minert=inertia(system)
    sigma=1.0
    natoms=len(system)
    t1=-0.5*kt*np.log(np.pi/sigma**2)
    t2=-0.5*kt*np.log((2.*kt)**3.*minert[0]*minert[1]*minert[2]) # may be missing pi here CHECK
    rot=t1+t2
    return rot/natoms

def vibrations(w, kt, cut):
    zpe = np.sum(np.array([i/2.0 for i in w[cut:] ]))
    therm=np.sum(np.array( [kt*np.log(1-np.exp(-i/kt)) for i in w[cut:] ] ))
    vib=zpe+therm
    return vib

def get_frequencies(fname):
    """ reads frequencies in cm-1 and returns them in atomic units """
    ftemp=open(fname)
    fw=np.array([float(line.split()[0]) for line in ftemp])*invcmtoau
    ftemp.close()
    return fw    

def get_energies(fname):
    """ file has a very specific format, energies are in eV, returns them in au """
    ftemp=open(fname)
    energies = [float(line.split()[1]) for line in ftemp]
    ftemp.close()
    # here should check length of file
    if (len(energies) != 3):
        print("too many or too few energies!") 
    return energies[0]*evtohartree, energies[1]*evtohartree, energies[2]*evtohartree

def main(geometry, w0, wd, wref, energies, dn):
    """ geometry: fhi-aims format
        w0: single column file with vib freq in cm-1 for pristine sys
        wd: single column file with vib freq in cm-1 for defect sys
        wref: single column file with vib freq in cm-1 for reference
        energies: 2 column file with energy of ref, energy pristine, energy defect in eV
        dn: integer variation of number of atoms"""
    molc=read(geometry,format='aims')
    temperatures = np.array([i for i in np.arange(0.00001, 800, 50)]) # this is in Kelvin
    ktarray=temperatures*kbttohartree 
    pressures = np.array([1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,1., 10.]) # this is in atm
    parray = pressures*atm2ap
    p0=1*atm2ap
    deltan=float(dn)
    fwref=get_frequencies(wref)
    fw0=get_frequencies(w0)
    fwd=get_frequencies(wd) 
    eref, e0, ed = get_energies(energies)
    eref=eref/len(molc)
    deltae=ed-e0
    for kt in ktarray:
        tcont=translationat(molc, kt)
        rcont=rotationsat(molc, kt)
        vcontref=vibrations(fwref, kt, 6)/len(molc)
        vcontprist=vibrations(fw0, kt, 3)
        vcontdefect=vibrations(fwd, kt, 3)
        deltaf=vcontdefect-vcontprist
        mu0=eref+tcont+rcont+vcontref
        total=(deltae+deltaf-deltan*mu0)
#        print(kt/kbttohartree, tcont, rcont, vcontref, (deltae+deltaf+1.0*mu0)/evtohartree)
#        print(kt/kbttohartree, (deltae+deltaf+1.0*mu0)/evtohartree)
        for p in parray:
            mupt=kt*np.log(p/p0)/len(molc)
            totalpt=total-deltan*mupt
            print(p/atm2ap, kt/kbttohartree, totalpt/evtohartree)
        print(" ") 

if __name__ == '__main__':
   main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])

