# Calculation of phonon vibrational spectrum

import  pandas as pd
import matplotlib
import matplotlib as mpl
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import numpy as np
from numpy import array, reshape, zeros, append, arange
import math
from math import log, e, pi
import numpy
import sys
from sys import argv

pi=numpy.pi
convert=29979245800.0*2*pi # cm^-1 to Hz



d0=pd.read_csv('w0', sep='\s+',header=0)
d1=pd.read_csv('addS', sep='\s+',header=0)
d2=pd.read_csv('VS', sep='\s+',header=0)
d3=pd.read_csv('VS2', sep='\s+',header=0)
d4=pd.read_csv('VS22', sep='\s+',header=0)
d5=pd.read_csv('VW', sep='\s+',header=0)



w0=d0['w0']
w1=d1['addS']
w2=d2['VS']
w3=d3['VS2']
w4=d4['VS22']
w5=d5['VW']  #THZ to Hz



def _main(N,sigma):

    def VDOS(omega):
       F=[]

       omega = numpy.array(omega)
       for w in arange(0,450,1):

         temp3=numpy.array([math.exp(-((w-i)**2)/(2*sigma**2)) for i in omega])
         temp3=temp3/N
         F.append(numpy.sum(temp3))

       return F
    F0=[]
    F1=[]
    F2=[]
    F3=[]
    F4=[]
    F5=[]
    F6=[]
    F7=[]
    F7=[]
    F8=[]
    w=arange(0,450,1)

    F0=VDOS(w0)
    F1=VDOS(w1)
    F2=VDOS(w2)
    F3=VDOS(w3)
    F4=VDOS(w4)
    F5=VDOS(w5)
    #print(F0)
    fig, axs = plt.subplots(6, 1, sharex=True)
    axs[0].plot(w ,F0,'r', label='pristine')
    axs[0].legend()
    axs[1].plot(w ,F1,'b', label='adatom S')
    axs[1].legend()
    axs[2].plot(w ,F2,'g', label='VS')
    axs[2].legend()
    axs[3].plot(w ,F3,'k', label='VS2')
    plt.ylabel('VDOS', size=16)
    axs[3].legend()
    axs[4].plot(w ,F4,'y', label='VS22')
    axs[4].legend()
    axs[5].plot(w ,F5, 'c', label='VMo')
    axs[5].legend()

    plt.xlabel('w (cm-1)', size=16)
    plt.show()



if __name__ == "__main__":
      _main(int(sys.argv[1]),float(sys.argv[2]))
