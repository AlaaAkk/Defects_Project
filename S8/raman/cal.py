import  math
from math import *
import  pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from numpy import array, reshape, zeros, append, arange
import math
from math import log, e, pi
import numpy


d0=pd.read_csv('freq', sep='\s+',header=0)
T=500  #k
p0=1013250  #
k=1.380649e-16
h=6.62607015e-27
hbar=1.054571817e-27
sigma=8
pi=3.14159265359
m=3.1941378e-22
IA=IB=1.314051643394595e-37
IC=2.42660958899724e-37
R=8.31446261815324

convert=29979245800.0*2*pi #cm^-1 to Hz
w0=d0['x']*convert
A=np.log((((2*pi*m)**(3/2))*((k*T)**(5/2)))/(p0*(h**3)))
B=np.log(sqrt(pi)/sigma)
C=np.log(((8*pi*k*T/(h**2))**(3/2))*IA*sqrt(IC))
temp=numpy.array([(np.log(1-math.exp(-(hbar*i)/(k*T)))) for i in w0])
F=numpy.sum(temp)
print('F',F)
print('A', A)
print('B', B)
print('C', C)
final=A+B+C
print(R*F)
print(R*final)
print(R*(final-F))
