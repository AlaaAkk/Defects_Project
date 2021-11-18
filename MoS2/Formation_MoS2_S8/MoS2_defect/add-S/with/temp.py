import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
import numpy as np
import math

E_defect=-3335811.405118595
E_pristine=-3324932.942491689
mu_S=-10879
w=8.106968900897939e+13 # s-1
hbar=6.58211956e-16  #ev.s
#deltaH=mu_MoS2-mu_Mo

T_i=10
T_f=1000

k=8.617333262145e-5 #ev,k-1
Ef=[]
jj=0
for i in arange(T_i,T_f,100):
    T=i
    temp=1-math.exp((hbar*w)/(k*T))
   # print(np.log(-temp))
    F=k*T*np.log(1-math.exp((-hbar*w)/(k*T)))
   # print(F)
    Ef=E_defect-E_pristine-(1*mu_S)+0.5*(hbar*w)+F
 #   print(str(T) +' ' + str(Ef))
    plt.plot(T, Ef, 'bo')   
    plt.xlabel('Temperature')
    plt.ylabel('Formation Energy')
    jj=jj+1
#plt.yticks(arange(-5000, 2, 200.0))    
plt.show()

p_i=1e-4
p_f=1
p0=1
T=550
for i in arange(p_i,p_f,0.001):
    p=i
    temp=1-math.exp((hbar*w)/(k*T))
   # print(np.log(-temp))
    F=k*T*np.log(1-math.exp((-hbar*w)/(k*T)))
   # print(F)
    Ef=E_defect-E_pristine-(1*mu_S)+(0.5*k*T*np.log(p/p0))+0.5*(hbar*w)+F
    print(str(p) +' ' + str(Ef))
    plt.plot(p, Ef, 'ro')   
    plt.xlabel('pressure')
    plt.ylabel('Formation Energy')
    jj=jj+1
#plt.yticks(arange(-5000, 2, 200.0))    
plt.show()
