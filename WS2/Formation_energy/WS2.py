#!/usr/bin/env python
# coding: utf-8

# In[1526]:


import  pandas as pd
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import numpy as np
from numpy import array, reshape, zeros, append, arange
import math
from math import log, e, pi
import numpy
from sys import argv


#  #### Total Energies in eV

# In[1527]:


E0=-11861477.311070437 # pristine
E1=-11872355.158936081# addon S
E2=-11850595.236526014 # mono S vacancy
E3=-11839713.678986968 # di S vacancy up&down
E4=-11839713.145423930 # di S vacancy neighboring
E5=-11408771.970668823 # mono W vacancy
E_WS2=-474458.138971034 # primtitive
ES8=-87031.629750345    #8 atoms in unitcell of S8 ring
EW=-905397.333160509 # W BCC


# #### Constants

# In[1528]:


pi=numpy.pi
convert=29979245800.0*2*pi # cm^-1 to Hz


# #### Frequencies

# In[1529]:


d0=pd.read_csv('w0', sep='\s+',header=0)
d1=pd.read_csv('addS', sep='\s+',header=0)
d2=pd.read_csv('VS', sep='\s+',header=0)
d3=pd.read_csv('VS2', sep='\s+',header=0)
d4=pd.read_csv('VS22', sep='\s+',header=0)
d5=pd.read_csv('VW', sep='\s+',header=0)
dS=pd.read_csv('w', sep='\s+',header=0)
dW=pd.read_csv('W_BCC', sep='\s+',header=0)

wW=dW['x']*convert
wS8=dS['x']*convert
w0=d0['w0']*convert
w1=d1['addS']*convert
w2=d2['VS']*convert
w3=d3['VS2']*convert
w4=d4['VS22']*convert
w5=d5['VW']*convert


# In[1530]:


p=1914.57 # 1.89e-3 atm (atm to cgs * 1013250))
p0=1013250  # atm to g/(cm s^2)
kk=1.380649e-16 # erg/k (cm^2.g/ks^2)
k=8.617333262145e-05 # ev/k
h=6.62607015e-27  # erg.s
hb=6.582119569e-16 # eV.s
hbar=1.054571817e-27 # erg.s
sigma=8
m=4.258952992e-22 # 32.06*8 in amu changed to g
IA=1.314051643394595e-37 # g.cm^2
IB=1.314146807283309e-37
IC=2.42660958899724e-37


# In[1531]:


A=[]
B=[]
C=[]
D=[]
E=[]
mu_0=[]
I=np.sqrt(IA)*np.sqrt(IB)*np.sqrt(IC)
for T in range(100,2300,100):
#    if(T!=0):
        A = np.append(np.log((((2*pi*m)**(3/2))*((kk*T)**(5/2)))/(p0*(h**3))),A)
        B=np.append(np.log(np.sqrt(pi)/sigma)+ np.log((((8*pi*kk*T)/(h**2))**(3/2))*I),B)
        temp=numpy.array([(np.log(1-math.exp(-(hbar*i)/(kk*T)))) for i in wS8])
        C=np.append(np.sum(temp),C)
        E= np.append(k*T*np.log(p/p0),E)
mu_0=-k*T*(A+B-C)


# In[1532]:


temp2=numpy.array([(hb*i)/(2) for i in wS8])
D=np.sum(temp2)


# In[1533]:


T=arange(100,2300,100)
#mu_0=-k*T*(A+B-C)
plt.plot(T,mu_0, 'red', label='For S8 ring')
plt.xlabel('T in K', fontsize=12)
plt.ylabel(r'$\mu_0$ [eV]', fontsize=12)
plt.legend(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig('mu0_T.pdf')
plt.show()


# In[1534]:


T=arange(100,2300,100)
mu_S8=mu_0 + E + D + ES8
plt.plot(T,mu_S8, 'red', label='For S8 ring')
plt.xlabel('T in K', fontsize=12)
plt.ylabel(r'$\mu_{S8}$ eV', fontsize=12)
plt.legend(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()


# In[1535]:


mu_S=mu_S8/8
plt.plot(T,mu_S, 'red', label='For S')
plt.ylabel(r'$\mu_S$ eV', fontsize=12)
plt.legend(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig('muS_T.png',dpi=400)


