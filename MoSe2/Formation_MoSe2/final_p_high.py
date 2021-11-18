#!/usr/bin/env python
# coding: utf-8

# # MoSe2

# In[174]:


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
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'notebook')
#get_ipython().run_line_magic('matplotlib', 'widget')


# # Formation energy as  chemical potential vary
#

# In[175]:


pi=numpy.pi
convert=29979245800.0*2*pi # cm^-1 to Hz


#  #### Total Energies in eV

# In[176]:


E0=-6108318.476328624 # pristine
E1=-6174864.306138340 # addon S
E2=-6041768.483695555 # mono S vacancy
E3=-5975218.839354995 # di S vacancy up&down
E4=-5975218.551673065 # di S vacancy neighboring
E5=-5997073.704 # mono Mo vacancy
E_MoSe2=-244332.785924554 # primtitive
ESe8=-532375.657399940   #8 atoms in unitcell
EMo=-222473.348568306  # Total energy in eV


# In[177]:


d0=pd.read_csv('w0', sep='\s+',header=0)
d1=pd.read_csv('addS', sep='\s+',header=0)
d2=pd.read_csv('VS', sep='\s+',header=0)
d3=pd.read_csv('VS2', sep='\s+',header=0)
d4=pd.read_csv('VS22', sep='\s+',header=0)
d5=pd.read_csv('VMo', sep='\s+',header=0)
dW=pd.read_csv('Mo_BCC', sep='\s+',header=0)
dSe=pd.read_csv('w_Se8', sep='\s+',header=0)
wSe=dSe['x']*convert
wW=dW['x']*convert # kj/mol to eV
w0=d0['w0']*convert
w1=d1['addS']*convert
w2=d2['VS']*convert
w3=d3['VS2']*convert
w4=d4['VS22']*convert
w5=d5['VMo']*convert #THZ to Hz


# In[178]:


p01=1

p=80928.077 # 1.89e-3 atm (atm to cgs * 1013250))
p0=1013250  # atm to g/(cm s^2)
kk=1.380649e-16 # erg/k (cm^2.g/ks^2)
k=8.617333262145e-05 # ev/k
h=6.62607015e-27  # erg.s
hb=6.582119569e-16 # eV.s
hbar=1.054571817e-27 # erg.s
sigma=8
m=1.04907603472e-21
IA=4.149045888664045e-37
IB=4.149088585718857e-37
IC=7.612594467348306e-37


# In[179]:


pi=numpy.pi
convert=29979245800.0*2*pi # cm^-1 to Hz


# In[180]:


xr=np.array(np.arange(199.85,1326.85,10))

xr=1/xr

yr=(8.0886-(4989.5*xr))

pS=np.exp(yr)*0.00131579 #mm


# In[181]:


def free_energy(omega,T):
   F=[]

   omega = numpy.array(omega)
   #for p in arange(1.89e-3,1.89e+3,e-1):

   temp3=numpy.array([(hb*i/2 + k*T*np.log(1-math.exp(-(hb*i)/(k*T)))) for i in omega])
   F.append(numpy.sum(temp3))


   return F
def DeltaF(X,Y,T):
   deltaF=[]
   F1=free_energy(X,T)
   F2=free_energy(Y,T)
   zip_object = zip(F1, F2)
   for i, j in zip_object:
      deltaF.append(i-j)
   return deltaF


# In[183]:


fig = plt.figure()
ax = plt.axes(projection='3d')

for T in arange(199.85,1326.85,10):

    D=[]
    E=[]
    mu_0=[]
    mu_S8=[]
    mu_Se=[]
    addSe=[]
    VSe=[]
    VSe2=[]
    VSe22=[]
    I=np.sqrt(IA)*np.sqrt(IB)*np.sqrt(IC)
    A=np.log((((2*pi*m)**(3/2))*((kk*T)**(5/2)))/(p0*(h**3)))
    B=np.log(np.sqrt(pi)/sigma)+ np.log((((8*pi*kk*T)/(h**2))**(3/2))*I)
    temp=numpy.array([(np.log(1-math.exp(-(hbar*i)/(kk*T)))) for i in wSe])
    C=np.sum(temp)
    mu_0=-k*T*(A+B-C)
    for i in pS:

        E.append(k*T*np.log(i/p01))


    temp2=numpy.array([(hb*i)/(2) for i in wSe])
    D=np.sum(temp2)
#print('D',D)
    mu_S8=mu_0 + np.array(E) + D + ESe8
    mu_Se=np.array(mu_S8)/8
    #T=arange(199.85,1326.85,10)

    z1 = pS
    addSe = [E1-E0-a + DeltaF(w1,w0,T) for a in (mu_Se)]
    VSe =  [E2-E0+a + DeltaF(w2,w0,T) for a in mu_Se]
    VSe2 =  [E3-E0+2*a + DeltaF(w3,w0,T) for a in mu_Se]
    VSe22 = [E4-E0+2*a + DeltaF(w4,w0,T) for a in  mu_Se]

    a = np.empty(np.size(z1))
    a.fill(T)
    TT=a
 #  ax = plt.axes(projection='3d')
# Data for a three-dimensional line

    addSe=numpy.concatenate( addSe, axis=0 )
    VSe=numpy.concatenate( VSe, axis=0 )
    VSe2=numpy.concatenate( VSe2, axis=0 )
    VSe22=numpy.concatenate( VSe22, axis=0 )

    ax.plot3D(TT, z1, np.array(addSe),'r')
#    ax.plot3D(TT, z1, np.array(VSe),'b')
#    ax.plot3D(TT, z1, np.array(VSe2),'g')
#    ax.plot3D(TT, z1, np.array(VSe22),'k')
#
ax.set_xlabel(' Tempreture [k]')
ax.set_zlabel(' Formation Energy [eV]')
ax.set_ylabel(' Pressure [atm]')
#print('mu_S8',mu_S8)
plt.show()


# In[ ]:




