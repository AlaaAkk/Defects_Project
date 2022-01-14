#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
from sys import argv



# # Formation energy as  chemical potential vary
#

# In[2]:


pi=numpy.pi
convert=29979245800.0*2*pi # cm^-1 to Hz


#  #### Total Energies in eV
######################################################

d0=pd.read_csv('w0', sep='\s+',header=0)
d1=pd.read_csv('addS', sep='\s+',header=0)
d2=pd.read_csv('VS', sep='\s+',header=0)
d3=pd.read_csv('VS2', sep='\s+',header=0)
d4=pd.read_csv('VS22', sep='\s+',header=0)
d5=pd.read_csv('VW', sep='\s+',header=0)
dW=pd.read_csv('W_BCC', sep='\s+',header=0)
dSe=pd.read_csv('w_Se8', sep='\s+',header=0)
wSe=dSe['x']*convert
wW=dW['x']*convert # kj/mol to eV
w0=d0['w0']*convert
w1=d1['addS']*convert
w2=d2['VS']*convert
w3=d3['VS2']*convert
w4=d4['VS22']*convert
w5=d5['VW']*convert #THZ to


# # Formation energy as  chemical potential vary
#

def conc(r,g):
   rho=[]
   T=arange(199.85,1326.85,100)
   t=k*T
   for i,j in zip(r,t):
       temp=(g)*math.exp(-i/j)
       rho.append(temp)

   return rho
#  #### Total Energies in eV

# In[10]:


E0=-14644858.463124419 # pristine
E1=-14711404.113635940 # addon S
E2=-14578308.639903018 # mono S vacancy
E3=-14511759.369944071 # di S vacancy up&down
E4=-14511758.944663303 # di S vacancy neighboring
E5=-14192151.211700918 # -14192151.211700918 mono Mo vacancy
E_WSe2=-585793.321771426 # primtitive
#ESe8=-66548.05582020049*8   #8 atoms in unitcell
ESe8_2=-532375.657399940   #8 atoms in unitcell
#EMo=-111238.871574569*2 # Total energy in eV
EMo=-905397.333160509
ESe8=4*(E_WSe2-EMo/2)

#####################################################
xr=np.array(np.arange(199.85,1326.85,100))

xr=1/xr

yr=(8.0886-(4989.5*xr))

pS=np.exp(yr)*0.00131579 #mmhg to atm


# In[ ]:


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


# In[5]:
##################################################################
#  #### Total Energies in eV

# In[2]:


# ## $$ \mu=\mu_{0}+ kT \ln\frac{p}{p_{0}}+ E_{DFT}+ \sum_{i}\frac{\hbar \omega_{i}}{2}$$

# In[7]:
mpl.rcParams['legend.fontsize'] = 10

def free_energy(omega):
   F=[]

   omega = numpy.array(omega)
   for T in arange(199.85,1326.85,100):

     temp3=numpy.array([(hb*i/2 + k*T*np.log(1-math.exp(-(hb*i)/(k*T)))) for i in omega])
     F.append(numpy.sum(temp3))

   return F


# In[8]:


def DeltaF(X,Y):
   deltaF=[]
   F1=free_energy(X)
   F2=free_energy(Y)
   zip_object = zip(F1, F2)
   for i, j in zip_object:
       deltaF.append(i-j)
   return deltaF


# In[9]:


#fig = plt.figure()
#ax = plt.axes(projection='3d')

conc_addS_2=np.zeros((len(xr),len(xr)))
conc_VS_2=np.zeros((len(xr),len(xr)))
conc_VS2_2=np.zeros((len(xr),len(xr)))
conc_VS22_2=np.zeros((len(xr),len(xr)))
conc_VW_2=np.zeros((len(xr),len(xr)))


conc_addS_1=np.zeros((len(xr),len(xr)))
conc_VS_1=np.zeros((len(xr),len(xr)))
conc_VS2_1=np.zeros((len(xr),len(xr)))
conc_VS22_1=np.zeros((len(xr),len(xr)))
conc_VW_1=np.zeros((len(xr),len(xr)))




addS=np.zeros((len(xr),len(xr)))
VS=np.zeros((len(xr),len(xr)))
VW=np.zeros((len(xr),len(xr)))
VS2=np.zeros((len(xr),len(xr)))
VS22=np.zeros((len(xr),len(xr)))
addS_2=np.zeros((len(xr),len(xr)))
VW_2=np.zeros((len(xr),len(xr)))
VS_2=np.zeros((len(xr),len(xr)))
VS2_2=np.zeros((len(xr),len(xr)))
VS22_2=np.zeros((len(xr),len(xr)))
for pindx,p in enumerate(list(pS*p0)):
   D=[]
   E=[]
   mu_0=[]
   mu_S8=[]
   mu_S=[]
   mu_W=[]
   mu_S8_2=[]
   mu_S_2=[]
   mu_W_2=[]
   I=np.sqrt(IA)*np.sqrt(IB)*np.sqrt(IC)
   for Tindx,T in enumerate(list(arange(199.85,1326.85,100))):
        A=np.log((((2*pi*m)**(3/2))*((kk*T)**(5/2)))/(p0*(h**3)))
        B=np.log(np.sqrt(pi)/sigma)+ np.log((((8*pi*kk*T)/(h**2))**(3/2))*I)
        temp=numpy.array([(np.log(1-math.exp(-(hbar*i)/(kk*T)))) for i in wSe])
        C=np.sum(temp)
        E.append(k*T*np.log(p/p0))
        mu_0.append(-k*T*(A+B-C))


#print('E',E)
#print('mu_0', mu_0)
   temp2=numpy.array([(hb*i)/(2) for i in wSe])
   D=np.sum(temp2)
   mu_S8=np.array(mu_0) + np.array(E) + D + ESe8
   mu_S8_2=np.array(mu_0) + np.array(E) + D + ESe8_2
   mu_S=np.array(mu_S8)/8
   mu_W=E_WSe2-2*mu_S
   mu_S_2=np.array(mu_S8_2)/8
   mu_W_2=E_WSe2-2*mu_S_2
   addS[pindx,:] = [E1-E0-a + b for a, b in zip(mu_S, DeltaF(w1,w0))]
   VS [pindx,:] = [E2-E0+a + b for a, b in zip(mu_S, DeltaF(w2,w0))]
   VS2 [pindx,:] = [E3-E0+2*a + b for a, b in zip(mu_S, DeltaF(w3,w0))]
   VS22 [pindx,:] = [E4-E0+2*a + b for a, b in zip(mu_S, DeltaF(w4,w0))]
   VW [pindx,:] = [E5-E0+a + b for a, b in zip(mu_W, DeltaF(w5,w0))]
#VMo = [E5-E0+a + b for a, b in zip(mu_Mo, DeltaF(w5,w0))]

   addS_2[pindx,:] = [E1-E0-a + b for a, b in zip(mu_S_2, DeltaF(w1,w0))]
   VS_2[pindx,:] = [E2-E0+a + b for a, b in zip(mu_S_2, DeltaF(w2,w0))]
   VS2_2 [pindx,:] = [E3-E0+2*a + b for a, b in zip(mu_S_2, DeltaF(w3,w0))]
   VS22_2 [pindx,:] = [E4-E0+2*a + b for a, b in zip(mu_S_2, DeltaF(w4,w0))]
   VW_2 [pindx,:] = [E5-E0+a + b for a, b in zip(mu_W_2, DeltaF(w5,w0))]
   #mu_S8=np.array(mu_0) + np.array(E) + D + ESe8

   conc_addS_1[pindx,:] = conc(addS[pindx,:],1)
   conc_VS_1[pindx,:] = conc(VS[pindx,:],1)
   conc_VS2_1[pindx,:] = conc(VS2[pindx,:],6)
   conc_VS22_1[pindx,:] = conc(VS22[pindx,:],6)
   conc_VW_1[pindx,:] = conc(VW[pindx,:],6)

   conc_addS_2[pindx,:] = conc(addS_2[pindx,:],1)
   conc_VS_2[pindx,:] = conc(VS_2[pindx,:],1)
   conc_VS2_2[pindx,:] = conc(VS2_2[pindx,:],6)
   conc_VS22_2[pindx,:] = conc(VS22_2[pindx,:],6)
   conc_VW_2[pindx,:] = conc(VW_2[pindx,:],6)


  # ax.plot3D(T, z1,np.array(addS),'r')
  # ax.plot3D(T, z1,np.array(VS),'b')
  # ax.plot3D(T, z1,np.array(VS2),'g')
  # ax.plot3D(T, z1,np.array(VS22),'k')
##print(y22)
##plt.savefig('3D.png',dpi=400)
# set up a figure twice as wide as it is tall
fig = plt.figure(figsize=plt.figaspect(0.5))

#===============
#  First subplot
#===============
# set up the axes for the first plot
ax = fig.add_subplot(1, 2, 1, projection='3d')
T=arange(199.85,1326.85,100)
X, Y = np.meshgrid(pS, T)
#ax = plt.axes(projection='3d')
ax.plot_surface(Y, X, addS.T, rstride=1, cstride=1,color='red',shade=False)
ax.plot_surface(Y, X, VS.T, rstride=1, cstride=1,color='blue',shade=False)
ax.plot_surface(Y, X, VS2.T, rstride=1, cstride=1,color='green',shade=False)
ax.plot_surface(Y, X, VS22.T, rstride=1, cstride=1,color='orange',shade=False)
ax.plot_surface(Y, X, VW.T, rstride=1, cstride=1,color='k',shade=False)
ax.set_xlabel(' Tempreture [k]')
ax.set_zlabel(' Formation Energy [eV]')
ax.set_ylabel(' Pressure [atm]')
ax.set_title('WSe2 low Se environment')
#===============
# Second subplot
#===============
# set up the axes for the second plot
ax = fig.add_subplot(1, 2, 2, projection='3d')

# plot a 3D wireframe like in the example mplot3d/wire3d_demo

ax.plot_surface(Y, X, addS_2.T, rstride=1, cstride=1,color='red', label='addon Se',shade=False)
#ax.set_legend(fontsize=7)
ax.plot_surface(Y, X, VS_2.T, rstride=1, cstride=1,color='blue',shade=False)
ax.plot_surface(Y, X, VS2_2.T, rstride=1, cstride=1,color='green',shade=False)
ax.plot_surface(Y, X, VS22_2.T, rstride=1, cstride=1,color='orange',shade=False)
ax.plot_surface(Y, X, VW_2.T, rstride=1, cstride=1,color='k',shade=False)
ax.set_xlabel(' Tempreture [k]')
ax.set_zlabel(' Formation Energy [eV]')
ax.set_ylabel(' Pressure [atm]')
ax.set_title('WSe2 high Se environment')

b1 = plt.Rectangle((0, 0), 1, 1, fc="red")
b2 = plt.Rectangle((0, 0), 1, 1, fc="blue")
b3 = plt.Rectangle((0, 0), 1, 1, fc="green")
b4 = plt.Rectangle((0, 0), 1, 1, fc="orange")
b5 = plt.Rectangle((0, 0), 1, 1, fc="k")

ax.legend([b1, b2,b3,b4,b5], ['adatom', 'Mono Se', 'di up$\&$down','di neigh','Mono W'])
plt.show()




fig = plt.figure(figsize=plt.figaspect(0.5))
#===============
#  First subplot
#===============
# set up the axes for the first plot
ax = fig.add_subplot(1, 2, 1, projection='3d')
T=arange(473,1600,100)
X, Y = np.meshgrid(pS, T)
ax.plot_surface(Y, X, conc_addS_1.T, rstride=1, cstride=1,color='r',shade=False)
ax.plot_surface(Y, X, conc_VS_1.T, rstride=1, cstride=1,color='b',shade=False)
ax.plot_surface(Y, X, conc_VS2_1.T, rstride=1, cstride=1,color='g',shade=False)
ax.plot_surface(Y, X, conc_VS22_1.T, rstride=1, cstride=1,color='orange',shade=False)
ax.plot_surface(Y, X, conc_VW_1.T, rstride=1, cstride=1,color='k',shade=False)
ax.set_xlabel(' Tempreture [k]')
ax.set_zlabel(' Concentration')
ax.set_ylabel(' Pressure [atm]')
ax.set_title('WSe2 low Se environment')
b5 = plt.Rectangle((0, 0), 1, 1, fc="k")

ax.legend([b1, b2,b3,b4,b5], ['adatom', 'Mono Se', 'di up$\&$down','di neigh','Mono W'])
#===============
# Second subplot
#===============
# set up the axes for the second plot
ax = fig.add_subplot(1, 2, 2, projection='3d')
T=arange(473,1600,100)
X, Y = np.meshgrid(pS, T)
ax.plot_surface(Y, X, conc_addS_2.T, rstride=1, cstride=1,color='r',shade=False)
ax.plot_surface(Y, X, conc_VS_2.T, rstride=1, cstride=1,color='b',shade=False)
ax.plot_surface(Y, X, conc_VS2_2.T, rstride=1, cstride=1,color='g',shade=False)
ax.plot_surface(Y, X, conc_VS22_2.T, rstride=1, cstride=1,color='k',shade=False)
ax.plot_surface(Y, X, conc_VW_2.T, rstride=1, cstride=1,color='m',shade=False)
ax.set_xlabel(' Tempreture [k]')
ax.set_zlabel(' Concentration')
ax.set_ylabel(' Pressure [atm]')
ax.set_title('WSe2 high Se environment')
ax.legend([b1, b2,b3,b4,b5], ['adatom', 'Mono Se', 'di up$\&$down','di neigh','Mono W'])
plt.show()




# In[ ]:




