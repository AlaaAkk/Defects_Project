import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-3303169.508929158
E_pristine=-3324932.942491689
mu_MoS2=-3324932.942491689/25
mu_Mo=-111237
E_Mo=-1392514.497532745
mu_S=-10879
deltaH=mu_MoS2-mu_Mo
mu_i=deltaH/2

mu_f=mu_S
print(mu_i)
#mu=[]
Ef=[]
jj=0
for i in arange(mu_i,mu_f,0.05):
    mu=i
    Ef=E_defect-E_pristine+(2*mu)
    #plt.plot(mu, Ef, 'bo')   
    print(str(mu) +' ' + str(Ef))
   # print(Ef)
   # print(mu)
    jj=jj+1
#plt.yticks(arange(-5000, 2, 200.0))    
plt.show()
