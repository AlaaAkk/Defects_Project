import pandas as pd
import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt

d=pd.read_csv( "spectralWeights.dat" ,sep='\s+', header=0)
x= d['k']
y= d['Energies']
z= d['spectral_weights']

#pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
#matplotlib.rcParams.update(pgf_with_rc_fonts)
#rc('text', usetex = True)

SMALL_SIZE = 16
MEDIUM_SIZE = 18
BIGGER_SIZE = 24

#plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
#plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
#plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
#plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
#plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
#plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
#plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
#

plt.scatter(x, y, c=z, cmap='hot')
plt.xlabel('K points')
plt.ylabel('Energy')
cbar = plt.colorbar()
cbar.set_label('Spectral Weights', rotation=270, labelpad=20, y=0.5,)


plt.tight_layout()
plt.savefig('Unfolded_Band.png')
plt.show()

