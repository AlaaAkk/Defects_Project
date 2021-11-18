import os #like in the terminal
import  pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

SMALL_SIZE = 6
MEDIUM_SIZE = 8
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
#k-grid convergence

d1=pd.read_csv( 'Temperature' ,sep='\s+', header=0) #os.path other than writing the full path
d2=pd.read_csv( 'Pressure' ,sep='\s+', header=0) #os.path other than writing the full path
fig, axs = plt.subplots(2)

axs[0].plot(d1['x'], d1['y'], 'r-', label='Constant pressure 1 atm')
axs[0].legend(fontsize=6)
axs[0].set_title('Formation energy as function of Tempreture', fontsize=9)
axs[1].plot(d2['x'], d2['y'], 'b-', label='Constant Tempreture 550K')
axs[1].legend(fontsize=6)
axs[1].set_title('Formation energy as function of Pressure', fontsize=9)
plt.show()

