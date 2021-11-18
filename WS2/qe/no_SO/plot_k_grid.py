import os #like in the terminal
import  pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#k-grid convergence

d1=pd.read_csv( 'aims.dat' ,sep='\s+', header=0) #os.path other than writing the full path
d2=pd.read_csv( 'qe.dat' ,sep='\s+', header=0) #os.path other than writing the full path
#header means first line
#['#k', 'HF', 'pbe', 'lda']
plt.plot(d1['x'], d1['y'],'b-', label='Fhi-aims')
plt.plot(d2['x']-0.7, d2['y'], 'r-', label='QE')
plt.xlim([-16, 0])
plt.legend(fontsize=10)
#plt.title('Convergence of Formation energy of Mo mono vacancy', fontsize=16)
plt.xlabel('Energies eV', fontsize=16)
plt.ylabel('DOS', fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout() #better
plt.savefig(os.path.join( 'no_So_comp.pdf'))
plt.show()
