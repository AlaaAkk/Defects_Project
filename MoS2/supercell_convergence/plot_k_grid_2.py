import os #like in the terminal
import  pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#k-grid convergence

d=pd.read_csv( 'data_2' ,sep='\s+', header=0) #os.path other than writing the full path
#header means first line
#['#k', 'HF', 'pbe', 'lda']
plt.plot(d['x'], d['y']-min(d['y']), 'ro-')
plt.legend(fontsize=10)
plt.title('Convergence of Formation energy of Mo mono vacancy', fontsize=16)
plt.xlabel('supercell size', fontsize=11)
plt.ylabel('Delta Formation energy [eV]', fontsize=11)
new_names=['3x3x1','4x4x1','5x5x1','6x6x1','8x8x1']
plt.xticks(range(1, 6, 1), new_names, fontsize=16)
plt.yticks(fontsize=10)
plt.tight_layout() #better
plt.savefig(os.path.join( 'lat_conv.png'), dpi=400)
plt.show()
