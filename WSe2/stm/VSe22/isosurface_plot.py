#!/usr/bin/env python

import os, sys, re
import pandas as pd
from optparse import OptionParser
import ase
import numpy as np
from ase.io.cube import read_cube_data, write_cube
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="write report to FILE", metavar="FILE")
parser.add_option("-d", "--delta", dest="delta", default = 0.00002,
                  help="delta for isovalue")
parser.add_option("-i", "--isovalue", dest="isovalue", default = -0.00005,
                  help="Isovalue for isosurface")
parser.add_option("-n", "--name", dest="name", default = 'Molecule',
                  help="Name of file")

(options, args) = parser.parse_args()

import matplotlib.pyplot as plt
from matplotlib import offsetbox
import matplotlib
from IPython.display import display, Math, Latex
from sklearn import (manifold, datasets, decomposition, ensemble,
                     discriminant_analysis, random_projection)
from optparse import OptionParser
import pandas as pd
import sys, os
import plotly.offline as pyo
import plotly.figure_factory as ff
import plotly.graph_objs as go
import plotly.io as pio
import matplotlib.gridspec as gridspec
import seaborn as sns
import plotly.express as px
import itertools
import shutil
import re
from matplotlib import rc
from scipy.interpolate import interp2d
from scipy.interpolate import griddata
import numpy.ma as ma
from numpy.random import uniform, seed

from matplotlib import rc
pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
matplotlib.rcParams.update(pgf_with_rc_fonts)
rc('text', usetex = True)

SMALL_SIZE = 16
MEDIUM_SIZE = 18
BIGGER_SIZE = 24

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

xmax = 14.69643837
ymax = 15.27298677

Voleum = pd.read_csv(os.path.join(os.getcwd(), options.filename), sep='\s+', header=0)
delta = np.abs(float(options.isovalue))*0.5
Isovalue = float(options.isovalue)

Isosurface = Voleum[(Isovalue-delta<Voleum['V'])&(Voleum['V']<Isovalue+delta)]
MinValue = min(Isosurface['Z'])
Isosurface['Z'] = Isosurface['Z'] - MinValue




plt.scatter(Isosurface['X'], Isosurface['Y'],
            c= Isosurface['Z'], cmap='hot')

plt.xlabel(r'Distance, \AA')
plt.ylabel(r'Distance, \AA')
plt.title(r'di Se Vacancy Neighboring')
cbar = plt.colorbar()
cbar.set_label('Height, \AA', rotation=270, labelpad=20, y=0.5,)


plt.tight_layout()
plt.savefig(os.path.join(os.getcwd(), 'STM_{}.png'.format(options.name)))
#plt.show()

plt.clf()

plt.scatter(Isosurface['X'], Isosurface['Y'],
            c= Isosurface['Z'], cmap='hot')

plt.xlabel(r'Distance, \AA')
plt.ylabel(r'Distance, \AA')
plt.title(r'add atom S')
cbar = plt.colorbar()
cbar.set_label('Height, \AA', rotation=270, labelpad=20, y=0.5,)
plt.vlines(x=np.linspace(0, 14.69643837, 6), ymin=0, ymax=15.27298677, colors='white', linestyles='dashed')
plt.hlines(y=np.linspace(0, 15.27298677, 6), xmin=0, xmax=14.69643837, colors='white', linestyles='dashed')
plt.tight_layout()
plt.savefig(os.path.join(os.getcwd(), 'STM_{}_lines.png'.format(options.name)))
#plt.show()


















#plt.clf()

#ax = plt.subplot(111)
#im=ax.imshow(np.array(Isosurface['Z']-min(Isosurface['Z'])).reshape(101, 101))
             ##)
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#ax.set_xlabel('r.l.u', fontsize=14)
#ax.set_ylabel('r.l.u', fontsize=14)
#plt.colorbar(im, cax=cax, label='Intensity (arb. units)')
#plt.show()


