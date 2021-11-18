# This is a script written by S. K. Wallace for plotting the total DOS and species DOS from an FHI-aims output for up to 4 species
#
### Basic usage
# To run the script use:
# python plot_DOS.py species1 species2 species3
# Where species1, species2, ... should be the element symbol, e.g. Cu
#
#
### Additional options
#
# Setting range for plot axes:
#
# By default the script will set the x-axis limits as the minimum and maximum energies of the total DOS for all plots so that they are aligned, however the user can also set custom limits to plot a limited energy region of the DOS (search for 'custom limits' in the script and comment out as you see fit)
# NOTE: currently the user also needs to uncomment or comment setting ylim at end of script for the plot, may try to automate this later
#
# For semiconductors:
#
# If the user also wishes to shift the zero of the energy of the plot so that it is coincident with the VBM, the user must search the aims output file for energy eigenvalue for the VBM and then replace 'VBM = 0' in this script with that value.
# The user can also add the band gap from the aims output file to mark it on the plot, by replacing 'Eg = 0' with the final value for the band gap found by searching for 'band' in the aims output file.
#
#############################################################################################################
#############################################################################################################
#############################################################################################################

import sys
import numpy as np
import math
import matplotlib.pyplot as plt

########### Reading in total DOS file, setting xlimits for all subplots and plotting total DOS ###############

# Reading in total DOS data from FHI-aims output
total_dos_file = "KS_DOS_total_raw_tetrahedron.dat"
energy_tot, dos_tot = np.loadtxt(total_dos_file, unpack=True)
#print(energy_tot)

#############################################################################################################
### USER INPUTS
#############################################################################################################
# Setting plot limits so all plots are set to same energy range as total dos
xmin=min(energy_tot)
xmax=max(energy_tot)
# Or custom limits
#xmin=-9
#xmax=9
ymin=0
ymax=1.5

species=['S0035','S0036','S0037','S0038','S0052','S0053']
# Option to shift zero of energy to VBM and mark on the band gap (by default these will be set to zero for no shifting)
# User inputs to shift zero of energy
#VBM = 0
#Eg=0
VBM = -5.73443138
Eg = -5.44599899
#Eg =-5.41264238

#############################################################################################################
#############################################################################################################

# Assigning number of subplots based on number of user-inputted arguments
#plt.subplot(, 1, 1)
#plt.subplot(len(species)+1, 1, 1)
#plt.fill_between(energy_tot, dos_tot, color='black')
#plt.xlim(xmin,xmax)
#plt.ylim(ymin,100*ymax)
#plt.ylabel("Total DOS", size=6)

if (VBM !=0) or (Eg !=0):
  # Adding lines and labels for VBM and CBM
  text_y_position = 2*ymax
  plt.axvline(x=VBM, color='black', linestyle='--')
  plt.axvline(x=Eg, color='black', linestyle='--')
  plt.text(VBM-0.3, text_y_position-7, 'VBM', fontsize=3)
  plt.text(Eg+0.3, text_y_position-7, 'CBM', fontsize=3)
  plt.title('VS Supercell',size=6)
#
#############################################################################################################
########## Starting to plot species proj dos requested by user ##############################################

for i in range(len(species)):
#  print i

  # Reading in arguments and assigning as species
  # Define data files from species requested by user
  #species=str(sys.argv[i])
  data_file = "atom_proj_dos_tetrahedron_"+str(species[i])+"_raw.dat"
  #data_file2 = "atom_proj_dos_Mo0034_raw.dat"
  #data_file3 = "atom_proj_dos_Mo0037_raw.dat"
  #data_file4 = "atom_proj_dos_Mo0049_raw.dat"
  #data_file5 = "atom_proj_dos_Mo0051_raw.dat"
  #data_file6 = "atom_proj_dos_Mo0066_raw.dat"

  # Determining how many orbitals for species_proj_dos are to be plotted and dividing up data file ready for plotting
  # Opening input file and storing all data to a 2D array
  with open(data_file) as data: # open file for reading
      all_data = [] # initialise empty list
      for line in data:
          all_data.append(line.strip().split())  # storing element and coordinates into 2D array of strings (or list appended onto another list?)
  #        print all_data
  data.close()
  #if (all_data[3][-1]=="4"):
   #print "4 -Pb"
  energy, total_species_dos, l0, l1, l2, l3, l4 = np.loadtxt(data_file, unpack=True)
  #  print(energy)
  #else:
   # print "Error in data file read in!"

  plt.subplot(len(species)+1, 1, 1+i)
  #plt.subplot(2, 1, 2)
  plt.fill_between(energy, total_species_dos, color='gray', alpha=0.09)
  #if (all_data[3][-1]=="4"):
  plt.plot(energy, l0, label='s', color='purple', lw=1 )
  plt.plot(energy, l1, label='p', color='blue', lw=1)
  plt.plot(energy, l2, label='d', color='green', lw=1)
  plt.plot(energy, l3, label='f', color='red', lw=1)
  plt.plot(energy, l4, label='g', color='orange', lw=1)
  #if (VBM !=0) or (Eg !=0):
  # Adding lines to mark VBM and CBM
  plt.axvline(x=VBM, color='black', linestyle='--')
  plt.axvline(x=Eg, color='black', linestyle='--')

#############################################################################################################
### USER INPUT: option to scale some pDOS axes differently
#############################################################################################################
#  if (str(species)=="As"):
#    ymax=6
#  else:
#    ymax=12
#############################################################################################################

  plt.yticks([])
  plt.xticks([])
  plt.xlim(xmin,xmax)
  plt.ylim(ymin,ymax)
  #plt.ylabel(' DOS', size=6)
  plt.ylabel(str(species[i]), size=6)

plt.xlabel('Energy (eV)', size=6) # Just for final plot
plt.subplots_adjust(bottom=0.0002)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.4), fancybox=True, shadow=True, ncol=6,fontsize='x-small')
plt.savefig('ldos_S.pdf')

plt.show()
