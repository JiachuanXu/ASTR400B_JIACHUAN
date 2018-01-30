# GalaxyMass.py
# To give the statistic of the components' mass of given galaxy,
# total mass, baryon fraction.
# If more than one galaxy are given, they should form a group, then
# this script also return the total mass and total baryon fraction
# the galaxy group.
# Usage: $ python GalaxyMass.py <file 1> <file 2> <filen 3> ...
# The file name should include path
# Output: A table in pdf format
# By Jiachuan Xu, Jan 23, 2018
import numpy as np
import astropy.units as u
from astropy.table import QTable, Column
from astropy.io import ascii
import sys, string, math
from ReadFile import Read

# The ComponentMass function returns the total mass of specific
# component of galaxy
# Input: filename, particle type
# particle type: 1-Dark Matter, 2-Disk Stars, 3-Bulge Stars
# Output: total mass [M_sun]
def ComponentMass(filename, particle_type):
	skip, skip, data = Read(filename)
	index = np.where(data['type'] == particle_type)
	tot_mass = np.around(math.fsum(data['m'][index])*1e10, 3)*u.solMass
	return tot_mass

# Main program

# Confirm the arguments
file_num = len(sys.argv)
if file_num == 1:
	print("Usage: $ python GalaxyMass.py <file 1> <file 2> <filen 3> ...")
	exit()

# Generate the output table, set the labels, data types and units
# | name | halo mass | disk mass | bulge mass | total mass | baryon fraction |
gm_list = QTable(\
	names=('Name', r'$M_{Halo}$', r'$M_{Disk}$', r'$M_{Bulge}$',\
		r'$M_{Total}$', r'$f_{bar}$'),\
	meta={'name':'first table'}, \
	dtype=('|U16', np.float, np.float, np.float, np.float, np.float))
gm_list[r'$M_{Halo}$'] = []*u.solMass
gm_list[r'$M_{Disk}$'] = []*u.solMass
gm_list[r'$M_{Bulge}$'] = []*u.solMass
gm_list[r'$M_{Total}$'] = []*u.solMass
gm_list[r'$f_{bar}$'] = []*u.dimensionless_unscaled

# Initialize the total mass of components
group_halo, group_disk, group_bulge = 0.*u.solMass, 0.*u.solMass, 0.*u.solMass

# Loop through the galaxies
for file in range(1, file_num):
	# Extract the galaxy name
	galaxy_name = sys.argv[file].split('/')[-1]
	galaxy_name = galaxy_name.split('_')[0]
	# Return components' mass
	halo_mass = ComponentMass(sys.argv[file], 1.)
	disk_mass = ComponentMass(sys.argv[file], 2.)
	bulge_mass = ComponentMass(sys.argv[file], 3.)
	# Add to local group mass 
	group_halo += halo_mass
	group_disk += disk_mass
	group_bulge += bulge_mass
	# Calculate galaxy total mass and baryon fraction
	galaxy_mass = halo_mass + disk_mass + bulge_mass
	baryon_frac = np.around((disk_mass + bulge_mass)/galaxy_mass, 3)
	# Add a row to the table
	galaxy = [galaxy_name, halo_mass, disk_mass, bulge_mass, galaxy_mass, baryon_frac]
	gm_list.add_row(galaxy)

# Edit the row of local group: total mass and baryon fraction
group_mass = group_halo + group_disk + group_bulge
group_fbar = np.around((group_disk + group_bulge)/group_mass, 3)
galaxy_group = ["Galaxy Group", group_halo, group_disk, group_bulge, group_mass, group_fbar]
gm_list.add_row(galaxy_group)

# Convert the units into 1e12 M_sun 
gm_list[r'$M_{Halo}$'] = np.around(gm_list[r'$M_{Halo}$'].to(1e12*u.solMass), 3)
gm_list[r'$M_{Disk}$'] = np.around(gm_list[r'$M_{Disk}$'].to(1e12*u.solMass), 3)
gm_list[r'$M_{Bulge}$'] = np.around(gm_list[r'$M_{Bulge}$'].to(1e12*u.solMass), 3)
gm_list[r'$M_{Total}$'] = np.around(gm_list[r'$M_{Total}$'].to(1e12*u.solMass), 3)

# Print the table, and save in tex format
# NOTE: IF THE TERMINAT SAYS "TypeError: unhashable type: 'MaskedConstant'"
# WHILE RUNNING, PLEASE UPDATE YOUR NUMPY TO THE LATEST VERSION FROM
# https://github.com/numpy/numpy.git
print(gm_list)
ascii.write(gm_list, 'mass_break_down.tex', format='latex')