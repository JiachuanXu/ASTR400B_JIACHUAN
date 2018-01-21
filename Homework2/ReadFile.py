# ReadFile.py
# To read the time, particle number and particle info
# (including position, velocity, mass) from simulation data
# By Jiachuan Xu, Jan 21, 2018

import numpy as np
import astropy.units as u

# The Read function reads and returns
# time, total particle number of a snapshot,
# position, velocity and mass of each particle
# Input: <filename> #include path, e.g. ../MW_000.txt
# Output: (time [10Myr], particle bumber [int], data[label][:] [float])
# label is among: type, m, x, y, z, vx, vy, vz
def Read(filename):

	# Read header
	file = open(filename,'r')
	line1 = file.readline()
	label, value = line1.split()
	time = float(value)*10.0*u.Myr # time in units of 10Myr
	line2 = file.readline()
	label, value = line2.split()
	tot_num = int(value) # total particles number for this snapshot
	file.close()

	# Read info for each particle 
	data = np.genfromtxt(filename,dtype = None,names = True,skip_header = 3)

	# Return info readed
	return time, tot_num, data