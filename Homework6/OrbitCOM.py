# OrbitCOM.py 
# To calculate com position and velocity and time through snap numbers
# Usage: $ python OrbitCOM.py
# By Jiachuan Xu, Feb, 20, 2018

# import modules
import numpy as np 
import astropy.units as u 
from ReadFile import Read 
from CenterOfMass import CenterOfMass

# define OrbitCOM function
# input: galaxy: galaxy name
# 		 start: start snapnumber
#		 end: end number 
#		 n: interval for output COM info
# output: a ascii list of COM info:
#	t(Gyr) x y z(kpc) vx vy vz(km/s)
def OrbitCOM(galaxy, start, end, n):
	# define output filename, initialize
	fileout = "Orbit_%s.txt" % galaxy
	rows = int((end-start)/n)+1
	Orbit = np.zeros((rows, 7))
	
	## DEFINE TOLERANCE AND VOLUMN DECREASE FACTOR
	delta = 0.5
	VolDec = 4.

	# loop through snapnumbers
	for i in np.arange(start, end+1, n):
		print("Looping at snap number = %d\n"%i)
		# index in array
		ct = int((i-start)/n)
		# construct filename
		filename = "%s/%s_%03d.txt"%("../MW_VLowRes", galaxy, i)
		# construct COM object using disk particles, get COM info
		center = CenterOfMass(filename, 2)
		rc = center.COM_P(delta, VolDec)
		vc = center.COM_V(30, rc)
		# write data into list
		Orbit[ct][0] = float(center.time/u.Gyr)
		Orbit[ct][1] = float(rc[0]/u.kpc)
		Orbit[ct][2] = float(rc[1]/u.kpc)
		Orbit[ct][3] = float(rc[2]/u.kpc)
		Orbit[ct][4] = float(vc[0]/(u.km/u.s))
		Orbit[ct][5] = float(vc[1]/(u.km/u.s))
		Orbit[ct][6] = float(vc[2]/(u.km/u.s))
	# print list
	np.savetxt(fileout, Orbit, header='#t x y z vx vy vz', \
		comments='# time in Gyr, position in kpc, velocity in km/s\n',\
		fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])

# Main program
# Generate list of MW, M31 and M33
OrbitCOM("MW", 0, 800, 5)
OrbitCOM("M31", 0, 800, 5)
OrbitCOM("M33", 0, 800, 5)


