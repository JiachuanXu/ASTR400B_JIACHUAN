import numpy as np
import astropy.units as u

def Read(filename):

	file = open(filename,'r')
	line1 = file.readline()
	label, value = line1.split()
	time = float(value)*10.0*u.Myr
	line2 = file.readline()
	label, value = line2.split()
	tot_num = int(value)
	file.close()

