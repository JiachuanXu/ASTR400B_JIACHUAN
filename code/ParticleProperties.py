# ParticleProperties.py 
# To retrive the information of certain particle in certain snapshot
# Including position(lyr), velocity(km/s), mass(1e10M_sun)
# Usage: $ python ParticleProperties <filename> <type> <sequence>
# By Jiachuan Xu, Jan 21, 2018
import numpy as np
import astropy.units as u
import sys, string
from ReadFile import Read

# The ParticleInfo function returns the info of certain particle 
# Input: particle type, particle sequence in that type, particles data
# particle type: 1-Dark Matter, 2-Disk Stars, 3-Bulge Stars
# particles data: readed from simulation data via Read in ReadFile
# Output: 3D position[kpc], 3D velocity[km/s], particle mass[M_sun]
def ParticleInfo(particle_type, particle_seq, data):

	# Select the particle
	index = np.where(data['type'] == particle_type)
	target = index[0][particle_seq]

	# Read the info (3 decimal places)
	x_tar = np.around(data['x'][target], 3)
	y_tar = np.around(data['y'][target], 3)
	z_tar = np.around(data['z'][target], 3)
	vx_tar = np.around(data['vx'][target], 3)
	vy_tar = np.around(data['vy'][target], 3)
	vz_tar = np.around(data['vz'][target], 3)
	_mass = np.around(data['m'][target]*1e10, 3)*u.solMass
	# Position and velocity triplet
	_position = u.Quantity([x_tar, y_tar, z_tar], u.kpc)
	_velocity = u.Quantity([vx_tar, vy_tar, vz_tar], u.km/u.s)

	return _position, _velocity, _mass

# Main program
# Confirm the arguments
if len(sys.argv) != 4:
	print("Usage: $ python ParticleProperties <filename> <type> <sequence>")
	exit()

# Read filename, particle type, particle sequence from arguments
filename = sys.argv[1]
particle_type = int(sys.argv[2])
particle_seq = int(sys.argv[3])

# Extract the sequence number and define particle type 
snap_seq = filename.split("_")[-1]
snap_seq = int(snap_seq.split(".")[0])
typename = ["", "Dark Matter", "Disk Star", "Bulge Stars"]

# Read file
time, tot_num, data = Read(filename)

# Extract particle info 
position, velocity, mass = ParticleInfo(particle_type, particle_seq, data)

# Print the results
print("You are inquiring about the No.%d %s type particle in SnapNumber %d" \
	  % (particle_seq, typename[particle_type], snap_seq))
print("Position: ", np.around(position.to(u.lyr),3))
print("Velocity: ", velocity)
print("Mass:     ", mass)