# RotationCurve.py
# To Calculate mass profile and rotation curve according to M(r)
# assuming spherical symmetry and circular orbit
# Usage: $ python RotationCurve.py <galaxy name> <snapnumber>
# Output: MassProfile.png, RotationCurve.png
# By Jiachuan Xu, Feb 13, 2018

import numpy as np 
import astropy.units as u 
import sys
import matplotlib.pyplot as plt 
from astropy.constants import G
from ReadFile import Read 
from CenterOfMass import CenterOfMass
from GalaxyMass import ComponentMass
# Gravity Constant 
G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
# Set data dir
def DATADIR():
	return ".."
# Set tolerance for COM_P
def TOLERANCE():
	return 100
# Set index-galaxy component mapping
def INDEX2STR(ptype):
	if ptype==1:
		return "Dark Matter"
	if ptype==2:
		return "Disk Star"
	if ptype==3:
		return "Bulge Star"
	else:
		return -1

# Mass Profile class
# Containing info about mass profile (component and total) and 
# rotation curve(component and total)
class MassProfile:

	# initialization of the class object
    # input: Galaxy Name, Snapnumber
	def __init__(self, galaxyname, snapnumber):
	 	# save the galaxy name as a global
		self.gname = galaxyname
    	# Compose the filename from galaxy name and snapnumber 
		self.filename = "%s/%s_%03d.txt"%(DATADIR(), galaxyname, snapnumber)
        # read in time, total particle number, and other data
		self.time,self.total,self.data = Read(self.filename)
        # store the mass, positions of all particles, and set the units
		self.m = self.data['m']
		self.x = self.data['x']
		self.y = self.data['y']
		self.z = self.data['z']
		self.position = np.column_stack([self.x, self.y, self.z])

		self.x *= u.kpc
		self.y *= u.kpc
		self.z *= u.kpc
		self.position *= u.kpc
		
	# calculate mass enclosed at r kpc (COM frame)
	# input: particle type, radii array [kpc]
	# output: enclosed mass [Msun], stat
	# stat: if galaxy doesn't have the component, stat=0
	# if not, stat = 1
	def MassEnclosed(self, ptype, radii):
		# find ptype particles. If not finded, print caveat and return
		ptype_index, = np.where(self.data['type']==ptype)
		if len(ptype_index)==0:
			print("CAVEAT: NO %s component in %s!"%(INDEX2STR(ptype), self.gname))
			return 0, 0
		# to ensure radii is sorted
		radii = np.sort(radii)
		# get COM position and calculate R with regard to COM 
		center = CenterOfMass(self.filename, ptype)
		COMP = center.COM_P(TOLERANCE())
		R = np.sqrt(np.sum((self.position[ptype_index]-COMP)*(self.position[ptype_index]-COMP), 1))
		# loop through radii to deliver mass into bins
		self.massary = np.zeros(len(radii))
		mass = self.m[ptype_index]
		for i in np.arange(0,len(radii),1):
			ind = np.where(R<radii[i])
			self.massary[i] = np.sum(mass[ind])
		# set units
		self.massary *= 1e10
		self.massary *= u.Msun

		return self.massary, 1
	
	# calculate total mass enclosed at r [kpc]
	# input: radii[kpc]
	# output: total mass enclosed [Msun]
	def MassEnclosedTotal(self, radii):
		# initialize mass bin, and get components' mass 
		self.totmassary = np.zeros(len(radii))
		m1, skip = self.MassEnclosed(1, radii)
		m2, skip = self.MassEnclosed(2, radii)
		m3, skip = self.MassEnclosed(3, radii)
		# sum, return
		self.totmassary = m1+m2+m3

		return self.totmassary
	
	# calculate mass profile according to Hernquist sphere
	# input: radius [kpc], a [kpc], Mhalo [Msun]
	# output: mass encluded [Msun]
	def HernquistMass(self, radius, a, Mhalo):
		# calculate according to forulae and return
		frac = (radius/(radius+a))**2.
		Menclosed = Mhalo*frac*u.Msun

		return Menclosed
	
	# calculate circular velocity assuming circular orbit using mass
	# included at r
	# input: particle type, radii [kpc]
	# output: circular velocity [km/s], stat
	# stat: if do not have this component, stat=0, if not, stat=1
	def CircularVelocity(self, ptype, radii):
		# retrive mass enclosed 
		Menclosed, stat = self.MassEnclosed(ptype, radii)
		# depending on stat, calculate circular velocity
		if stat:vcir = np.around(np.sqrt(G*Menclosed/(radii*u.kpc)),2)
		else: vcir = 0.

		return vcir, stat

	# calculate total circular velocity 
	# input: radii [kpc]
	# output: circular velocity [km/s]
	def CircularVelocityTotal(self, radii):
		Menclosed = self.MassEnclosedTotal(radii)
		vcir = np.around(np.sqrt(G*Menclosed/(radii*u.kpc)))

		return vcir
	
	# calculate circular velocity according to Hernquist profile
	# assuming circular orbit
	# input: radius [kpc], scale factor a [kpc], Mhalo [Msun]
	# output: circular velocity [km/s]
	def HernquistVCirc(self, radius, a, Mhalo):
		vcir = np.sqrt(G*Mhalo*radius/(a+radius)**2.)

		return vcir 


# Main program
# Confirm the arguments
if len(sys.argv) != 3:
    print("Usage: $ RotationCurve.py <galaxy name> <snap number>")
    exit()

# initialize the figure for mass profile 
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# initialize a MassProfile class object
galaxy = MassProfile(sys.argv[1], int(sys.argv[2]))
# set radii array
radii = np.arange(0.1, 30.1, 0.1)*u.kpc
# calculate all components' mass profile, and plot them
MP_DM, DMSTAT = galaxy.MassEnclosed(1, radii)
MP_DS, DSSTAT = galaxy.MassEnclosed(2, radii)
MP_BS, BSSTAT = galaxy.MassEnclosed(3, radii)
MP_Tot = galaxy.MassEnclosedTotal(radii)

if DMSTAT:plt.semilogy(radii, MP_DM, c='r', linestyle="--", label="Dark Matter")
if DSSTAT:plt.semilogy(radii, MP_DS, c='b', linestyle="--", label="Disk Star")
if BSSTAT:plt.semilogy(radii, MP_BS, c='g', linestyle="--", label="Bulge Star")
plt.semilogy(radii, MP_Tot, c='black', linestyle="-", label="Total Mass")

# presume a, and calculate Hernquist profile 
# according to total dark matter halo mass
Mhalo = ComponentMass(galaxy.filename, 1)
a = 20.*u.kpc
MP_Hernquist = galaxy.HernquistMass(radii, a, Mhalo)
plt.semilogy(radii, MP_Hernquist, c='pink', linestyle=":", label="Hernquist: a=%.3fkpc"%a.value)

# set x, y axis label
plt.xlabel("Radii From COM [kpc]", fontsize=22)
plt.ylabel("Mass Enclosed [Msun]", fontsize=22)

# add a legend with some customizations.
legend = ax.legend(loc='lower right',fontsize='x-large')

# add figure title
fig.suptitle('Mass Profile of %s at Snapnumber = %d'%(sys.argv[1], int(sys.argv[2])), fontsize=22)
# save figure
fig.savefig("MassProfile_%s_%03d.png"%(sys.argv[1], int(sys.argv[2])), dpi=350)

# initialize the figure for rotation curve
# for the comments, it's the same with mass profile part 
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

RC_DM, DMSTAT = galaxy.CircularVelocity(1, radii)
RC_DS, DSSTAT = galaxy.CircularVelocity(2, radii)
RC_BS, BSSTAT = galaxy.CircularVelocity(3, radii)
RC_Tot = galaxy.CircularVelocityTotal(radii)

if DMSTAT:plt.semilogy(radii, RC_DM, c='r', linestyle="--", label="Dark Matter")
if DSSTAT:plt.semilogy(radii, RC_DS, c='b', linestyle="--", label="Disk Star")
if BSSTAT:plt.semilogy(radii, RC_BS, c='g', linestyle="--", label="Bulge Star")
plt.semilogy(radii, RC_Tot, c='black', linestyle="-", label="Total Mass")

Mhalo = ComponentMass(galaxy.filename, 1)

RC_Hernquist = galaxy.HernquistVCirc(radii, a, Mhalo)
plt.semilogy(radii, RC_Hernquist, c='pink', linestyle=":", label="Hernquist: a=%.3fkpc"%a.value)

plt.xlabel("Radii From COM [kpc]", fontsize=22)
plt.ylabel("Velocity [km/s]", fontsize=22)

legend = ax.legend(loc='lower right',fontsize='x-large')

fig.suptitle('Rotation Curve of %s at Snapnumber = %d'%(sys.argv[1], int(sys.argv[2])), fontsize=22)

fig.savefig("RotationCurve_%s_%03d.png"%(sys.argv[1], int(sys.argv[2])), dpi=350)