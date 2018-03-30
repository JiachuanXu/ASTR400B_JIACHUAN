# OrbitIntegration.py
# To predict the future trajectory of M33 and M31
# By Jiachuan Xu, Mar,27, 2018
import numpy as np 
import astropy.units as u
from CenterOfMass import CenterOfMass
from GalaxyMass import ComponentMass


# Define gravity constant, in kpc^3/M_sun/Gyr^2
G = 4.498768E-6
# tolerance and volumn decreasing factor for CenterOfMass
delta = 10.
VolDec = 2.
# snap0
mwsnap0 = '../MW_000.txt'
m31snap0 = '../M31_000.txt'
m33snap0 = '../M33_000.txt'

# Define M33 analytical orbit class 
# 
class M33AnalyticOrbit:
	# init: filename: file to store the integrated orbit
	# Also init the phase space, size & mass of different components
	# using snap0 as init condition
	def __init__(self, filename):
		self.filename = filename
		M31 = CenterOfMass(m31snap0,2) # realization of M31
		M33 = CenterOfMass(m33snap0,2) # realization of M33
		R_M31 = M31.COM_P(delta, VolDec) # Absolute position, kpc
		V_M31 = M31.COM_V(15.,R_M31) # Absolute velocity, km/s
		R_M33 = M33.COM_P(delta, VolDec)
		V_M33 = M33.COM_V(15.,R_M33)
		self.x = (R_M33 - R_M31)[0] # relative position of M33 to M31
		self.y = (R_M33 - R_M31)[1] # kpc
		self.z = (R_M33 - R_M31)[2]
		self.vx = (V_M33 - V_M31)[0] # relative velocity '''
		self.vy = (V_M33 - V_M31)[1] # km/s
		self.vz = (V_M33 - V_M31)[2]
		self.rd = 5.*u.kpc # disk scale length, kpc
		self.Mdisk = ComponentMass(m31snap0,2) # disk mass, M_sun
		self.rbulge = 1.*u.kpc # bulge scale length, kpc
		self.Mbulge = ComponentMass(m31snap0,3) # bulge mass, M_sun
		self.rhalo = 60.*u.kpc # determined from Assignment 3, kpc
		self.Mhalo = ComponentMass(m31snap0,1) # halo mass, M_sun

	# Define acceleration form a Hernquist profile
	# used for halo and bulge components
	# input: (input as pure number, rather than astropy quantity)
	# 	M: mass of gravity source >(in M_sun)<
	#	ra: Hernquist profile scale length >(in kpc)<
	#	x,y,z: relative distance bet. M31 & M33 >(in kpc)<
	# direc: dummy variable, specify the direction
	#			0-x 1-y 2-z
	# output: >(as pure number)<
	#	a: acceleration along direc >(in kpc/Gyr^2)<
	def HernquistAccel(self, M, ra, x, y, z, direc):
		r = np.sqrt(x**2.+y**2.+z**2.)
		if direc==0:
			accel = -G*M/(r*(ra+r)**2.)*x 
		elif direc==1:
			accel = -G*M/(r*(ra+r)**2.)*y
		elif direc==2:
			accel = -G*M/(r*(ra+r)**2.)*z
		else:
			print("OrbitIntegration.py: \
				Error: illegal dummy variable direc")
			exit()
		return accel

	# Define acceleration from a Miyamoto-Nagai 1975 profile
	# used for disk component
	# input & output: see HernquistAccel()
	# Except that
	#	rd: scale length of disk >(in kpc)<
	def MiyamotoNagaiAccel(self, M, rd, x, y, z, direc):
		zd = rd/5. # disk scale height
		B = rd+np.sqrt(z**2.+zd**2.)
		R = np.sqrt(x**2.+y**2.)
		if direc==0:
			accel = -G*M/((R**2.+B**2.)**1.5)*x 
		elif direc==1:
			accel = -G*M/((R**2.+B**2.)**1.5)*y
		elif direc==2: # note that z direction is different
			accel = -G*M/((R**2.+B**2.)**1.5)*z*\
			B/np.sqrt(z**2.+zd**2.) 
		else:
			print("OrbitIntegration.py: \
				Error: illegal dummy variable direc")
			exit()
		return accel

	# Define acceleration caused by M31 
	# halo + disk + bulge
	# input: 
	#	x,y,z >(in kpc, but no units, pure number)<
	#	direc: 0-x 1-y 2-z
	# output: 
	#	total acceleration along direc >(in kpc/Gyr^2)<
	def M31Accel(self, x, y, z, direc):
		accel_tot = \
		self.HernquistAccel(self.Mhalo.value, self.rhalo.value, \
			x, y, z, direc)+\
		self.HernquistAccel(self.Mbulge.value, self.rbulge.value, \
			x, y, z, direc)+\
		self.MiyamotoNagaiAccel(self.Mdisk.value, self.rd.value, \
			x, y, z, direc)
		return accel_tot

	# Integration "brick" function: Leap Frog
	# Scheme: x_n,v_n->x_{n+1/2}->a_{n+1/2}->v_{n+1}->x_{n+1}
	# Input: (no units, pure number)
	#	Dt: time interval >(in Gyr)<
	#	x_i,y_i,z_i: initial position >(in kpc)<
	#	vx_i,vy_i,vz_i: initial velocity >(in km/s)<
	# Output:
	#	x_pp, y_pp, z_pp: updated position after Dt >(in kpc)<
	#	vx_pp, vy_pp, vz_pp: updated velocity after Dt >(in km/s)<
	def LeapFrog(self, x_i, y_i, z_i, vx_i, vy_i, vz_i, Dt):
		# factor for velocity: 1*km/s = norm*kpc/Gyr
		norm = 1.022
		# "springboard" position, kpc
		xp = x_i + vx_i*Dt/2.*norm
		yp = y_i + vy_i*Dt/2.*norm
		zp = z_i + vz_i*Dt/2.*norm
		# "springboard" acceleration, kpc/Gyr^2
		ax_p = self.M31Accel(xp, yp, zp, 0)
		ay_p = self.M31Accel(xp, yp, zp, 1)
		az_p = self.M31Accel(xp, yp, zp, 2)
		# updated velocity, kpc/Gyr
		vx_pp = vx_i*norm + ax_p*Dt 
		vy_pp = vy_i*norm + ay_p*Dt 
		vz_pp = vz_i*norm + az_p*Dt 
		# updated position, kpc
		x_pp = x_i + Dt/2.*(vx_i*norm + vx_pp)
		y_pp = y_i + Dt/2.*(vy_i*norm + vy_pp)
		z_pp = z_i + Dt/2.*(vz_i*norm + vz_pp)

		return x_pp, y_pp, z_pp, \
		vx_pp/norm, vy_pp/norm, vz_pp/norm

	# Integration loop
	# Since this is a func in CLASS M33AnalyticOrbit
	# The function would naturally take initializing snapshot as IC
	# And integrates from t_o to t_max
	# With the help of leap frog 
	# We would finall reunite at the eternity
	# Input: (no units, pure number)
	#	t_o: start time, >(Gyr)<
	#	Dt: time interval >(Gyr)<
	#	t_max: end time >(Gyr)<
	# Output:
	#	<orbit file>: a file, named after self.filename,
	#				  recording the trajectory
	# Writing comments
	# makes me feel like
	# a poet
	def OrbitIntegrator(self, t_o, Dt, t_max):
		# construct a stack for t, position and velocity
		# initialize with IC
		x_temp = np.array([self.x.value,]) # in kpc
		y_temp = np.array([self.y.value,]) # in kpc
		z_temp = np.array([self.z.value,]) # in kpc
		vx_temp = np.array([self.vx.value,]) # in km/s
		vy_temp = np.array([self.vy.value,]) # in km/s
		vz_temp = np.array([self.vz.value,]) # in km/s
		t_temp = np.array([t_o,])
		# initialize t, we are gonna to start our journey
		i = 0
		while t_temp[i]<t_max:
			x_l, y_l, z_l, vx_l, vy_l, vz_l = \
			self.LeapFrog(x_temp[i],y_temp[i],z_temp[i],\
				vx_temp[i],vy_temp[i],vz_temp[i],Dt)
			t_l = t_temp[i]+Dt
			# update posi & velo 
			x_temp = np.append(x_temp, x_l)
			y_temp = np.append(y_temp, y_l)
			z_temp = np.append(z_temp, z_l)
			vx_temp = np.append(vx_temp, vx_l)
			vy_temp = np.append(vy_temp, vy_l)
			vz_temp = np.append(vz_temp, vz_l)
			# update time
			t_temp = np.append(t_temp, t_l)
			# update indice
			i += 1
		# reach our destination, write our memory
		size = i+1
		fp = open(self.filename,'w')
		header = "# time: Gyr\n# x,y,z: kpc\n# vx,vy,vz: km/s\n"
		nameline = "# t x y z vx vy vz\n"
		fp.write(header+nameline)
		for i in np.arange(size):
			fp.write("%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n"\
				%(t_temp[i],x_temp[i],y_temp[i],z_temp[i],\
					vx_temp[i],vy_temp[i],vz_temp[i]))
		fp.close()