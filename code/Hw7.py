# Hw7
# By Jiachuan Xu, Mar 29, 2018
import numpy as np 
from OrbitIntegration import M33AnalyticOrbit
import matplotlib.pyplot as plt

# Calculate analytic orbit 
filename = "M33-M31_Analytic_Orbit.txt"
M33AOrb = M33AnalyticOrbit(filename)
M33AOrb.OrbitIntegrator(0.,0.01,10)

# Read data
analy_data = np.genfromtxt(filename,\
	skip_header=3,names=True,dtype=None)
data_M31 = np.genfromtxt("../Orbit_M31.txt",\
	skip_header=1,names=True,dtype=None)
data_M33 = np.genfromtxt("../Orbit_M33.txt",\
	skip_header=1,names=True,dtype=None)

# Pre-process of simulation data
# Get |R| and |v|
R_M332M31_x = data_M31['x']-data_M33['x']
R_M332M31_y = data_M31['y']-data_M33['y']
R_M332M31_z = data_M31['z']-data_M33['z']
V_M332M31_x = data_M31['vx']-data_M33['vx']
V_M332M31_y = data_M31['vy']-data_M33['vy']
V_M332M31_z = data_M31['vz']-data_M33['vz']
R_M332M31 = np.sqrt(R_M332M31_x**2.+R_M332M31_y**2.+R_M332M31_z**2.)
V_M332M31 = np.sqrt(V_M332M31_x**2.+V_M332M31_y**2.+V_M332M31_z**2.)
anal_R = \
np.sqrt(analy_data['x']**2.+analy_data['y']**2.+analy_data['z']**2.)
anal_V = \
np.sqrt(analy_data['vx']**2.+analy_data['vy']**2.+analy_data['vz']**2.)

# Make |R| figure
fig = plt.figure(figsize=(10,10))
fig.suptitle(\
	"Separation Between M31 And M33 As A Function of Time",\
	fontsize=22)

plt.plot(data_M31['t'],R_M332M31,\
	c='r',linestyle="--",label="Simulation Orbit")
plt.plot(analy_data['t'],anal_R,\
	c='b',linestyle="--",label="Analytic Orbit")
# Set label, legend, then save figure
plt.xlabel("Age (Gyr)",fontsize=22)
plt.ylabel("Distance (kpc)",fontsize=22)
plt.legend(loc="upper right",fontsize=16)

fig.savefig("Distance.png", dpi=350)

# Make |V| figure
fig = plt.figure(figsize=(10,10))
fig.suptitle(\
	"Relative Velocity Between M31 And M33 As A Function of Time",\
	fontsize=22)

plt.plot(data_M31['t'],V_M332M31,\
	c='r',linestyle="--",label="Simulation Orbit")
plt.plot(analy_data['t'],anal_V,\
	c='b',linestyle="--",label="Analytic Orbit")

# Set label, legend, then save figure
plt.xlabel("Age (Gyr)",fontsize=22)
plt.ylabel("Relative Velocity (km/s)",fontsize=22)
plt.legend(loc="upper right", fontsize=16)

fig.savefig("RelativeVelocity.png", dpi=350)

Answer = "\
2. At the first 1-2 Gyr, the analytic result fits well with simulation\
 result. But after 2 Gyr, the distance in simulation become smaller,\
 and M33 start to orbit M31 more quickly, and the distance and velocity\
 diagram \"oscillate\".\n\
3. The effect of MW, because at this time, MW is no further than 600 kpc\
 away from the system, this could bring extra gravity for M33, which pull\
 it back. Also, the tidal field strips M33, which could result in damping\
 of the M33 orbit.\n\
4. Treat the MW-M31 system as two-body system, first calculate the \
phase space change using leap frog sheme but with both MW and M31 \
providing gravity, then calculate the motion of M31 and MW using leap frog.\
 But this could at best apply for pre-merger situation. During merger, I\
 prefer to resort to simulation. For post-merger case, trying to take \
 reminant as a integral part.\
"
print(Answer)