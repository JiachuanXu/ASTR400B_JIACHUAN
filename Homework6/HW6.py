import numpy as np 
import matplotlib.pyplot as plt

data_MW = np.genfromtxt("Orbit_MW.txt",skip_header=1,names=True,dtype=None)
data_M31 = np.genfromtxt("Orbit_M31.txt",skip_header=1,names=True,dtype=None)
data_M33 = np.genfromtxt("Orbit_M33.txt",skip_header=1,names=True,dtype=None)

R_MW2M31_x = data_MW['x']-data_M31['x']
R_MW2M31_y = data_MW['y']-data_M31['y']
R_MW2M31_z = data_MW['z']-data_M31['z']
V_MW2M31_x = data_MW['vx']-data_M31['vx']
V_MW2M31_y = data_MW['vy']-data_M31['vy']
V_MW2M31_z = data_MW['vz']-data_M31['vz']
R_M332M31_x = data_M31['x']-data_M33['x']
R_M332M31_y = data_M31['y']-data_M33['y']
R_M332M31_z = data_M31['z']-data_M33['z']
V_M332M31_x = data_M31['vx']-data_M33['vx']
V_M332M31_y = data_M31['vy']-data_M33['vy']
V_M332M31_z = data_M31['vz']-data_M33['vz']

R_MW2M31 = np.sqrt(R_MW2M31_x**2.+R_MW2M31_y**2.+R_MW2M31_z**2.)
R_M332M31 = np.sqrt(R_M332M31_x**2.+R_M332M31_y**2.+R_M332M31_z**2.)

V_MW2M31 = np.sqrt(V_MW2M31_x**2.+V_MW2M31_y**2.+V_MW2M31_z**2.)
V_M332M31 = np.sqrt(V_M332M31_x**2.+V_M332M31_y**2.+V_M332M31_z**2.)

fig = plt.figure(figsize=(10,10))
fig.suptitle("Magnitude of Separation As A Function of Time")

plt.plot(data_M31['t'],R_MW2M31,c='b',linestyle="--",label="MW v.s. M31")
plt.plot(data_M31['t'],R_M332M31,c='r',linestyle="--",label="M31 v.s. M31")

plt.xlabel("Age (Gyr)",fontsize=22)
plt.ylabel("Distance (kpc)",fontsize=22)

plt.legend(loc="upper right")
fig.savefig("Distance.png", dpi=350)

fig = plt.figure(figsize=(10,10))
fig.suptitle("Magnitude of Relative Velocity As A Function of Time")

plt.plot(data_M31['t'],V_MW2M31,c='b',linestyle="--",label="MW v.s. M31")
plt.plot(data_M31['t'],V_M332M31,c='r',linestyle="--",label="M31 v.s. M31")

plt.xlabel("Age (Gyr)",fontsize=22)
plt.ylabel("Relative Velocity (km/s)",fontsize=22)

plt.legend(loc="upper right")
fig.savefig("RelativeVelocity.png", dpi=350)

# find decay rate 
localM_R = np.array([])
localM_t = np.array([])
index = np.where(data_M31['t']>6)
index = index[0]
for i in index[0:len(index)-1]:
	if (R_M332M31[i]>R_M332M31[i-1])and(R_M332M31[i]>R_M332M31[i+1]):
		localM_t = np.append(localM_t, data_M31['t'][i])
		localM_R = np.append(localM_R, R_M332M31[i])
		print("local max founded")
v_decay = -(localM_R[-1]-localM_R[0])/(localM_t[-1]-localM_t[0])
print("mean decay velocity of M33: %f kpc/Gyr"%v_decay)

# For HW 6 specific
print("Answers:\n\t1. 3 close encounters before merger\n\
	2. where distance decrease, velocity increase\n\
	3. Around 6.2 Gyr, M33's orbit decrease when they merge \
	(a little non-continuous in first derivative), maybe some\
	collision happens between M33 and merger reminant which \
	change the kinetic and/or morphology of M33\n\
	 4. Around 15 kpc/Gyr, for 75 kpc distance, it takes about\
	  5 Gyr to merge.")

