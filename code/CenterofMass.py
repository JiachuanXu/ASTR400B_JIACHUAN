# CenterofMass.py 
# To get the Center of Mass (COM) Position and Velocity
# Usage: $python CenterOfMass.py <snapnumber> <ptype> <tolerance> <rad>
# snapnumber: the number of snapshot
# ptype: 1-dark matter 2-disk star 3-bulge star
# tolerance: the maximum difference between RCOM and RCOM2 in COM 
#            position iteration, in kpc
# rad: within the radius in COM frame would this program calculate the
#      COM velocity, in kpc
# By Jiachuan Xu, Feb, 6, 2018

# import modules
import numpy as np
import astropy.units as u
from ReadFile import Read
import sys


# The class that manage the info of COM of specific particle type
class CenterOfMass:

    # initialization of the class object
    # input: (path)filename, particle type
    def __init__(self, filename, ptype):
        # read in the file and particle type
        self.time, self.total, self.data = Read(filename)
            
        #create an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)
    
        # store the mass, positions, velocities of particles
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
        self.position = np.column_stack([self.x, self.y, self.z])
        self.velocity = np.column_stack([self.vx, self.vy, self.vz])
     
    # calculate the total mass 
    def total_mass(self):
        return np.sum(self.m)*u.Msun*1e10
    
    # a general function to calculate COM 
    # input: weight, coordinate
    # weight: mass array
    # coordinate: position or velocity array
    def COMdefine(self, weight, coord):
        # x, y, z component of COM 
        r_cx = np.sum(weight*coord[:,0])
        r_cy = np.sum(weight*coord[:,1])
        r_cz = np.sum(weight*coord[:,2])
        r_c = np.array([r_cx, r_cy, r_cz])
        r_c /= np.sum(weight)

        # return as ndarray
        return r_c

    # a more sophisticated function to calculate COM position 
    # input: tolerance
    def COM_P(self, delta, VolDec):
        # first estimate of COM position
        [XCOM, YCOM, ZCOM] = self.COMdefine(self.m, self.position)
        # Calculate the magnitude of COM position
        RCOM = np.sqrt(XCOM**2.+YCOM**2.+ZCOM**2.)

        # transfer into COM frame
        COMframe_position = self.position-np.array([XCOM, YCOM, ZCOM])
        # magnitude of all the particles' position in COM frame
        RNEW = np.sqrt(np.sum(COMframe_position*COMframe_position, 1))
        # select the maximum if particle position in COM frame
        RMAX = RNEW.max()
        # devided by VolDec(volumn decrease), shrink the range to do refine calculation
        RMAX /= VolDec
        index_refine = np.where(RNEW<RMAX)
        # calculate refined COM position and its magnitude
        [XCOM_NEW, YCOM_NEW, ZCOM_NEW] = \
            self.COMdefine(self.m[index_refine],self.position[index_refine])
        RCOM2 = np.sqrt(XCOM_NEW**2.+YCOM_NEW**2.+ZCOM_NEW**2.)

        # if the difference between two RCOM is smaller than tolerance,
        # we regard the COM position to be converged.
        # if not, continue the while loop till converged
        while np.abs(RCOM2-RCOM)>delta:
            # replace the old COM with the new one
            RCOM = RCOM2
            XCOM, YCOM, ZCOM = XCOM_NEW, YCOM_NEW, ZCOM_NEW
            # repeat the former procedure to iterate
            COMframe_position = \
                self.position[index_refine]-np.array([XCOM, YCOM, ZCOM])
            RNEW = np.sqrt(np.sum(COMframe_position*COMframe_position, 1))
            RMAX = RNEW.max()
            RMAX /= VolDec
            index_refine = np.where(RNEW<RMAX)
            [XCOM_NEW, YCOM_NEW, ZCOM_NEW] = self.COMdefine\
                (self.m[index_refine], self.position[index_refine])
            RCOM2 = np.sqrt(XCOM_NEW**2.+YCOM_NEW**2.+ZCOM_NEW**2.)

        # return converged COM position in kpc
        return np.array([XCOM, YCOM, ZCOM])*u.kpc
    
    # a function to calculate the COM velocity within specific rigion
    # with regard to COM
    # input: radius, COM position 
    # radius: the radius within will the program calculate v of COM
    # COM position: ndarray of COM position
    def COM_V(self, radius, RC):
        # transfer into COM frame
        COMframe_position = self.position*u.kpc-RC
        # calculate magnitude of position in COM frame
        RNEW = np.sqrt(np.sum(COMframe_position*COMframe_position, 1))
        # calculate COM velocity within the region
        index = np.where(RNEW < radius*u.kpc)
        VXCOM, VYCOM, VZCOM = \
            self.COMdefine(self.m[index], self.velocity[index])
        # return as ndarray in km/s
        return np.array([VXCOM, VYCOM, VZCOM])*u.km/u.s

# Main program
# Confirm the arguments
if len(sys.argv) != 5:
    print("Usage: $ CenterOfMass.py <snapnumber> <ptype> <tolerance> <rad>")
    exit()

# Create a Center of mass object for the MW, M31, M33 
MWCOM = CenterOfMass("../MW_%03d.txt"%int(sys.argv[1]), int(sys.argv[2]))
M31COM = CenterOfMass("../M31_%03d.txt"%int(sys.argv[1]), int(sys.argv[2]))
M33COM = CenterOfMass("../M33_%03d.txt"%int(sys.argv[1]), int(sys.argv[2]))
# Calculate Center of Mass position for the MW, M31, M33 
rc_MW = MWCOM.COM_P(float(sys.argv[3]), 2.)
rc_M31 = M31COM.COM_P(float(sys.argv[3]), 2.)
rc_M33 = M33COM.COM_P(float(sys.argv[3]), 2.)
# Calculate Center of Mass velocity for the MW, M31, M33 
vc_MW = MWCOM.COM_V(float(sys.argv[4]), rc_MW)
vc_M31 = M31COM.COM_V(float(sys.argv[4]), rc_M31)
vc_M33 = M33COM.COM_V(float(sys.argv[4]), rc_M33)
# Calculate disk mass for MW, M31, M33 data
MW_mass = MWCOM.total_mass()
M31_mass = M31COM.total_mass()
M33_mass = M33COM.total_mass()
print("MW:\n\tDisk Mass:", np.around(MW_mass.to(1e10*u.Msun), 3), \
    "\n\tCOM Position:",np.around(rc_MW, 3), \
    "\n\tCOM Velocity:", np.around(vc_MW, 3))
print("M31:\n\tDisk Mass:", np.around(M31_mass.to(1e10*u.Msun), 3), \
    "\n\tCOM Position:", np.around(rc_M31, 3), \
    "\n\tCOM Velocity:", np.around(vc_M31, 3))
print("M33:\n\tDisk Mass:", np.around(M33_mass.to(1e10*u.Msun), 3), \
    "\n\tCOM Position:", np.around(rc_M33, 3), \
    "\n\tCOM Velocity:", np.around(vc_M33, 3))
# Calculate current separation between MW and M31 
MWtoM31_r = rc_MW - rc_M31
MWtoM31_v = vc_MW - vc_M31
MWtoM31_rmag = np.sqrt(np.sum(MWtoM31_r*MWtoM31_r))
MWtoM31_vmag = np.sqrt(np.sum(MWtoM31_v*MWtoM31_v))
MWtoM31_vrad = np.sum(MWtoM31_v*MWtoM31_r)/MWtoM31_rmag
MWtoM31_vtan = np.sqrt(MWtoM31_vmag**2.-MWtoM31_vrad**2.)
print("MW v.s. M31:\n\tDistance:", np.around(MWtoM31_rmag, 3), \
    "\n\tRadial Velocity:", np.around(MWtoM31_vrad, 3), \
    "\n\tTangential Velocity:", np.around(MWtoM31_vtan, 3))
# Calculate current separation between M33 and M31
M33toM31_r = rc_M33 - rc_M31
M33toM31_v = vc_M33 - vc_M31
M33toM31_rmag = np.sqrt(np.sum(M33toM31_r*M33toM31_r))
M33toM31_vmag = np.sqrt(np.sum(M33toM31_v*M33toM31_v))
M33toM31_vrad = np.sum(M33toM31_v*M33toM31_r)/M33toM31_rmag
M33toM31_vtan = np.sqrt(M33toM31_vmag**2.-M33toM31_vrad**2.)
print("M33 v.s. M31:\n\tDistance:", np.around(M33toM31_rmag, 3), \
    "\n\tRadial Velocity:", np.around(M33toM31_vrad, 3), \
    "\n\tTangential Velocity:", np.around(M33toM31_vtan, 3))

# for homework 4 specific
#print("Answer to question 4: When two galaxie are about to merge, their \
#outliers could be highly extended and irregular due to the interaction\
# between the two galaxies. In that case, the center of a galaxy and \
# its outliers do not perform like an integral entity, so a center of \
# mass describing the more clustered inner part of the galaxy is more \
# useful, since the COM concept works well for a well-gavity-bounded \
# system. After iteration, we could constraint our focus on a more \
# concentrated galaxy")