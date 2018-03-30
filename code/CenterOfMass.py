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
    # input: tolerance, volumn decreasing factor
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
        # devided by VolDec, shrink the range to do refine calculation
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
