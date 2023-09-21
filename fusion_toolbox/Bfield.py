import numpy as np
import matplotlib.pyplot as plt

class Coil:
    def __init__(self,pts,I):
        """
        Coil Class.

        Parameters
        ----------
        pts : np.ndarray of shape (N,3)
            List of x,y,z triplets which define a closed current loop. The last entry should be
            equal to the first.
        I : float
            Coil current, in Amps
        """
        if not np.allclose(pts[0],pts[-1]): #If loop is not properly closed, connect first/last pts
            pts = np.append(pts,pts[0][None,:],axis=0)
        self.pts = pts
        self.I = I

    def B(self,xyz):
        """
        Returns the magnetic field 

        Parameters
        ----------
        pts : np.ndarray with last dimension of size 3

        Returns
        -------
        np.ndarray of 3D vectors at each point of xyz. Shape is (*xyz.shape,3)
        """
        B_out = np.zeros((len(xyz),3)) #Return array of B vectors in same shape as input

        drs = self.pts[1:]-self.pts[:-1]

        for i,pt in enumerate(self.pts[:-1]): #Skip last point because it repeats
            B_out += self.BGreen(xyz,pt,drs[i])

        return B_out

    def BGreen(self, xyz_samples, xyz_center, dl):
        """Evaluates the B field at all sample locations due to a wire segment. Uses the formula:
        dB = mu_0/4pi * I * dl x (r-r_0)/|r-r_0|^3

        Parameters
        ----------
        xyz_samples : np.ndarray of shape (N,3)
            Locations to evaluate B at

        xyz_center : np.ndarray of shape (3)
            Location of wire segment

        dl : np.ndarray of shape (3)
            Vector for wire segment. Should have magnitude equal to length of wire segment.

        """
        _r = xyz_samples - xyz_center[None,:]
        Bvecs = np.cross(dl,_r)/np.linalg.norm(_r,axis=1)[:,None]**3
        return 1e-7 * self.I * Bvecs #mu_0/4pi = 1e-7 H/m
    def plot(self):
        """
        Displays a 3D plot of the coil
        """
        ax = plt.figure().add_subplot(projection='3d')
        x,y,z = self.pts.T
        ax.plot(x,y,z)
        bbox_val = np.max( (np.max(np.abs(self.pts.T[0])),
                            np.max(np.abs(self.pts.T[1])),
                            np.max(np.abs(self.pts.T[2])))
                         )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_xlim(-bbox_val,bbox_val)
        ax.set_ylim(-bbox_val,bbox_val)
        ax.set_zlim(-bbox_val,bbox_val)
        limits = np.array([getattr(ax, f'get_{axis}lim')() for axis in 'xyz'])
        ax.set_box_aspect(np.ptp(limits, axis = 1))

        plt.show()

class PFCoil(Coil):
    def __init__(self,R,Z0,I,npts=100):
        thetas = np.linspace(0,2*np.pi,npts)
        X = R*np.cos(thetas)
        Y = R*np.sin(thetas)
        Z = Z0*np.ones(npts)
        self.pts = np.vstack((X,Y,Z)).T
        self.I = I

class TFCoil(Coil):
    def __init__(self,R,r,phi,I,npts=100):
        thetas = np.linspace(0,2*np.pi,npts)
        X = np.cos(phi)*(R+r*np.cos(thetas))
        Y = np.sin(phi)*(R+r*np.cos(thetas))
        Z = r*np.sin(thetas)
        self.pts = np.vstack((X,Y,Z)).T
        self.I = I

'''
thetas = np.linspace(0,2*np.pi,100)
X = 1*np.cos(thetas)
Y = 1*np.sin(thetas)
Z = 0*X
XYZ = np.vstack((X,Y,Z)).T

a = Coil(XYZ,1)
field = a.B(np.array([[0,0,0],[0,1,2]]))
print(field)
'''

a = PFCoil(1,0,1)
field = a.B(np.array([[0,0,0],[0,1,2]]))
print(field)
