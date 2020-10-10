#!/usr/bin/env python
# coding: utf-8


from numpy import cos, sin, pi, sqrt
from scipy.integrate import quad


# Computation of the magnetic field distribution of a diametrical 
# magnetized cylinder with real parameters

class B_flux():
    """
    This class computes the three components of the magnetic field distribution
    of a diametrically magnetised cylinder with the cylindrical coordinate system
    where the origin is at the revolution axis and middle of its thickness.
    Inputs for this class include:
    J (T) - the remance (input unit Tesla to yield an output in Tesla
    or J*10**4 to yield an output in Gauss)
    R (m) - the radius of the cylinder
    h (m) - the thickness of the cylinder
    r (m) - the radial coordinate of the computed point
    z (m) - the z coordinate of the computed point
    alpha (rad) - the azimuthal angle of the computed point
    Axial, Azimuthal and Radial denote the axial, azimuthal and radial components
    of the magnetic field
    """
    def __init__(self, J, R, h, r, z, alpha):
        self.J = J
        self.R = R
        self.h = h
        self.r = r
        self.z = z
        self.alpha = alpha
    def PK_squared(self,beta):
        return self.R**2+self.r**2 - 2*self.R*self.r*cos(beta)
    def main_func_axial(self,beta):
        func_1 = 1/sqrt((self.h/2-self.z)**2+self.PK_squared(beta))
        func_2 = 1/sqrt((self.h/2+self.z)**2+self.PK_squared(beta))
        return (func_1-func_2)*cos(beta)
    def Axial(self):
        return self.J*self.R*quad(self.main_func_axial,-pi,pi)[0]*sin(self.alpha)/4/pi
    def main_func_azimuthal(self,beta):
        func_1 = (self.h/2-self.z)/(self.PK_squared(beta)*sqrt((self.h/2-self.z)**2+self.PK_squared(beta)))
        func_2 = (self.h/2+self.z)/(self.PK_squared(beta)*sqrt((self.h/2+self.z)**2+self.PK_squared(beta)))
        return (func_1 + func_2)*-self.R*sin(beta)**2
    def Azimuthal(self):
        return self.J*self.R*quad(self.main_func_azimuthal,-pi,pi)[0]*cos(self.alpha)/4/pi
    def main_func_radial(self,beta):
        func_1 = (self.h/2-self.z)/(self.PK_squared(beta)*sqrt((self.h/2-self.z)**2+self.PK_squared(beta)))
        func_2 = (self.h/2+self.z)/(self.PK_squared(beta)*sqrt((self.h/2+self.z)**2+self.PK_squared(beta)))
        return (func_1 + func_2)*(self.r*cos(beta)-self.R*cos(beta)**2)
    def Radial(self):
        return self.J*self.R*quad(self.main_func_radial,-pi,pi)[0]*sin(self.alpha)/4/pi

        
        


# In[ ]:





# In[89]:


10**4*B_flux(1, 2.5/1000,5/1000,3/1000,2/1000,pi/3).Radial()


# In[ ]:




