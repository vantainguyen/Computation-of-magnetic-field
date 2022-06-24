#!/usr/bin/env python
# coding: utf-8

# ## Computation of the magnetic field of a cone with axial magnetization

# In[26]:


# Import necessary library
from numpy import sqrt, sin, cos, pi,log,tan
from scipy.integrate import quad, dblquad, tplquad
# B_axial = B_KUZ + B_HCZ
mu0 = 4*pi*10**-7
# Computation of H_KUZ

# Auxilliary function 1
"""
Parameter definition

The magnetization is axially oriented (V. T. Nguyen, “Magnetic field 
distribution of a conical permanent magnet with an application in magnetic resonance imaging,” 
J. Magn. Magn. Mater. 498, 166136 (2019).)

The function to compute the magnetic field is denoted as: 
B_cone_axial_axial where the meaning of the abbreviations as B - magnetic flux density; 
cone - conical geometry; axial - axial orientation of the magnetization; 
axial - the axial component of the magnetic flux density

J: the remanence in Tesla of the lower cone
J_1: the remanance in Tesla of the upper cone
r_b: the base radius of the lower cone
L: the height of the lower cone
2*alpha: the outer apex angle of the lower annular cone 
2*beta: the inner apex angle of the lower annular cone
2*phi: the apex angle of the upper cone
K: the height of the upper cone
(r, z, delta): the coordinates of the computed point in the Cylindrical coordinate system
pxi: the distance along Z axis between the two apexes 
(auxiliary parameter khi denotes the slant height of the cone)
"""

def main_func_KUZ(theta, r_b, L, r, z, delta):
    pxi = 2*r*cos(delta-theta)
    khi = r**2 + (z-L)**2
    main_1 = 2*(pxi*r_b - 2*khi)/(4*khi - pxi**2)/sqrt(r_b*(r_b - pxi)+khi)
    main_2 = 4*sqrt(khi)/(4*khi - pxi**2)
    return (main_1+main_2)*(z-L)

# B_KUZ

def B_KUZ(J,r_b, L, r, z, delta):
    return J*quad(lambda theta: main_func_KUZ(theta, r_b, L, r, z, delta),-pi,pi)[0]/4/pi

# Auxilliary function 2

def main_func_HCZ(theta, alpha, r_b, L, r, z, delta):
    var_1 = sqrt(L**2+r_b**2)
    var_2 = 2*r*sin(alpha)*cos(delta - theta)+2*z*cos(alpha)
    var_3 = r**2 + z**2
    main_1 = (2*(z*var_2+2*var_3*cos(alpha)-var_2**2*cos(alpha))*var_1-
              4*var_3*z + 2*var_2*var_3*cos(alpha))/(4*var_3-var_2**2)/sqrt(var_1*(var_1-var_2)+var_3)
    main_2 = (-4*var_3*z+2*var_2*var_3*cos(alpha))/(4*var_3-var_2**2)/sqrt(var_3)
    main_3 = cos(alpha)*log((2*(sqrt(var_1*(var_1-var_2)+var_3)+var_1)-var_2)/(2*sqrt(var_3)-var_2))
    return main_1 - main_2 - main_3

# B_HCZ

def B_HCZ(J, alpha, r_b, L, r, z, delta):
    return -J*sin(alpha)**2*quad(lambda theta: main_func_HCZ(theta, alpha, r_b, L, r, z, delta),
                                 -pi,pi)[0]/4/pi

# B_cone_axial
# phi is the half convex angle of the cone
# delta is the coordinate angle of the computed point
def B_cone_axial_axial(J, alpha, L, r, z, delta):
    r_b = L*tan(alpha)
    return B_KUZ(J,r_b, L, r, z, delta) + B_HCZ(J, alpha, r_b, L, r, z, delta)




# Import necessary library
from numpy import sqrt, sin, cos, pi, log, tan
from scipy.integrate import quad, dblquad, tplquad

mu0 = 4*pi*10**-7

"""
Parameter definition

The magnetization is axially oriented (V. T. Nguyen, “Analytical computation of the magnetic field
of a conical permanent magnet with arbitrarily
uniform magnetization,” AIP Advances 10, 045208 (2020).)

The function to compute the magnetic field is denoted as: 
B_cone_diam_axial where the meaning of the abbreviations as B - magnetic flux density; 
cone - conical geometry; diam - diametrical orientation of the magnetization; 
axial - the axial component of the magnetic flux density

J: the remanence in Tesla of the cone
h: the height of the cone
2*phi: the apex angle of the cone 
(r, z, beta): the coordinates of the computed point in the Cylindrical coordinate system
"""

# The axial component of the magnetic field

def main_func_diam_axial(h, R, phi, r, z, beta, alpha):
    
    # declaring the auxilliary parameters
    
    khi = sqrt(h**2 + R**2)
    psi = 2*r*sin(phi)*cos(beta - alpha) + 2*z*cos(phi)
    theta = r**2 + z**2
    
    # denoting some components in the main functions
    # first component
    main_1_num = 2*(z*psi + 2*theta*cos(phi) - psi**2*cos(phi))*khi - 4*theta*z + 2*psi*theta*cos(phi)
    main_1_deno = (4*theta - psi**2)*sqrt(khi*(khi - psi) + theta)
    main_1 = main_1_num/main_1_deno
    
    # second component
    main_2_num = -4*theta*z + 2*psi*theta*cos(phi)
    main_2_deno = (4*theta - psi**2)*sqrt(theta)
    main_2 = main_2_num/main_2_deno
    
    # third component
    main_3_num = 2*(sqrt(khi*(khi - psi) + theta) + khi) - psi
    main_3_deno = 2*sqrt(theta) - psi
    main_3 = cos(phi)*log(main_3_num/main_3_deno)
    
    return (main_1 - main_2 - main_3)*cos(alpha)

def B_cone_diam_axial(J, h, phi, r, z, beta):
    R = h*tan(phi)
    integral_B = quad(lambda alpha: main_func_diam_axial(h, R, phi, r, z, beta, alpha), 0, 2*pi)[0]
    return J*sin(2*phi)/8/pi*integral_B




# ## Computation of the magnetic field from cylinder with axial magnetization

# In[3]:


# Importing relevant libraries

import matplotlib as plt
import math as m
from numpy import sqrt, sin, cos, pi
from scipy.integrate import quad, tplquad

"""
Parameter definition

The magnetization is axially oriented (V. T. Nguyen, T - F. Lu, W. Robertson, P. Grimshaw, “Magnetic Field Distribution of an Elliptical Permanent Magnet,” 
Progress In Electromagnetics Research C, Vol. 97, 69–82, 2019.)

The function to compute the magnetic field is denoted as: 
B_cylinder_axial_axial where the meaning of the abbreviations as B - magnetic flux density; 
cylinder - cylindrical geometry; axial - axial orientation of the magnetization; 
axial - the axial component of the magnetic flux density

J: the remanence in Tesla of the cylinder
h: the height of the cylinder
R: the radius of the cylinder
(r, z, alpha): the coordinates of the computed point in the Cylindrical coordinate system
"""

# Defining the constants
mu0 = 4*pi*10**-7 # the magnetic permeability
K = 1/4/pi/mu0 # auxiliary constant

# Delta funcion
def delta(theta,alpha,r):
    return 2*r*cos(alpha-theta)

# Pxi plus
def pxi_plus(r,z,h):
    return r**2 + (z-h)**2

# Pxi minus
def pxi_minus(r,z,h):
    return r**2 + z**2

# Mag function
def mag(J):
    return J/(4*pi)

# Auxiliary radius
def r0(theta,R):
    a = R
    b = R
     
    return a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)

# Computation of the magnetic field from the upper (plus) surface (dB_plus)
def dB_plus(theta,r, alpha, z, R, h):
    db_func = (2*(delta(theta,alpha,r)*r0(theta,R) - 
                  2*pxi_plus(r,z,h))/(4*pxi_plus(r,z,h) - 
                                      delta(theta,alpha,r)**2)/sqrt((r0(theta,R)*
                                                                     (r0(theta,R)-delta(theta,alpha,r))+
                                                                     pxi_plus(r,z,h))) + 
               4*sqrt(pxi_plus(r,z,h))/\
               (4*pxi_plus(r,z,h)-delta(theta,alpha,r)**2))*(z-h)
    return db_func
    
# Computation of the magnetic field from the lower (minus) surface (dB_minus)
def dB_minus(theta,r, alpha, z, R, h):
    db_func = -(2*(delta(theta,alpha,r)*r0(theta,R) - 
                  2*pxi_minus(r,z,h))/(4*pxi_minus(r,z,h) - 
                                      delta(theta,alpha,r)**2)/sqrt((r0(theta,R)*
                                                                     (r0(theta,R)-delta(theta,alpha,r))+
                                                                     pxi_minus(r,z,h))) + 
                4*sqrt(pxi_minus(r,z,h))/\
                (4*pxi_minus(r,z,h)-delta(theta,alpha,r)**2))*(z)
    return db_func

# Computation of the total magnetic field from both surfaces (dB)
def dB(theta,r, alpha, z, R, h):
    return dB_plus(theta,r, alpha, z, R, h) + dB_minus(theta,r, alpha, z, R, h)

def B_cylinder_axial_axial(J, h, R, r, z, alpha):
    integral_B = quad(lambda theta: dB(theta,r, alpha, z, R, h), -pi, pi)[0]
    return J*integral_B/4/pi




# ## Computation of magnetic field from elliptical cylinder with axial magnetization

# In[17]:


# Importing relevant libraries

import matplotlib as plt
from numpy import sqrt, sin, cos, pi, log
from scipy.integrate import quad, tplquad

"""
Parameter definition

The magnetization is axially oriented (V. T. Nguyen, T - F. Lu, W. Robertson, P. Grimshaw, “Magnetic Field Distribution of an Elliptical Permanent Magnet,” 
Progress In Electromagnetics Research C, Vol. 97, 69–82, 2019.)

The function to compute the magnetic field is denoted as: 
B_ellip_cylinder_axial_axial (azimuthal, radial) where the meaning of the abbreviations as B - magnetic flux density; 
ellip_cylinder - elliptical cylinder; axial - axial orientation of the magnetization; 
axial (azimuthal, radial) - the axial (azimuthal, radial) component of the magnetic flux density

J: the remanence in Tesla of the elliptical cylinder
h: the height of the elliptical cylinder
R: the radius of the elliptical cylinder
(r, z, alpha): the coordinates of the computed point in the Cylindrical coordinate system
"""

# Defining the constants
mu0 = 4*pi*10**-7 # the magnetic permeability


""" Computation of the axial component B_ellip_cylinder_axial_axial"""

# B_ellip_cylinder_axial_axial = B_plus + B_minus

# Computation of B_plus
def main_plus_ellip_axial_axial(h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    delta = 2*r*cos(alpha - theta)
    xi_plus = r**2 + (z - h)**2
    main_1_num = 2*(delta*r0 - 2*xi_plus)
    main_1_deno = (4*xi_plus - delta**2)*sqrt(r0*(r0 - delta) + xi_plus)
    main_1 = main_1_num/main_1_deno
    main_2_num = 4*sqrt(xi_plus)
    main_2_deno = 4*xi_plus - delta**2
    main_2 = main_2_num/main_2_deno
    return (main_1 + main_2)*(z - h)

def B_plus_ellip_axial_axial(J, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_plus_ellip_axial_axial(h, a, b, r, z, alpha, theta), -pi, pi)[0]
    return J*integral_B/4/pi 

# Computation of B_minus
def main_minus_ellip_axial_axial(h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    delta = 2*r*cos(alpha - theta)
    xi_minus = r**2 + z**2
    main_1_num = 2*(delta*r0 - 2*xi_minus)
    main_1_deno = (4*xi_minus - delta**2)*sqrt(r0*(r0 - delta) + xi_minus)
    main_1 = main_1_num/main_1_deno
    main_2_num = 4*sqrt(xi_minus)
    main_2_deno = 4*xi_minus - delta**2
    main_2 = main_2_num/main_2_deno
    return (main_1 + main_2)*(z)

def B_minus_ellip_axial_axial(J, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_minus_ellip_axial_axial(h, a, b, r, z, alpha, theta), -pi, pi)[0]
    return -J*integral_B/4/pi


def B_ellip_cylinder_axial_axial(J, a, b, h, r, z, alpha):
    return B_plus_ellip_axial_axial(J, a, b, h, r, z, alpha) + B_minus_ellip_axial_axial(J, a, b, h, r, z, alpha)

""" Computation of the azimuthal component B_ellip_cylinder_axial_azimuthal"""

# B_ellip_cylinder_axial_azimuthal = B_azi_plus + B_azi_minus

def main_plus_ellip_axial_azimuthal(h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    delta = 2*r*cos(alpha - theta)
    xi_plus = r**2 + (z - h)**2
    main_1 = log(2*(sqrt(r0*(r0 - delta) + xi_plus) + r0) - delta)
    main_2_nomi = 2*(2*xi_plus - delta**2)*r0 + 2*delta*xi_plus
    main_2_deno = (4*xi_plus - delta**2)*sqrt(r0*(r0 - delta) + xi_plus)
    main_2 = main_2_nomi/main_2_deno
    main_3 = log(2*sqrt(xi_plus) - delta)
    main_4_nomi = 2*delta*xi_plus
    main_4_deno = (4*xi_plus - delta**2)*sqrt(xi_plus)
    main_4 = main_4_nomi/main_4_deno
    return (main_1 - main_2 - main_3 + main_4)*sin(alpha - theta)

def B_plus_ellip_axial_azimuthal(J, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_plus_ellip_axial_azimuthal(h, a, b, r, z, alpha, theta), -pi, pi)[0]
    return J*integral_B/4/pi 


def main_minus_ellip_axial_azimuthal(h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    delta = 2*r*cos(alpha - theta)
    xi_minus = r**2 + z**2
    main_1 = log(2*(sqrt(r0*(r0 - delta) + xi_minus) + r0) - delta)
    main_2_nomi = 2*(2*xi_minus - delta**2)*r0 + 2*delta*xi_minus
    main_2_deno = (4*xi_minus - delta**2)*sqrt(r0*(r0 - delta) + xi_minus)
    main_2 = main_2_nomi/main_2_deno
    main_3 = log(2*sqrt(xi_minus) - delta)
    main_4_nomi = 2*delta*xi_minus
    main_4_deno = (4*xi_minus - delta**2)*sqrt(xi_minus)
    main_4 = main_4_nomi/main_4_deno

    return (main_1 - main_2 - main_3 + main_4)*sin(alpha - theta)

def B_minus_ellip_axial_azimuthal(J, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_minus_ellip_axial_azimuthal(h, a, b, r, z, alpha, theta), -pi, pi)[0]
    return -J*integral_B/4/pi

def B_ellip_cylinder_axial_azimuthal(J, a, b, h, r, z, alpha):
    return B_plus_ellip_axial_azimuthal(J, a, b, h, r, z, alpha) + B_minus_ellip_axial_azimuthal(J, a, b, h, r, z, alpha)

""" Computation of the radial component B_ellip_cylinder_axial_radial"""

# B_ellip_cylinder_axial_radial = B_radi_plus + B_radi_minus

def main_plus_ellip_axial_radial(h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    delta = 2*r*cos(alpha - theta)
    xi_plus = r**2 + (z - h)**2
    gamma = cos(alpha - theta)
    main_1_nomi = 2*(delta*r + 2*xi_plus*gamma - delta**2*gamma)*r0 - 4*xi_plus*r + 2*delta*xi_plus*gamma
    main_1_deno = (4*xi_plus - delta**2)*sqrt(r0*(r0 - delta) + xi_plus)
    main_1 = main_1_nomi/main_1_deno
    main_2 = gamma*log(2*(sqrt(r0*(r0 - delta) + xi_plus) + r0) - delta)
    main_3_nomi = -4*xi_plus*r + 2*delta*xi_plus*gamma
    main_3_deno = (4*xi_plus - delta**2)*sqrt(xi_plus)
    main_3 = main_3_nomi/main_3_deno
    main_4 = gamma*log(2*sqrt(xi_plus) - delta)
    return main_1 - main_2 - main_3 + main_4

def B_plus_ellip_axial_radial(J, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_plus_ellip_axial_radial(h, a, b, r, z, alpha, theta), -pi, pi)[0]
    return J*integral_B/4/pi

def main_minus_ellip_axial_radial(h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    delta = 2*r*cos(alpha - theta)
    xi_minus = r**2 + z**2
    gamma = cos(alpha - theta)
    main_1_nomi = 2*(delta*r + 2*xi_minus*gamma - delta**2*gamma)*r0 - 4*xi_minus*r + 2*delta*xi_minus*gamma
    main_1_deno = (4*xi_minus - delta**2)*sqrt(r0*(r0 - delta) + xi_minus)
    main_1 = main_1_nomi/main_1_deno
    main_2 = gamma*log(2*(sqrt(r0*(r0 - delta) + xi_minus) + r0) - delta)
    main_3_nomi = -4*xi_minus*r + 2*delta*xi_minus*gamma
    main_3_deno = (4*xi_minus - delta**2)*sqrt(xi_minus)
    main_3 = main_3_nomi/main_3_deno
    main_4 = gamma*log(2*sqrt(xi_minus) - delta)
    return main_1 - main_2 - main_3 + main_4
    
def B_minus_ellip_axial_radial(J, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_minus_ellip_axial_radial(h, a, b, r, z, alpha, theta), -pi, pi)[0]
    return -J*integral_B/4/pi

def B_ellip_cylinder_axial_radial(J, a, b, h, r, z, alpha):
    return B_plus_ellip_axial_radial(J, a, b, h, r, z, alpha) + B_minus_ellip_axial_radial(J, a, b, h, r, z, alpha)


# In[ ]:





# In[ ]:





# In[6]:





# In[ ]:






# ## Computation of magnetic field from elliptical cylinder with diametrical magnetization

# In[12]:


# Importing relevant libraries

import matplotlib as plt
import math as m
from numpy import sqrt, sin, cos, pi
from scipy.integrate import quad, tplquad

"""
Parameter definition

The magnetization is diametrically oriented (V. T. Nguyen, T - F. Lu, "Modelling of magnetic field distributions of elliptical cylinder permanent
magnets with diametrical magnetization", Journal of magnetism and magnetic materials, vol. 491, 2019)

The function to compute the magnetic field is denoted as: 
B_ellip_cylinder_diam_axial (azimuthal, radial) where the meaning of the abbreviations as B - magnetic flux density; 
ellip_cylinder - elliptical cylinder; diam - diametrical orientation of the magnetization; 
axial (azimuthal, radial) - the axial (azimuthal, radial) component of the magnetic flux density

J: the remanence in Tesla of the cylinder
h: the height of the elliptical cylinder
R: the radius of the elliptical cylinder
(r, z, alpha): the coordinates of the computed point in the Cylindrical coordinate system
"""

# Defining the constants
mu0 = 4*pi*10**-7 # the magnetic permeability

"""Computing the axial component of the magnetic field B_ellip_cylinder_diam_axial"""
# B_ellip_cylinder_diam_axial = B_diam_plus + B_diam_minus

# Computation of B_plus
def main_ellip_diam_axial(J_X, J_Y, h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    di = sqrt(r0**2 + r**2 - 2*r0*r*cos(alpha - theta))
    main_1 = 1/sqrt((z - h)**2 + di**2)
    main_2 = 1/sqrt(z**2 + di**2)
    main_3_num = a*b**3*J_X*cos(theta) + b*a**3*J_Y*sin(theta)
    main_3_deno = (b**2*cos(theta)**2 + a**2*sin(theta)**2)**(3/2)
    main_3 = main_3_num/main_3_deno
    
    return (main_1 - main_2)*main_3

def B_ellip_cylinder_diam_axial(J_X, J_Y, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_ellip_diam_axial(J_X, J_Y, h, a, b, r, z, alpha, theta), 0, 2*pi)[0]
    return integral_B/4/pi 

"""Computing the azimuthal component of the magnetic field B_ellip_cylinder_diam_azimuthal"""

def main_ellip_diam_azimuthal(J_X, J_Y, h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    di = sqrt(r0**2 + r**2 - 2*r0*r*cos(alpha - theta))
    main_1 = (h - z)/(di**2*sqrt(di**2 + (h - z)**2))
    main_2 = z/(di**2*sqrt(di**2 + z**2))
    main_3 = r0*sin(alpha - theta)
    main_4_num = a*b**3*J_X*cos(theta) + b*a**3*J_Y*sin(theta)
    main_4_deno = (b**2*cos(theta)**2 + a**2*sin(theta)**2)**(3/2)
    main_4 = main_4_num/main_4_deno
    return (main_1 + main_2)*main_3*main_4

def B_ellip_cylinder_diam_azimuthal(J_X, J_Y, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_ellip_diam_azimuthal(J_X, J_Y, h, a, b, r, z, alpha, theta), 0, 2*pi)[0]
    return integral_B/4/pi


"""Computing the radial component of the magnetic field B_ellip_cylinder_diam_radial"""

def main_ellip_diam_radial(J_X, J_Y, h, a, b, r, z, alpha, theta):
    r0 = a*b/sqrt(a**2*sin(theta)**2 + b**2*cos(theta)**2)
    di = sqrt(r0**2 + r**2 - 2*r0*r*cos(alpha - theta))
    main_1 = (h - z)/(di**2*sqrt(di**2 + (h - z)**2))
    main_2 = z/(di**2*sqrt(di**2 + z**2))
    main_3 = (r - r0*cos(alpha - theta))
    main_4_num = a*b**3*J_X*cos(theta) + b*a**3*J_Y*sin(theta)
    main_4_deno = (b**2*cos(theta)**2 + a**2*sin(theta)**2)**(3/2)
    main_4 = main_4_num/main_4_deno
    return (main_1 + main_2)*main_3*main_4

def B_ellip_cylinder_diam_radial(J_X, J_Y, a, b, h, r, z, alpha):
    integral_B = quad(lambda theta: main_ellip_diam_radial(J_X, J_Y, h, a, b, r, z, alpha, theta), 0, 2*pi)[0]
    return integral_B/4/pi





# ## Computation of magnetic field from a sphere with axial (Z coordinate) magnetization

# In[1]:


# Importing relevant libraries

import matplotlib as plt
import math as m
from numpy import sqrt, sin, cos, pi
from scipy.integrate import quad, dblquad

"""
Parameter definition

The magnetization is axially oriented 

The function to compute the magnetic field is denoted as: 
H_sphere_axial_axial (azimuthal, radial) where the meaning of the abbreviations as H - magnetic 
field strength; 
sphere - spherical geometry; axial - axial orientation of the magnetization; 
axial, azimuthal and radial - the axial, azimuthal and radial components of the magnetic 
field strength

J: the remanence in Tesla of the sphere
R: the radius of the sphere
(r, z, phi): the coordinates of the computed point in the Cylindrical coordinate system
"""

# Defining the constants
mu0 = 4*pi*10**-7 # the magnetic permeability


# Computation of H_sphere_axial_axial
def main_sphere_axial_axial(J, R, r, z, phi, zO, alpha):
    rO = sqrt(R**2 - zO**2)
    main_num = (z - zO)*zO
    main_deno = (rO**2 + r**2 - 2*rO*r*cos(phi - alpha) + (z - zO)**2)**(3/2)
    return main_num/main_deno

def H_sphere_axial_axial(J, R, r, z, phi):
    integral_H = dblquad(lambda alpha, zO: main_sphere_axial_axial(J, R, r, z, phi, zO, alpha), 
                         -R, R, lambda zO: 0, lambda zO: 2*pi)[0]
    return J*integral_H/4/pi/mu0 



