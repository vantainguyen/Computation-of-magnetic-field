This is an open-access module based on the license 4.0 which is free to use, adopt, change but not commercial use (Please contact the author for a permission). It is also a mandatory to cite the source of the equations based on which the codes are written as follows:
Van Tai Nguyen "Modelling of magnetic fields of permanent magnets with diametrical magnetization", MPhiL Thesis, the University of Adelaide, 2019.

# Computation-of-magnetic-field of a diamtrically magnetized cylinder
This source code is for the magnetic field computation of a diametrically magnetized cylinder
# The codes written herein are based on the equations derived by the author in his Master of Philosophy's Thesis (The University of Adelaide). The URL of the Thesis is attached as follows: https://hekyll.services.adelaide.edu.au/dspace/bitstream/2440/119972/1/Nguyen2019_Ma.pdf 
The parameters of the cylinder and the computed point in free space is as follows:
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
    of the magnetic field """  
    
    To use this Python module, simply download it and place it in the folder where is visible for the Python.
    # To call the module
    from Diametrical_cylinder import B_flux
    B_flux.Axial # To yield the axial component of the magnetic flux density
    B_flux.Azimuthal # To yield the azimuthal component of the magnetic flux density
    B_flux.Radial # To yield the radial component of the magnetic flux density
    
    
