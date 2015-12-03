# Written 2/12/15 by Duncan Forgan
# Creates a uniform density grid (useful control to compare with fractal grids)

import numpy as np
import fractalmodule as fr

# Read in input data:

print "-----"
print "UNIFORM DENSITY DISTRIBUTION IN A 3D CUBE"
print "(control for fractal simulations)"
print "-----"

# fractalD - fractal dimension
# rhozero - density scaling parameter
# ngrid - Number of grid cells per dimension (i.e. total ngrid^3)
# boxlength - Half length of box in code units (e.g. dimensions (-L,L))
# filename - output filename

fractalD, npower, deltarho, rhozero, seed, distkey, distunit, masskey, massunit,ngrid, boxlength, filename = fr.readInputs_Cube()

# Define the cubic Grid
x,y,z,dr = fr.createCube(ngrid, boxlength,distkey)

print "***"
print "***This script produces a control cube of uniform density: ignoring fractal calculations***"
print "***"

rho = np.zeros((ngrid,ngrid,ngrid))

rho[:,:,:] = 1.0

# Write density grid to file

fr.writeCubicGridToFile(x,y,z,ngrid,rho,dr,massunit,distunit,masskey,distkey,filename)

