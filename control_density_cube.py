# Written 2/12/15 by Duncan Forgan
# Creates a uniform density grid (useful control to compare with fractal grids)

import numpy as np
import fractalmodule as fr

# Read in input data:

# fractalD - fractal dimension
# rhozero - density scaling parameter
# ngrid - Number of grid cells per dimension (i.e. total ngrid^3)
# boxlength - Half length of box in code units (e.g. dimensions (-L,L))
# filename - output filename

ngrid = input("How many grid cells per dimension? ")
boxlength = input("What is the box length (i.e. [-L,L])? ")

filename =  raw_input("What is the output filename? ")

if filename=='':
    filename = 'grid.dat'
    print "No filename read: setting to default filename ",filename
    
print "--"

# Define the cubic Grid
x,y,z,dr = fr.createCube(ngrid, boxlength)

rho = np.zeros(ngrid,ngrid,ngrid)

rho[:,:,:] = 1.0

# Write density grid to file

fr.writeCubicGridToFile(x,y,z,ngrid,rho,dr,filename)

