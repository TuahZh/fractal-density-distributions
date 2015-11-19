# Written 16/11/15 by Duncan Forgan
# Given an input power spectrum, generates fractal density grids in 3D

import numpy as np
from scipy.fftpack import ifftn

# Read in fractal dimension


#fractalD = input("What is the fractal dimension? ")
fractalD = 2.4

npower = 2.0*(3.0-fractalD)+2

print "This corresponds to a power spectrum with index ",npower

#rhozero = input("What is rhozero?")
rhozero = 1.0

# Grid size and dimensions

#ngrid = input("How many grid cells are required?")
#boxlength = input("What is the box length (in pc)? ")

ngrid = 8
boxlength = 1.0

dr = boxlength/float(ngrid)
ncells = ngrid*ngrid*ngrid

# Output filename

#filename = raw_input("What is the output filename? ")
filename = 'trial'

# Construct wavenumbers

k = np.zeros(ngrid+1) # kx, ky, kz identical for a given cell


for ix in range(1,ngrid):
    k[ix]  = -np.pi + 2.0*np.pi/(ix*dr)


# Construct Power Spectrum

powerspec = np.zeros((ngrid,ngrid,ngrid))

for ix in range(1,ngrid):
    for iy in range(1,ngrid):
        for iz in range(1,ngrid):
            kmag= np.sqrt(k[ix]*k[ix]+ k[iy]*k[iy] + k[iz]*k[iz])
            
            powerspec[ix,iy,iz] = np.power(kmag,-npower)
            

# Use random phases to reconstruct the Fourier Transform of the density
# Need to be careful about use of complex numbers here

rho_fft = np.sqrt(powerspec)
j = np.complex(0,1.0) # imaginary number

phases = j*2.0*np.pi*np.random.random((ngrid,ngrid,ngrid))
phases = np.exp(phases)

#rho_fft = rho_fft*phases
            
# Compute the inverse fourier transform of rho_fft
# TODO - what is the appropriate 3D function?
rho = ifftn(rho_fft)

print rho.shape, rho[1,1,1], rho_fft[1,1,1]

# Scale density by rhozero

rho = np.exp(rho/rhozero)
print rho[1,1,1]

# Write density grid to file

f_obj = open(filename, 'w')

# First write number of x, y and z cells to header

line = str(ngrid)+' '+str(ngrid)+ ' ' +str(ngrid) 

f_obj.write(line+'\n')

# Now write in format
# xcell ycell zcell rho  (x,y,z co-ordinates being leftward cell face)

for ix in range(ngrid):
    x = ix*dr
    for iy in range(ngrid):
        y = iy*dr
        for iz in range(ngrid):
            z = iz*dr
            line = str(x) + ' ' + str(y) + ' ' + str(z) + ' ' + str(rho[ix,iy,iz])
            f_obj.write(line+'\n')

            