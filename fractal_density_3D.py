# Written 16/11/15 by Duncan Forgan
# Given an input power spectrum, generates fractal density grids in 3D

import numpy as np
from scipy.fftpack import ifftn

# Read in fractal dimension


#fractalD = input("What is the fractal dimension? ")
fractalD = 2.4

npower = 2.0*(3.0-fractalD)+2

print "Fractal Dimension is ", fractalD
print "This corresponds to a power law power spectrum in k-space with index ",npower

#rhozero = input("What is rhozero?")
rhozero = 1.0

# Grid size and dimensions

#ngrid = input("How many grid cells are required?")
#boxlength = input("What is the box length (in pc)? ")

ngrid = 16
boxlength = 1.0

dr = 2.0*boxlength/float(ngrid)
ncells = ngrid*ngrid*ngrid
j = np.complex(0,1.0) # imaginary number

# Output filename

#filename = raw_input("What is the output filename? ")
filename = 'trial'

# Construct wavenumbers

k = np.zeros(ngrid) # kx, ky, kz identical for a given cell
x = np.zeros(ngrid)
y = np.zeros(ngrid)
z = np.zeros(ngrid)

print "Defining Grid"
print "Boxlength: ", boxlength
print "Cell Size: ", dr

for ix in range(ngrid):
    x[ix] = - boxlength + ix*dr
    y[ix] = x[ix]
    z[ix] = x[ix]
    
# Construct k such that it can be fed into ifft    
# k[0] = zero frequency
# k[1:ngrid/2+1] = positive frequencies (ascending)
# k[ngrid/2+1:-1] = negative frequencies (decreasingly negative)

print "Constructing Wavenumbers"
k[0] = 0.0

# positive frequencies first
for ix in range(1,ngrid/2+1):
    #k[ix] = ix
    k[ix] = 2.0*np.pi/(2.0*boxlength-ix*dr)
    #print k[ix], boxlength-ix*dr
    
# Now negative frequencies
for ix in range(1,ngrid/2):
    k[ngrid/2+ix] = -k[ngrid/2-ix]

# Construct Power Spectrum

print '--'
print "Computing Power Spectrum"

powerspec = np.zeros((ngrid,ngrid,ngrid))

for ix in range(1,ngrid):
    for iy in range(1,ngrid):
        for iz in range(1,ngrid):
            kmag= np.sqrt(k[ix]*k[ix]+ k[iy]*k[iy] + k[iz]*k[iz])
            
            powerspec[ix,iy,iz] = np.power(kmag,-npower)
            

# Use random phases to reconstruct the Fourier Transform of the density
# Need to be careful about use of complex numbers here

rho_fft = np.sqrt(powerspec)

print "Computing Random Phases"
phases = j*2.0*np.pi*np.random.random((ngrid,ngrid,ngrid))
phases = np.exp(phases)

rho_fft = rho_fft*phases
            
# Compute the inverse fourier transform of rho_fft

print "--"
print "Carrying out Inverse FFT"
rho = ifftn(rho_fft)


# Scale density by rhozero

print "Scaling Density"
rho = np.real(np.exp(rho/rhozero))

# Write density grid to file

print "--"
print "Writing to file ", filename
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

print "File Write Complete"  