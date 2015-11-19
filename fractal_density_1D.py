# Written 16/11/15 by Duncan Forgan
# Given an input power spectrum, generates fractal density grids in 3D

import numpy as np
from scipy.fftpack import fft, ifft, fftshift
import matplotlib.pyplot as plt

#fractalD = input("What is the fractal dimension? ")
fractalD = 2.4

npower = 2.0*(3.0-fractalD)+2

print "This corresponds to a power spectrum with index ",npower

#rhozero = input("What is rhozero?")
rhozero = 10.0

# Grid size and dimensions

#ngrid = input("How many grid cells are required?")
#boxlength = input("What is the box length (in pc)? ")

ngrid = 256
boxlength = np.pi

dr = 2.0*boxlength/float(ngrid)
ncells = ngrid*ngrid*ngrid
j = np.complex(0,1.0) # imaginary number

# Output filename

#filename = raw_input("What is the output filename? ")
filename = 'trial'

# Construct wavenumbers

k = np.zeros(ngrid) # kx, ky, kz identical for a given cell
x = np.zeros(ngrid)

for ix in range(ngrid):
    x[ix] = - boxlength + ix*dr
    
# Construct k such that it can be fed into ifft    
# k[0] = zero frequency
# k[1:ngrid/2+1] = positive frequencies (ascending)
# k[ngrid/2+1:-1] = negative frequencies (decreasingly negative)

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

powerspec = np.zeros(ngrid)

for ix in range(ngrid):
    kmag= np.sqrt(k[ix]*k[ix])      
    
    if(kmag==0.0): powerspec[ix]=0.0
    else: powerspec[ix] = np.power(kmag,-npower)

        
# Use random phases to reconstruct the Fourier Transform of the density
# Need to be careful about use of complex numbers here

rho_fft = np.sqrt(powerspec)

print "FFT: ",fft(rho_fft)

phases = j*2.0*np.pi*np.random.random(ngrid)
phases = np.exp(phases)

rho_fft = rho_fft*phases
            
# Compute the inverse fourier transform of rho_fft
# TODO - what is the appropriate 3D function?
rho = ifft(rho_fft)

print "x: ", x
print "rho: ", rho
# Scale density by rhozero

rho = np.exp(rho/rhozero)

print "rho scaled: ", rho

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.scatter(fftshift(k), fftshift(powerspec))
ax1.scatter(k, fft(rho_fft), color = 'green')
#ax1.plot(x,rho)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(x, rho)
plt.show()

# Write density grid to file

f_obj = open(filename, 'w')

# First write number of x, y and z cells to header

line = str(ngrid)+' '+str(ngrid)+ ' ' +str(ngrid) 

f_obj.write(line+'\n')



# Now write in format
# xcell rho  (x,y,z co-ordinates being leftward cell face)

for ix in range(ngrid):
    line = str(x[ix]) +' '+ str(rho[ix])
    f_obj.write(line+'\n')

f_obj.close()        