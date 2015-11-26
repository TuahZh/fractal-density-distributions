# Written 16/11/15 by Duncan Forgan
# Given an input power spectrum, generates fractal density grids in 1D
# Taken from Walch et al, arXiv:1109.3478

import numpy as np
from scipy.fftpack import ifft
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib.mlab import normpdf
from datetime import datetime

# Read in input data:

# fractalD - fractal dimension
# rhozero - density scaling parameter
# ngrid - Number of grid cells per dimension (i.e. total ngrid^3)
# boxlength - Half length of box in code units (e.g. dimensions (-L,L))

fractalD = input("What is the fractal dimension? ")

npower = 2.0*(3.0-fractalD)+2

print "Fractal Dimension is ", fractalD
print "This corresponds to a power law power spectrum in k-space with index ",npower

rhozero = input("What is rhozero?")

# Grid size and dimensions

ngrid = input("How many grid cells are required?")
boxlength = input("What is the box length (in pc)? ")

try:
    seed = input("Input a random number seed, or press enter to generate it from system clock: ")
except SyntaxError:
    date = datetime.today()
    print "Generating Random Number seed from current timestamp ", date
    seed = np.int(date.hour + date.minute + date.second + date.day + date.month + date.year)
    print "Seed is ", seed
    

dr = 2.0*boxlength/float(ngrid)
ncells = ngrid*ngrid*ngrid
j = np.complex(0,1.0) # imaginary number

# Output filename

filename = raw_input("What is the output filename? ")


# Define the Grid
print "Defining Grid"
print "Boxlength: ", boxlength
print "Cell Size: ", dr

k = np.zeros(ngrid)
x = np.zeros(ngrid)

for ix in range(ngrid):
    x[ix] = - boxlength + ix*dr

print "Constructing Wavenumbers"

    
# Construct k such that it can be fed into ifft    
# k[0] = zero frequency
# k[1:ngrid/2+1] = positive frequencies (ascending)
# k[ngrid/2+1:-1] = negative frequencies (decreasingly negative)

k[0] = 0.0

# positive frequencies first
for ix in range(1,ngrid/2+1):
    k[ix] = 2.0*np.pi/(2.0*boxlength-ix*dr)
    
    
# Now negative frequencies
for ix in range(1,ngrid/2):
    k[ngrid/2+ix] = -k[ngrid/2-ix]
    
# Construct Power Spectrum

print '--'
print "Computing Power Spectrum"

powerspec = np.zeros(ngrid)

for ix in range(ngrid):
    kmag= np.sqrt(k[ix]*k[ix])      
    
    if(kmag==0.0): powerspec[ix]=0.0
    else: powerspec[ix] = np.power(kmag,-npower)

# The power spectrum has no phase data - must generate it before Fourier Transform
# Use random phases to reconstruct the Fourier Transform of the density

rho_fft = np.sqrt(powerspec)

randgen = np.random.RandomState(seed)

print "Computing Random Phases, seed: ",seed

phases = j*2.0*np.pi*np.random.random(ngrid)
phases = np.exp(phases)

rho_fft = rho_fft*phases
            
# Compute the inverse fourier transform of rho_fft

print "--"
print "Carrying out Inverse FFT"

rho = ifft(rho_fft)

# Scale by rhozero
print "Scaling Density"
rho = np.real(np.exp(rho/rhozero))

meanrho = np.mean(rho.flatten())
rhotest = np.log10(rho.flatten()/meanrho)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(x, rhotest)

(mu,sigma) = norm.fit(rhotest)

# Test that the density histogram is lognormal

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
n, bins, patches = ax3.hist(rhotest, bins = ngrid)

# Construct a lognormal PDF with the fitted parameters 
rhofit = normpdf(bins, mu,sigma)

# Scale the normal distribution to fit
histmax = np.amax(n)
fitmax = np.amax(rhofit)

rhofit = histmax*rhofit/fitmax

ax3.plot(bins, rhofit)
ax3.set_xlabel(r'log $\frac{\rho}{\bar{\rho}}$', fontsize = 16)
ax3.set_ylabel(r'PDF', fontsize = 16)
plt.show()

print "Density PDF is fitted by a lognormal with parameters: "
print "Mean: ",mu
print "Sigma: ", sigma

# Write density grid to file
print "--"
print "Writing to file ", filename

f_obj = open(filename, 'w')

# First write number of x cells to header

line = str(ngrid)

f_obj.write(line+'\n')

# Now write in format
# xcell rho  (x co-ordinates being leftward cell face)

for ix in range(ngrid):
    line = str(x[ix]) +' '+ str(rho[ix])
    f_obj.write(line+'\n')

f_obj.close()    

print "File Write Complete"      