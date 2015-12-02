# Written 16/11/15 by Duncan Forgan
# Given an input power spectrum, generates fractal density grids in 3D
# Taken from Walch et al, arXiv:1109.3478

import numpy as np
import fluctuations as fl

# Read in input data:

# fractalD - fractal dimension
# rhozero - density scaling parameter
# ngrid - Number of grid cells per dimension (i.e. total ngrid^3)
# boxlength - Half length of box in code units (e.g. dimensions (-L,L))
# filename - output filename

fractalD, npower, rhozero, seed, ngrid, boxlength, filename = fl.readInputs_Cube()

# Define the cubic Grid
x,y,z,dr = fl.createCube(ngrid, boxlength)

# Create Wavenumbers
k = fl.constructWavenumbers(ngrid, boxlength, dr)

# Compute a 3D power spectrum (power law k^-npower)
powerspec = fl.computePowerLawPowerSpectrum(npower, ngrid, ngrid, ngrid, k, k, k)

# The power spectrum has no phase data - must generate it before Fourier Transform
# Use random phases to reconstruct the Fourier Transform of the density

rho_fft = fl.PowerSpectrumToFT(seed, powerspec)
            
# Compute the inverse fourier transform of rho_fft

rho = fl.InverseFFTN(rho_fft)

# Scale density by rhozero

print "Scaling Density"

rho = np.real(rho)
rho = np.real(np.exp(rho/rhozero))

# Test the density to check its PDF is lognormal (optional - feel free to comment out)
fl.testForLognormalPDF(rho, ngrid)

# Write density grid to file

fl.writeCubicGridToFile(x,y,z,ngrid,rho,dr,filename)

