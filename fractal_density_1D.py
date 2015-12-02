# Written 16/11/15 by Duncan Forgan
# Given an input power spectrum, generates fractal density grids in 1D
# Taken from Walch et al, arXiv:1109.3478

import numpy as np
import fractalmodule as fr

# Read in input data:

# fractalD - fractal dimension
# rhozero - density scaling parameter
# ngrid - Number of grid cells 
# boxlength - Half length of box in code units (e.g. dimensions (-L,L))
# filename - output filename

fractalD, npower, rhozero, seed, ngrid, boxlength, filename = fr.readInputs_Cube()

# Define the cubic Grid
x,dr = fr.create1DGrid(ngrid,boxlength)

# Create Wavenumbers
k = fr.constructWavenumbers(ngrid, boxlength, dr)

# Compute a 3D power spectrum (power law k^-npower)
powerspec = fr.compute1DPowerLawPowerSpectrum(npower, ngrid, k)

# The power spectrum has no phase data - must generate it before Fourier Transform
# Use random phases to reconstruct the Fourier Transform of the density

rho_fft = fr.PowerSpectrumToFT(seed, powerspec)
            
# Compute the inverse fourier transform of rho_fft

rho = fr.InverseFFT(rho_fft)

# Scale density by rhozero

print "Scaling Density"

rho = np.real(rho)
rho = np.real(np.exp(rho/rhozero))

# Test the density to check its PDF is lognormal (optional - feel free to comment out)
nbins = ngrid/10
fr.testForLognormalPDF(rho, nbins)


