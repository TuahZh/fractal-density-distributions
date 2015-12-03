# Written 16/11/15 by Duncan Forgan
# Given an input power spectrum, generates fractal density grids in 3D
# Taken from Walch et al, arXiv:1109.3478

import numpy as np
import fractalmodule as fr

# Read in input data:

print "-----"
print "FRACTAL DENSITY DISTRIBUTION IN A 3D CUBE"
print "-----"


# fractalD - fractal dimension
# deltarho - strength of density fluctuations
# rhozero - Mean density in the box
# ngrid - Number of grid cells per dimension (i.e. total ngrid^3)
# boxlength - Half length of box in code units (e.g. dimensions (-L,L))
# seed - random number seed
# masskey, massunit - string describing mass unit, value
# distkey, distunit - string describing distance unit, value
# filename - output filename

fractalD, npower, deltarho, rhozero, seed, distkey, distunit, masskey, massunit,ngrid, boxlength, filename = fr.readInputs_Cube()

# Define the cubic Grid
x,y,z,dr = fr.createCube(ngrid, boxlength,distkey)

# Create Wavenumbers
k = fr.constructWavenumbers(ngrid, boxlength, dr, distkey)

# Compute a 3D power spectrum (power law k^-npower)
powerspec = fr.computePowerLawPowerSpectrum(npower, ngrid, ngrid, ngrid, k, k, k)

# The power spectrum has no phase data - must generate it before Fourier Transform
# Use random phases to reconstruct the Fourier Transform of the density

rho_fft = fr.PowerSpectrumToFT(seed, powerspec)
            
# Compute the inverse fourier transform of rho_fft

rho = fr.InverseFFTN(rho_fft)

# Scale density by rhozero

print "Scaling Density"

rho = np.real(rho)
rho = rhozero*np.real(np.exp(rho*deltarho))

# Test the density to check its PDF is lognormal (optional - feel free to comment out)
fr.testForLognormalPDF(rho, ngrid)

# Write density grid to file

fr.writeCubicGridToFile(x,y,z,ngrid,rho,dr,massunit,distunit,filename)

