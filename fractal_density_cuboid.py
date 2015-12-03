# Written 16/11/15 by Duncan Forgan
# Given an input power spectrum, generates fractal density grids in 3D
# Taken from Walch et al, arXiv:1109.3478

import numpy as np
import fractalmodule as fr

# Read in input data:

print "-----"
print "FRACTAL DENSITY DISTRIBUTION IN A 3D CUBOID"
print "-----"

# fractalD - fractal dimension
# rhozero - density scaling parameter
# ngrid - Number of grid cells per dimension (i.e. total ngrid^3)
# boxlength - Half length of box in code units (e.g. dimensions (-L,L))
# filename - output filename

fractalD, npower, deltarho, rhozero, seed, distkey, distunit, masskey, massunit, ngrid, xlength, ylength, zlength, filename = fr.readInputs_Cuboid()

# Define the cuboid Grid
x,y,z,dx,dy,dz = fr.createCuboid(ngrid,ngrid,ngrid,xlength,ylength,zlength,distkey)

# Create Wavenumbers
kx = fr.constructWavenumbers(ngrid, xlength, dx,distkey)
ky = fr.constructWavenumbers(ngrid, ylength, dy,distkey)
kz = fr.constructWavenumbers(ngrid, zlength, dz,distkey)

# Compute a 3D power spectrum (power law k^-npower)
powerspec = fr.computePowerLawPowerSpectrum(npower, ngrid, ngrid, ngrid, kx, ky, kz)

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

fr.writeCuboidGridToFile(x,y,z,ngrid,ngrid,ngrid,rho,dx,dy,dz,massunit,distunit,filename)

