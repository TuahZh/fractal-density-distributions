Fractal Density Distributions on Cartesian Grids
================================================

This is a set of python scripts designed to create cartesian density grids
(in 1D or 3D) with fractal properties.  It was originally designed to be used
in hydrodynamic simulations, in particular as a template to generate initial 
conditions for smoothed particle hydrodynamics (SPH) simulations.  A good place
to start reading on this would be [here](http://arXiv.org/abs/1109.3478)

The code uses the fact that the density power spectrum in a fractal space
with dimension D is given by 

P(k) ~ k^{-n}

Where the spectral index (Stutzki et al 2008, Federrath et al 2009)

n = 2(3-D)+2

The code initialises a density grid, and computes the above power spectrum. It
then adds random phase data to create a set of Fourier modes for the density, and
finally computes a Fourier Transform to create a normalised density distribution 
in configuration space (rho_FFT).  This density distribution is then scaled by the mean
density to give the true density

rho = rhozero exp(deltarho rho_FFT)

The deltarho variable is an input parameter which can increase or decrease the magnitude
of fluctuations in the density distribution.  Given a set of units as input, this 
density distribution is then written to an ASCII file.

Input parameters:
----------------

`fractalD` - fractal dimension

`deltarho` - strength of density fluctuations

`rhozero` - Mean density in the box

`ngrid` - Number of grid cells per dimension (i.e. total ngrid^3)

`boxlength` - Half length of box in code units (e.g. dimensions (-L,L))

`seed` - random number seed

`masskey, massunit` - string describing mass unit, variable for its value

`distkey, distunit` - string describing distance unit, variable for its value

`filename` - output filename

Code Dependencies:
------------------

This repository contains a module which handles most of the I/O and calculation. It
relies on the following Python modules:

numpy - numpy arrays, sqrt, int and real casting, complex numbers

datetime - to generate random number seeds from system clock

scipy.fftpack - inverse fourier transforms ifft, ifftn

These final three modules are used to check the resulting density distribution
has a probability density function (PDF) that is lognormal

`scipy.stats` - to fit the density PDF with a lognormal

`matplotlib.pyplot` - to construct and plot the density PDF

`matplotlib.mlab` - to plot the lognormal fit

