# Written 2/12/15 by dh4gan
# Helper module for fractal density calculations

import numpy as np
from datetime import datetime
from scipy.fftpack import ifftn
from scipy.stats import norm
from matplotlib.mlab import normpdf
import matplotlib.pyplot as plt


def readInputs():
    '''Reads input data common to any fractal density calculation'''
    fractalD = input("What is the fractal dimension? ")

    npower = 2.0*(3.0-fractalD)+2

    print "Fractal Dimension is ", fractalD
    print "This corresponds to a power law power spectrum in k-space with index ",npower

    rhozero = input("What is rhozero? ")
    
    try:
        seed = input("Input a random number seed, or press enter to generate it from system clock: ")
    except SyntaxError:
        date = datetime.today()
        print "Generating Random Number seed from current timestamp ", date
        seed = np.int(date.hour + date.minute + date.second + date.day + date.month + date.year)
        print "Seed is ", seed
    
    return fractalD,npower,rhozero,seed

def readInputs_Cube():
        
    fractalD, npower, rhozero,seed = readInputs()

    # Grid size and dimensions

    ngrid = input("How many grid cells per dimension? ")
    boxlength = input("What is the box length (in pc)? ")
    
    # Output filename

    filename =  raw_input("What is the output filename? ")

    if filename=='':
        filename = 'grid.dat'
        print "No filename read: setting to default filename ",filename
    
    print "--"

    return fractalD, npower, rhozero, seed, ngrid, boxlength,filename


def readInputs_Cuboid():
    
    fractalD, npower, rhozero,seed = readInputs()
    
    # Grid size and dimensions
    
    ngrid = input("How many grid cells per dimension? ")
    xmin = input("Minimum x co-ordinate")
    xmax = input("Maximum x co-ordinate")
    
    ymin = input("Minimum y co-ordinate")
    ymax = input("Maximum y co-ordinate")
    
    zmin = input("Minimum z co-ordinate")
    zmax = input("Maximum z co-ordinate")
    
    # Output filename


    filename =  raw_input("What is the output filename? ")

    if filename=='':
        filename = 'grid.dat'
        print "No filename read: setting to default filename ",filename
    
    print "--"
    
    return fractalD,npower,rhozero,seed,ngrid,xmin,xmax,ymin,ymax,zmin,zmax,filename


def constructWavenumbers(ngrid,boxlength,dr):
    
    k = np.zeros(ngrid)
    k[0] = 0.0

    # positive frequencies first
    for ix in range(1,ngrid/2+1):
        #k[ix] = ix
        k[ix] = 2.0*np.pi/(2.0*boxlength-ix*dr)
        #print k[ix], boxlength-ix*dr
    
        # Now negative frequencies
    for ix in range(1,ngrid/2):
        k[ngrid/2+ix] = -k[ngrid/2-ix]
    
    return k


def PowerSpectrumToFT(seed, powerspec):
    '''Generates Fourier Transform of a function given its power spectrum
    Assigns random phases to each Fourier mode'''
    
    
    # The power spectrum has no phase data - must generate it before Fourier Transform
    # Use random phases to reconstruct the Fourier Transform of the density

    ft = np.sqrt(powerspec)

    randgen = np.random.RandomState(seed)
    j = np.complex(0,1.0) # imaginary number

    print "Computing Random Phases, seed: ",seed
    phases = j*2.0*np.pi*np.random.random(powerspec.shape)
    phases = np.exp(phases)

    ft = ft*phases
    
    return ft

def InverseFFT(ft):
    '''Wrapper function to return Inverse FFT'''
    return ifftn(ft)
    
    
def testForLognormalPDF(rho, ngrid):
    '''Tests a density grid's PDF for lognormality'''
    
    meanrho = np.mean(rho.flatten())
    rhotest = np.log(rho.flatten()/meanrho)

    (mu,sigma) = norm.fit(rhotest)

    # Test that the density histogram is lognormal

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    n, bins, patches= ax3.hist(rhotest, bins = ngrid)

    # Construct a lognormal PDF with the fitted parameters 
    rhofit = normpdf(bins, mu,sigma)

    # Scale the normal distribution to fit
    histmax = np.amax(n)
    fitmax = np.amax(rhofit)

    rhofit = histmax*rhofit/fitmax

    print "Density PDF is fitted by a lognormal with parameters: "
    print "Mean: ",mu
    print "Sigma: ", sigma

    ax3.plot(bins, rhofit)
    ax3.set_xlabel(r'log $\frac{\rho}{\bar{\rho}}$', fontsize = 16)
    ax3.set_ylabel(r'PDF', fontsize = 16)
    plt.show()

    
    