# Written 2/12/15 by dh4gan
# Helper module for fractal density calculations

import numpy as np
from datetime import datetime
from scipy.fftpack import ifft, ifftn
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

    return fractalD, npower, rhozero, seed, ngrid, boxlength, filename

def readInputs_Cuboid():
    
    fractalD, npower, rhozero,seed = readInputs()
    
    # Grid size and dimensions
    
    ngrid = input("How many grid cells per dimension? ")
    xlength = input("x box half-length ([-L,L]): ")
    ylength = input("y box length: ")
    zlength = input("z box length: ")        
            
    # Output filename


    filename =  raw_input("What is the output filename? ")

    if filename=='':
        filename = 'grid.dat'
        print "No filename read: setting to default filename ",filename
    
    print "--"
    
    return fractalD,npower,rhozero,seed,ngrid,xlength,ylength,zlength,filename


def create1DGrid(ngrid,boxlength):
    dr = 2.0*boxlength/float(ngrid)
    
    print "Defining 1D Grid"
    print "Boxlength: ", boxlength
    print "Cell Size: ", dr
            
    x = np.zeros(ngrid)    
    
    for ix in range(ngrid):
        x[ix] = -boxlength + ix*dr
                    
    print '--'
    return x, dr 

    

def createCube(ngrid,boxlength):
    '''Create a cubic grid'''
    
    dr = 2.0*boxlength/float(ngrid)
    
    print "Defining Cubic Grid"
    print "Boxlength: ", boxlength
    print "Cell Size: ", dr
            
    x = np.zeros(ngrid)
    y = np.zeros(ngrid)
    z = np.zeros(ngrid)
    
    for ix in range(ngrid):
        x[ix] = -boxlength + ix*dr
        y[ix] = x[ix]
        z[ix] = x[ix]
                    
    print '--'
    return x,y,z,dr 

def createCuboid(ngridx,ngridy,ngridz,xlength,ylength,zlength):
    '''Create a cuboid grid '''
    
    dx = 2.0*xlength/float(ngridx)
    dy = 2.0*ylength/float(ngridy)
    dz = 2.0*zlength/float(ngridz)
    
    print "Defining Cuboid Grid"
    print "Boxlengths: ", xlength, ylength, zlength
    print "Cell Sizes: ", dx,dy,dz
    
    x = np.zeros(ngridx)
    y = np.zeros(ngridy)
    z = np.zeros(ngridz)
    
    for ix in range(ngridx):        
        x[ix] = -xlength + ix*dx
        
    for iy in range(ngridy):
        y[iy] = -ylength + iy*dy
        
    for iz in range(ngridz):
        z[iz] = -zlength + iz*dz
    
    print '--'
    return x,y,z,dx,dy,dz

def constructWavenumbers(ngrid,boxlength,dr):
    '''Construct wavenumbers so that they can be used in inverse FFT routines'''
    
    print "Constructing Wavenumbers: "
    print "Length of array: ",ngrid
    print "Maximum Scale: ", boxlength
    print "Minimum Scale: ", dr
    
    # Construct k such that it can be correctly fed into ifftn   
    # k[0] = zero frequency
    # k[1:ngrid/2+1] = positive frequencies (ascending)
    # k[ngrid/2+1:-1] = negative frequencies (decreasingly negative)

    
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
    
    print '--'
    return k

def computePowerLawPowerSpectrum(npower,ngridx,ngridy,ngridz, kx,ky,kz):
    '''Compute a power spectrum in 3D which is a simple power law'''
    
    
    print "Computing 3D Power Law Power Spectrum, index: ",npower
    
    powerspec = np.zeros((ngridx,ngridy,ngridz))

    for ix in range(1,ngridx):
        for iy in range(1,ngridy):
            for iz in range(1,ngridz):
            
                kmag= np.sqrt(kx[ix]*kx[ix]+ ky[iy]*ky[iy] + kz[iz]*kz[iz])
                powerspec[ix,iy,iz] = np.power(kmag,-npower)

    print '--'
    return powerspec

def PowerSpectrumToFT(seed, powerspec):
    '''Generates Fourier Transform of a function given its power spectrum
    Assigns random phases to each Fourier mode'''
    
    
    # The power spectrum has no phase data - must generate it before Fourier Transform
    # Use random phases to reconstruct the Fourier Transform of the density

    print "Generating Fourier Transform from Power Spectrum"
    print "Populating Random Fourier phases, with number seed ",seed

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
    
    print "Performing Inverse FFT"
    return ifft(ft)
    

def InverseFFTN(ft):
    '''Wrapper function to return Inverse FFT'''
    
    print "Performing Inverse multidimensional FFT"
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

def writeCubicGridToFile(x,y,z,ngrid,rho,dr,filename):
    
    print "--"
    print "Writing cubic grid to file ", filename
    f_obj = open(filename, 'w')

    # First write number of x, y and z cells to header

    line = str(ngrid)+' '+str(ngrid)+ ' ' +str(ngrid) 

    f_obj.write(line+'\n')

    # Now write in format
    # xcell ycell zcell dx dy dz rho  (x,y,z co-ordinates being initial cell face)

    for ix in range(ngrid):
        for iy in range(ngrid):
            for iz in range(ngrid):
                line = str(x[ix]) + ' ' + str(y[iy]) + ' ' + str(z[iz]) + ' ' + str(dr) + ' ' + str(dr)+ ' ' + str(dr) +' ' +str(rho[ix,iy,iz])
                f_obj.write(line+'\n')

    print "File Write Complete"
    
    