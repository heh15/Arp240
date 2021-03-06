import cube
from astropy.utils import data
from spectral_cube import SpectralCube, Projection
from radio_beam import Beam
from astropy import units as u
import time
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits
from astropy.wcs import WCS

from matplotlib import rcParams
rcParams['mathtext.default']='regular'

############################################################
# directory

Dir='/1/home/heh15/workingspace/Arp240/scatter/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
picDir=Dir+'picture/'

############################################################

# basic setting
G=6.67e-11
pc=3.1e16
k=1.38e-23

alpha=0.86
ratio=0.77

incl=0.45
beammaj=1.1
beammin=0.8
freq=230.54

bin=10
threshold_33GHz=4*1.2e-5

############################################################
# function

def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

############################################################
# main program

Galaxies=['NGC5257','NGC5258']
