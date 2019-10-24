'''
Feb. 4th, 2019

To calculate the optical depth of the 12CO and 13CO

'''

import numpy as np
import os,glob
import scipy.ndimage as sni
import sys
import re
import itertools
from shutil import copytree
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import make_lupton_rgb
import aplpy
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from matplotlib.patches import Arrow
from photutils import SkyCircularAnnulus
from photutils import SkyCircularAperture
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import SkyEllipticalAperture
from photutils import SkyEllipticalAnnulus
from photutils import EllipticalAperture
from photutils import EllipticalAnnulus
from photutils import aperture_photometry
import numpy.ma as ma
import math
import shutil
import pandas as pd


Dir='/1/home/heh15/workingspace/Arp240/NGC5257/ratio/'
ratioDir=Dir+'1213/contsub/'
ratiomap=ratioDir+'NGC5257_1213_ratio.fits'
imageDir=Dir+'image/'

############################################################
# basic parameters
beammaj=2.186
beammin=1.896
freq=112.84

############################################################
# function
def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data[0][0]
    position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')
    size=u.Quantity((54,42),u.arcsec)
    cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs
    data_masked=np.ma.masked_invalid(data_cut)
    return wcs_cut, data_masked

def Jy2K(jansky,beammaj,beammin,freq):
    kelvin=jansky/(0.0109*beammaj*beammin*(freq/115.271)**2)

    return kelvin

def colDen(opdepth,j):
    h=6.62e-34;k=1.38e-23;mu=0.122*3.16e-25
    gj=2*j+1
    colDen=opdepth*3*h*gj/(8*(math.pi)**3*mu**2*j*(math.exp(1)-1))

    return colDen

############################################################
# main program

ratio=12

fitsimage=imageDir+'NGC5257_13CO10_12m_contsub_smooth_2rms_mom0.fits'
data_masked=fits_import(fitsimage)[1]

tau_13=1.0/ratio
tau_12=tau_13*50.0

colDen=colDen(tau_12,1)
