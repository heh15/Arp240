'''
Nov 23rd,2018

Gas mass estimated from the 12CO10 map. 
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
from scipy.optimize import curve_fit

############################################################
# basic setting

PA= 110
incl=0.7072
xcenter=159.71
ycenter=162.11
ra=204.97066052
dec=50.406569999999995
radius_arcsec=np.array([0.75, 2.25, 3.75, 5.25, 6.75,8.25,9.75,11.25,12.75])
radius_kpc=radius_arcsec*0.48
size=radius_arcsec.shape[0]-1
position=SkyCoord(dec=dec*u.arcmin,ra=ra*u.degree,frame='icrs')
rings=dict.fromkeys((range(size)))
rings_mask=dict.fromkeys((range(size)))
pixel_area=0.3*0.3
pixel_sr=pixel_area/(60**2*180/math.pi)**2
D=99
majorbeam=2.021
minorbeam=1.610
beamarea=majorbeam*minorbeam*1.1331
beamarea_pix=beamarea/0.09


rms=1.6e-3

############################################################
# function

def mass_calc(flux,uncertainty):
    mass=1.05*10**4*(0.4/2)*flux*D**2
    mass_error=mass*uncertainty/flux

    return mass, mass_error

############################################################
Npixs=363
uncertainty=rms*math.sqrt(50)*10*sqrt(Npixs/beamarea_pix)
flux=3.83

Gas_mass=mass_calc(flux,uncertainty)[0]
Gas_error=mass_calc(flux,uncertainty)[1]

# star formation rate from herschel. 
time_dep=Gas_mass/2.059
time_error=time_dep*Gas_error/Gas_mass
