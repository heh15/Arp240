'''
Oct.30, 2018

draw the photometry points for the continuum
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

############################################################
# basic settings
Dir='/home/heh15/workingspace/Arp240/NGC5257/continuum/'
scriptDir=Dir+'script/'
logDir=Dir+'log/'
picDir=Dir+'picture/'

############################################################
# main program

filename=logDir+'south_flux.txt'
skiprows=[0,1,5,6,7]
data=pd.read_csv(filename,header=None,sep=r"\s*",skiprows=skiprows)
data_south=data[:][0:3]

frequency=np.array([107,95,33])
freq_log=np.log10(10**9*frequency)
flux_log=np.log10(data_south[1])
error=0.4343*flux_log*data_south[2]/data_south[1]
fig=plt.figure()
south=plt.errorbar(freq_log,flux_log,error,marker='o',color='blue',label='south')

data_center=data[:][3:]
flux_log=np.log10(data_center[1])
flux_log=np.log10(data_center[1])
error=0.4343*flux_log*data_center[2]/data_center[1]
center=plt.errorbar(freq_log,flux_log,error,marker='o',color='green',label='center')
plt.legend()
plt.xlabel('log10(frequency) (Hz)')
plt.ylabel('log10(flux) (Jy)')
