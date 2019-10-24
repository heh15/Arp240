'''
Oct. 25th, 2018

Calculate the total star formation rate through spitzer, herschel and 33 GHz image. 

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

############################################################
# basic settings

# general
Te=1e4

d=99

############################################################
# function

def sfr_uv(counts):
    flux=counts*1.4*10**(-15)*1528**2/(3*10**18)*0.68*10**(23)
    flux_erg=flux/10**23
    L_uv=4*math.pi*(d*10**6*3.086*10**18)**2*flux_erg
    SFR_tot=0.68*10**(-28)*L_uv
    
    return SFR_tot

def sfr_24um(f_mjsr):
    flux=f_mjsr*10**(-17)/(4.25*10**10)*2.45**2
    L_nu=4*math.pi*(100*10**6*3.086*10**18)**2*flux
    L=L_nu*1.2*10**13
    SFR_tot=10**(-42.69)*L
    
    return SFR_tot

def sfr_70um(f_her):
    flux_erg=f_her*10**(-23)
    L_nu=4*math.pi*(100*10**6*3.086*10**18)**2*flux_erg
    L=4.283*10**12*L_nu
    SFR_tot=10**(-43.23)*L

    return SFR_tot
    
############################################################
# main program

counts=21
SFR_UV=sfr_uv(counts)

f_mjsr=3675
SFR_24um=sfr_24um(f_mjsr)

f_her=5.68
SFR_70um=sfr_70um(f_her)


f_33GHz=4.14E-3
f_erg=f_33GHz*10**(-23)
L_nu=4*math.pi*(d*10**6*3.086*10**18)**2*f_erg
SFR_fre=4.6E-28*(Te/10**4)**(-0.45)*33**(0.1)*L_nu
SFR_tot=10**-27*(2.18*33**(-0.1)+15.1*33**(-0.7))**(-1)*L_nu
