'''
June 10th, 2019
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
from scipy.interpolate import interp1d
import pandas as pd

Dir='/home/heh15/workingspace/Arp240/NGC5257/SM/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
picDir=Dir+'picture/'
logDir=Dir+'log/'


############################################################
# basic setting

PA= 110
incl=0.45
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

############################################################
# main program

### Arp 240
mstar=np.transpose(np.loadtxt(logDir+'stellarmass.txt'))
mH2=np.transpose(np.loadtxt(logDir+'H2mass.txt'))
fraction=np.power(10,mH2[1])/np.power(10,mstar[1])

### NGC 3198
N3198_vel_star=np.transpose(np.loadtxt(logDir+'NGC3198_star_vel.txt',delimiter=','))
N3198_mass_star=2.33*10**5*N3198_vel_star[0]*N3198_vel_star[1]**2
N3198_star_ring=N3198_mass_star[:-1]-N3198_mass_star[1:]
N3198_vel_H2=np.transpose(np.loadtxt(logDir+'NGC3198_H2_vel.txt',delimiter=','))
N3198_mass_H2=2.33*10**5*N3198_vel_H2[0]*N3198_vel_H2[1]**2
f_out=interp1d(N3198_vel_H2[0],N3198_mass_H2)
N3198_mass_H2=f_out(N3198_vel_star[0])
N3198_H2_ring=N3198_mass_H2[:-1]-N3198_mass_H2[1:]
fraction1=N3198_mass_H2/N3198_mass_star

### Other galaxies. 
galaxies=['NGC0925','NGC2403','NGC2841','NGC2903','NGC2976','NGC3521','NGC4736','NGC5055','NGC6946','NGC7331']
fractions=dict.fromkeys(galaxies)
star_rads=dict.fromkeys(galaxies)

# Usually the stellar velocity range is greater than the H2 range. 
# default star1 range smaller than star2
def data_import(galaxy):
    H2file=logDir+'Frank_galaxy/'+galaxy+'_vel_H2.csv'
    H2=pd.read_csv(H2file,header=None)
    H2=H2.sort_values(0,ascending=True)
    starfiles=glob.glob(logDir+'Frank_galaxy/'+galaxy+'_vel_star*.csv')
    star1=pd.read_csv(starfiles[0],header=None)
    star1=star1.sort_values(0,ascending=True)
    if np.shape(starfiles)[0]==1:
        star2=pd.DataFrame()
        star2[0]=np.copy(star1[0])
        star2[1]=np.zeros(np.shape(star1)[0])
    else: 
        star2=pd.read_csv(starfiles[1],header=None)
        star2=star2.sort_values(0,ascending=True)
    star=pd.DataFrame()
    star[0]=star1[0]
    f_out=interp1d(star2[0],star2[1])
    try:
        star2_vel=f_out(star1[0])
    except ValueError:
        low=0;up=-1
        if star1[0].iloc[0]<star2[0].iloc[0]:
            low=np.argwhere(star1[0]>star2[0][0])[0][0]
        if star1[0].iloc[-1]>star2[0].iloc[-1]:
            up=np.argwhere(star1[0]<star2[0].iloc[-1])[-1][0]
        else:
            print('check your input data')
        star[0]=star1[0].iloc[low:up]
        star2_vel=f_out(star[0])
    star[1]=np.sqrt(star2_vel**2+star1[1]**2)
    star=star.dropna()
    return star, H2


def fraction_calc(H2,star):
    mass_H2=2.33*10**5*H2[0]*H2[1]**2
    star_rad=H2[0]
    f_out=interp1d(star[0],star[1])
    try:
        star_vel=f_out(star_rad)
    except ValueError:
        low=0;up=-1
        if H2[0].iloc[0]<star[0].iloc[0]:
            low=np.argwhere(H2[0]>star[0].iloc[0])[0][0]
        if H2[0].iloc[-1]>star[0].iloc[-1]:
            up=np.argwhere(H2[0]<star[0].iloc[-1])[-1][0]
        else:
            print('check your input data')
        mass_H2=mass_H2.iloc[low:up]
        star_rad=H2[0].iloc[low:up]
        star_vel=f_out(star_rad)
    mass_star=2.33*10**5*star_rad*star_vel**2
    fraction=mass_H2/mass_star
    return star_rad,fraction

fig=plt.figure()
for galaxy in galaxies: 
    star,H2=data_import(galaxy)
    fractions[galaxy]=fraction_calc(H2,star)[1] 
    star_rads[galaxy]=fraction_calc(H2,star)[0]
    plt.plot(star_rads[galaxy],fractions[galaxy],label=galaxy)

plt.plot(mstar[0],fraction,label='NGC 5257',linestyle='--')
plt.plot(N3198_vel_star[0],fraction1,label='NGC 3198',linestyle='--')
plt.legend()
plt.savefig(picDir+'fraction_comparison.png')

