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
# directory

Dir='/home/heh15/workingspace/Arp240/NGC5258/SM/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
pictureDir=Dir+'picture/'
logDir=Dir+'log/'
fitsDir=Dir+'Fits/'

############################################################
# basic settings
D=99
sr_arcsec=(180/math.pi*60**2)**2
arcsec_pc=480

PA= 213.3
incl=0.43
xcenter=159.09
ycenter=158.29
# corresponds to  13h39m57.692s, 0d49'50.838"
ra=13*15+39*15.0/60.0+57.692*15.0/3600.0
dec=49.0/60.0+50.838/3600.0
steps=(33.75-0.75)/1.5+1
radius_arcsec=np.linspace(0.75,33.75,steps)
radius_kpc=radius_arcsec*0.48
size=radius_arcsec.shape[0]-1
position=SkyCoord(dec=dec*u.degree,ra=ra*u.degree,frame='icrs')
pixel_area=0.3*0.3
pixel_sr=pixel_area/(60**2*180/math.pi)**2
D=99

############################################################
# functions

def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial

    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def aperture_ring(radius_arcsec,wcs):
    a_in=radius_arcsec-1.5
    a_out=radius_arcsec
    b_out=a_out*incl
    ring_sky=SkyEllipticalAnnulus(position,a_in=a_in*u.arcsec,a_out=a_out*u.arcsec,b_out=b_out*u.arcsec,theta=PA*u.degree)
    ring_pix=ring_sky.to_pixel(wcs=wcs)
    return ring_pix

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask(method='center')[0]
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

############################################################
# main program
rings=dict.fromkeys((range(size)))
rings_mask=dict.fromkeys((range(size)))
SB=np.zeros(size+1)

fitsimage=fitsDir+'mass_map_masked.fits'
wcs=fits_import(fitsimage)[0]
data_masked=fits_import(fitsimage)[1]
data=data_masked.data

# draw the center of the image
a=radius_arcsec[0]
b=radius_arcsec[0]*incl
center_sky=SkyEllipticalAperture(position,a=a*u.arcsec,b=b*u.arcsec,theta=PA*u.degree)
center_pix=center_sky.to_pixel(wcs=wcs)
center_mask=Apmask_convert(center_pix,data)
SB[0]=np.ma.mean(center_mask)

# draw the rings of the image
for i in range(int(steps-1)):
    rings[i]=aperture_ring(radius_arcsec[i+1],wcs)
    rings_mask[i]=Apmask_convert(rings[i],data)
    SB[i+1]=np.ma.mean(rings_mask[i])

# fig=plt.figure()
# ax=plt.subplot(projection=wcs)
# ax.imshow(data,origin='lower')
# center_pix.plot(color='red')
# for i in range(int(steps-1)):
#     rings[i].plot(color='red')


# fig=plt.figure()
# plt.plot(radius_kpc,SB)

SB_crit=SB[0]/math.exp(1)
index=np.argmax(SB<SB_crit)
length=radius_kpc[4]
print(length)
