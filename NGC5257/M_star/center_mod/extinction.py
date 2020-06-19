'''
Aug. 31st, 2017
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

beammaj=2.021
beammin=1.610
beamarea=beammaj*beammin*1.1331
beamarea_pix=beamarea/0.09
freq=115.27
alpha=4.86/5
pc_cm=3.0857*10**18

CO12_peak=12.94
h2=2*1.661*10**(-24)
beta=-1.75

############################################################
# function

def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
    position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')
    size=u.Quantity((54,42),u.arcsec)
    cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs
    return wcs_cut, data_cut

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

# def mass_calc(flux_36,flux_45,uncertainty_36,uncertainty_45):
#     mass=math.pow(10,5.65)*np.power(flux_36,2.85)*np.power(flux_45,-1.85)*(D/0.05)**2*0.7
#     error_36=np.power(2.85*np.divide(uncertainty_36,flux_36),2)
#     error_45=np.power(1.85*np.divide(uncertainty_45,flux_45),2)
#     error_total=np.sqrt(error_36+error_45)
#     error= error_total*mass

#     return mass, error

def mass_calc(flux_36,flux_45):
    mass=math.pow(10,5.65)*np.power(flux_36,2.85)*np.power(flux_45,-1.85)*(D/0.05)**2*0.7

    return mass

def mass_sum(mass_ring):
    mass_sum=np.empty(shape=(0,0))
    for i in range(1,mass_ring.shape[0]+1,1):
        mass=mass_ring[range(i)].sum()
        mass_sum=np.append(mass_sum,mass)
    
    return mass_sum

def exp_func(x, a, b):
    return a*np.exp(-b * x)

def SB_std(mass_SB,means_36,means_45,stds_36,stds_45,counts):
    # SB_36=flux_36/counts
    # SB_45=flux_45/counts
    error_36=np.power(2.85*np.divide(stds_36,means_36),2)
    error_45=np.power(1.85*np.divide(stds_45,means_45),2)
    error_total=np.sqrt(error_36+error_45)
    uncertainty= mass_SB*error_total
    
    return uncertainty

############################################################
# main program

# get the extinction parameter in the center

peak_K=CO12_peak/(0.0109*beammaj*beammin*(freq/115.27)**2)
H2_mass=alpha*peak_K*2*10**33/(pc_cm**2)
NH2=H2_mass/h2
A_08=2.96*10**(-22)*NH2
A_34=A_08*(3.4/0.802)**beta
tau_34=-np.log(10**(A_34/(-2.5)))


# measure the flux of spitzer 3.6um image

Image='spitzer_36um_regrid.fits'
wcs_cut=fits_import(Image)[0]
data_cut=fits_import(Image)[1]
shape=data_cut.shape
flux=np.empty(shape=(0,0))
counts=np.empty(shape=(0,0))

# draw the center of the image
a=radius_arcsec[0]
b=radius_arcsec[0]*incl
center_sky=SkyEllipticalAperture(position,a=a*u.arcsec,b=b*u.arcsec,theta=PA*u.degree)
center_pix=center_sky.to_pixel(wcs=wcs_cut)
center_mask=Apmask_convert(center_pix,data_cut)


# draw the ring of the image
for i in range(8):
    rings[i]=aperture_ring(radius_arcsec[i+1],wcs_cut)
    rings_mask[i]=Apmask_convert(rings[i],data_cut)



