'''
Sept 5th, 2017
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
ra=204.97051
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

############################################################
# function

def fits_import(fitsimage, item=0):
    hdr = fits.open(fitsimage)[item].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[item].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

# example position and size format. 
def cut_2d(data_masked,position,size,wcs):
    cut=Cutout2D(data=data_masked,position=position,size=size,wcs=wcs)
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
    apmask=aperture.to_mask(method='center')
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

def flux_aperture_get(data_masked,aperture,rms,chans,chan_width, beamarea_pix):
    data_cut=data_masked.data
    mask=data_masked.mask
    if np.shape(chans) == ():
        chans = np.full(data_cut.shape, chans) 
    flux=aperture_photometry(data_cut,apertures=aperture,mask=mask)['aperture_sum'][0]/beamarea_pix
    error=np.sqrt(chans)*rms*chan_width/np.sqrt(beamarea_pix)
    uncertainty=aperture_photometry(data_cut,apertures=aperture,mask=mask,error=error)['aperture_sum_err'][0]

    return flux, uncertainty

def mass_calc(flux,uncertainty, Xco=2.0):
    mass=1.05*10**4*(Xco/2)*flux*D**2*1.4
    mass_error=mass*uncertainty/flux

    return mass, mass_error

def mass_sum(mass_ring):
    mass_sum=np.empty(shape=(0,0))
    for i in range(1,mass_ring.shape[0]+1,1):
        mass=mass_ring[range(i)].sum()
        mass_sum=np.append(mass_sum,mass)
    
    return mass_sum

def exp_func(x, a, b):
    return a*np.exp(-b * x)

############################################################
# main program

Image='NGC5257_12CO10_combine_contsub_pbcor_mom0.fits'
wcs, data_masked=fits_import(Image)
data = data_masked.data

position = position; cutshape = u.Quantity((54, 42), u.arcsec)
wcs_cut, data_cut = cut_2d(data_masked, position, cutshape, wcs)
shape=data_masked.shape
flux=np.empty(shape=(0,0))
counts=np.empty(shape=(0,0))
uncertainties=np.empty(shape=(0,0))

# draw the center of the image
a=radius_arcsec[0]
b=radius_arcsec[0]*incl
center_sky=SkyEllipticalAperture(position,a=a*u.arcsec,b=b*u.arcsec,theta=PA*u.degree)
center_pix=center_sky.to_pixel(wcs=wcs_cut)
center_mask=Apmask_convert(center_pix,data_cut)

# draw the ring of the image
for i in range(8):
    rings[i]=aperture_ring(radius_arcsec[i+1],wcs)
    rings_mask[i]=Apmask_convert(rings[i],data)

# fig=plt.figure()
# ax=plt.subplot(projection=wcs)
# ax.imshow(data,origin='lower')
# # center_pix.plot(color='red')
# # for i in range(size):
# #     rings[i].plot(color='red')
# # plt.show() 
# rings[3].plot(color='red')
# fig.savefig('../picture/turnover_h2.png')

# measure the flux within each ring. 
flux_center= aperture_photometry(data_cut,apertures=center_pix)['aperture_sum'][0]
count=center_mask.count()
flux=np.append(flux,flux_center)
counts=np.append(counts,count)
for i in range(size):
    flux_ring=aperture_photometry(data_cut,apertures=rings[i])['aperture_sum'][0]
    flux=np.append(flux,flux_ring)
    count=rings_mask[i].count()
    counts=np.append(counts,count)

# measure the uncertainty within each ring. 
rms=0.0016
uncertainty=flux_aperture_get(data_cut,center_pix,rms, 50, 10, beamarea_pix)[1]
uncertainties=np.append(uncertainties, uncertainty)
for i in range(len(rings)):
    uncertainty= flux_aperture_get(data_cut,rings[i],rms, 50, 10, beamarea_pix)[1]
    uncertainties=np.append(uncertainties, uncertainty)

flux=flux/beamarea_pix
mass_ring, mass_error=mass_calc(flux,uncertainties)
mass_radius=np.log10(mass_sum(mass_ring))

# write the output to the directory
output=np.transpose(np.vstack([radius_kpc,mass_radius]))
fmt= '%10.3e %10.3e \n'
filename='../log/H2mass.txt'
np.savetxt(filename,output,fmt=fmt,delimiter='',newline='\n')

# write the molecular mass and error within each ring
output=np.transpose(np.vstack([radius_kpc,mass_ring,mass_error]))
fmt= '%10.3e %10.3e %10.3e \n'
filename='../log/H2mass_ring.txt'
np.savetxt(filename,output,fmt=fmt,delimiter='',newline='\n')

# areas of the rings
areas=counts*0.09*500**2
SB=mass_ring/areas

# fig=plt.figure()
# plt.plot(radius_kpc,SB)

# SB1=SB[2:]
# radius_ksub=radius_kpc[2:]
# popt, pcov = curve_fit(exp_func,radius_ksub,SB1,p0=(1, 0.1))
# SB_mod=exp_func(radius_kpc,*popt)  
# plt.plot(radius_kpc,SB_mod)


## measure the total H2 mass. 

radius = 17.53 # SDSS petrosian radius for NGC 5257 in arcsec. 
center = position 
aperture_whole = SkyCircularAperture(center, radius*u.arcsec)
aperture_whole_pix = aperture_whole.to_pixel(wcs=wcs)

chans = 50
chan_width = 10 
rms = 1.6e-3
flux, uncertainty = flux_aperture_get(data_masked, aperture_whole_pix, rms, chans, chan_width, beamarea_pix)
mass, uncertainty = mass_calc(flux, uncertainty, Xco=0.5)
