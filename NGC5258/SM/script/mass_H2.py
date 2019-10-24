'''
Sept 5th, 2017

Target directory: $workingDir+H2mass
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
# directory
Dir='/home/heh15/workingspace/Arp240/NGC5258/SM/'
logDir=Dir+'log/'
workDir=Dir+'H2mass/'
imageDir=Dir+'image/'

os.chdir(workDir)

############################################################
# basic setting

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

def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data[0][0]
    position=SkyCoord(dec=dec*u.degree,ra=ra*u.degree,frame='icrs')
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

def flux_aperture_get(data_cut,aperture,rms):
    flux=aperture_photometry(data_cut,apertures=aperture)['aperture_sum'][0]
    error_value=rms*math.sqrt(50)*10/sqrt(beamarea_pix)
    error=np.full((data_cut.shape[0],data_cut.shape[1]),error_value)
    uncertainty=aperture_photometry(data_cut,apertures=aperture,error=error)['aperture_sum_err'][0]

    return flux, uncertainty

def mass_calc(flux,uncertainty):
    mass=1.05*10**4*(0.6/2)*flux*D**2
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

Image=imageDir+'NGC5258_12CO10_combine_contsub_mom0.fits'
wcs_cut=fits_import(Image)[0]
data_cut=fits_import(Image)[1]
shape=data_cut.shape
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
for i in range(int(steps-1)):
    rings[i]=aperture_ring(radius_arcsec[i+1],wcs_cut)
    rings_mask[i]=Apmask_convert(rings[i],data_cut)

fig=plt.figure()
ax=plt.subplot(projection=wcs_cut)
ax.imshow(data_cut,origin='lower')
# center_pix.plot(color='red')
# for i in range(size):
#     rings[i].plot(color='red')
# plt.show() 
rings[3].plot(color='red')
fig.savefig('../picture/turnover_h2.png')

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
uncertainty=flux_aperture_get(data_cut,center_pix,rms)[1]
uncertainties=np.append(uncertainties, uncertainty)
for i in range(len(rings)):
    uncertainty= flux_aperture_get(data_cut,rings[i],rms)[1]
    uncertainties=np.append(uncertainties, uncertainty)

# select the region with enough S/N ratio
flux=np.copy(flux[0:17])
uncertainties=np.copy(uncertainties[0:17])
radius_kpc=np.copy(radius_kpc[0:17])

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
areas=counts*0.09*500**2/incl
areas=areas[0:17]
SB=mass_ring/areas
fig=plt.figure()
plt.plot(radius_kpc,SB,linestyle=None,marker='o')

SB1=SB[2:]
radius_ksub=radius_kpc[2:]
popt, pcov = curve_fit(exp_func,radius_ksub,SB1,p0=(1, 0.1))
SB_mod=exp_func(radius_kpc,*popt)  
plt.plot(radius_kpc,SB_mod)
plt.xlabel('radius(kpc)')
plt.ylabel('surface brightness (M_sun/pc^2)')
plt.title('surface brightness of the molecular gas')

