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

############################################################
# function

def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
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
Dir='/home/heh15/workingspace/Arp240/NGC5258/SM/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
pictureDir=Dir+'picture/'
logDir=Dir+'log/'

origDir=os.getcwd()
os.chdir(imageDir)

##############################
# measure the flux of spitzer 3.6um image

Image=imageDir+'spitzer_36um_regrid.fits'
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
for i in range(int(steps-1)):
    rings[i]=aperture_ring(radius_arcsec[i+1],wcs_cut)
    rings_mask[i]=Apmask_convert(rings[i],data_cut)


fig=plt.figure()
ax=plt.subplot(projection=wcs_cut)
ax.imshow(data_cut,origin='lower')
center_pix.plot(color='red')
# for i in range(11):
#     rings[i].plot(color='red')
rings[9].plot(color='red')
# rings[11].plot(color='green')
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_major_formatter('hh:mm:ss')
plt.show()
fig.savefig(pictureDir+'ellipse.png')


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

flux=flux*10**6*pixel_sr
area=counts*pixel_area
SB=np.divide(flux,area)

flux_36=flux

# measure the uncertainty of the surface brightness. 
stds_36=np.empty(shape=(0,0))
std_36=np.ma.std(center_mask)
stds_36=np.append(stds_36, std_36)
for i in range(len(rings)):
    std_36=np.ma.std(rings_mask[i])
    stds_36=np.append(stds_36,std_36)

stds_36=stds_36*10**6*pixel_sr

# measure the mean of the different region. 
means_36=np.empty(shape=(0,0))
mean_36=np.ma.mean(center_mask)
means_36=np.append(means_36, mean_36)
for i in range(len(rings)):
    mean_36=np.ma.mean(rings_mask[i])
    means_36=np.append(means_36,mean_36)

means_36=means_36*10**6*pixel_sr

# # measure the uncertainty within each ring.
# Image='SPITZER_I1_39933184_0000_2_A52342816_munc.fits'
# wcs_cut=fits_import(Image)[0]
# error_cut=fits_import(Image)[1]
# shape=data_cut.shape
# uncertainty=np.empty(shape=(0,0))

# uncertainty_center= aperture_photometry(data_cut,apertures=center_pix,error=error_cut)['aperture_sum_err'][0]
# uncertainty=np.append(uncertainty,uncertainty_center)

# for i in range(size):
#     uncertainty_ring=aperture_photometry(data_cut,apertures=rings[i],error=error_cut)['aperture_sum_err'][0]
#     uncertainty=np.append(uncertainty,uncertainty_ring)

# uncertainty=uncertainty*10**6*pixel_sr
# uncertainty_36=uncertainty

########################################
# measure the flux of spitzer 4.5um image. 

Image=imageDir+'spitzer_45um_regrid.fits'
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

# fig=plt.figure()
# ax=plt.subplot(projection=wcs_cut)
# ax.imshow(data_cut,origin='lower')
# center_pix.plot(color='red')
# for i in range(size):
#     rings[i].plot(color='red')
# plt.show()


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

flux=flux*10**6*pixel_sr
area=counts*pixel_area
SB=np.divide(flux,area)

flux_45=flux

# measure the uncertainty 
stds_45=np.empty(shape=(0,0))
std_45=np.ma.std(center_mask)
stds_45=np.append(stds_45, std_45)
for i in range(len(rings)):
    std_45=np.ma.std(rings_mask[i])
    stds_45=np.append(stds_45,std_45)

stds_45=stds_45*10**6*pixel_sr

# measure the mean 
means_45=np.empty(shape=(0,0))
mean_45=np.ma.mean(center_mask)
means_45=np.append(means_45, mean_45)
for i in range(len(rings)):
    mean_45=np.ma.mean(rings_mask[i])
    means_45=np.append(means_45,mean_45)

means_45=means_45*10**6*pixel_sr

# measure the uncertainty within each ring.
# Image='SPITZER_I2_39933184_0000_2_A52349581_munc.fits'
# wcs_cut=fits_import(Image)[0]
# error_cut=fits_import(Image)[1]
# shape=data_cut.shape
# uncertainty=np.empty(shape=(0,0))

# uncertainty_center= aperture_photometry(data_cut,apertures=center_pix,error=error_cut)['aperture_sum_err'][0]
# uncertainty=np.append(uncertainty,uncertainty_center)

# for i in range(size):
#     uncertainty_ring=aperture_photometry(data_cut,apertures=rings[i],error=error_cut)['aperture_sum_err'][0]
#     uncertainty=np.append(uncertainty,uncertainty_ring)

# uncertainty=uncertainty*10**6*pixel_sr
# uncertainty_45=uncertainty

##############################
# calculate the mass

## calculate the surface brightness

mass_ring=mass_calc(flux_36,flux_45)
error_ring=0.3*np.sqrt(10**6*mass_ring)
mass_ring[0]=mass_ring[0]/0.6

mass_output=np.copy(mass_ring[0:17])
error_output=np.copy(error_ring[0:17])
radius_output=np.copy(radius_kpc[0:17])

output=np.transpose(np.vstack([radius_output,mass_output,error_output]))
filename=logDir+'mstar_ring.txt'
fmt= '%10.3e %10.3e %10.3e \n'
np.savetxt(filename,output,fmt=fmt,delimiter='',newline='\n')

mass_radius_tmp=mass_sum(mass_ring)
# error_radius_tmp=mass_sum(error_ring)
error_radius_tmp=0.3*np.sqrt(10**6*mass_radius_tmp)

mass_radius=np.log10(mass_radius_tmp)
# error_radius=0.434*np.divide(error_radius_tmp,mass_radius_tmp)
error_radius=0.434*np.divide(error_radius_tmp,mass_radius_tmp)

# calculate the surface brightness of the mass. 
areas=np.empty(shape=(0,0))
area=math.pi*radius_kpc[0]**2*1000**2
areas=np.append(areas,area)

for i in range(mass_ring.shape[0]-1):
    area=math.pi*(radius_kpc[i+1]**2-radius_kpc[i]**2)*1000**2
    areas=np.append(areas,area)

mass_SB=np.divide(mass_ring,areas)
uncertainty_SB=error_ring/(areas)

fig=plt.figure()
plt.errorbar(radius_kpc,mass_SB,uncertainty_SB)

# popt, pcov = curve_fit(lambda t,a,b: a*numpy.exp(-b*t),radius_kpc,mass_SB,p0=(1, 0.1))
radius_ksub=radius_kpc[2:]
mass_Ssub=mass_SB[2:]
popt, pcov = curve_fit(exp_func,radius_ksub,mass_Ssub,p0=(1, 0.1))
mass_mod=exp_func(radius_kpc,*popt)    
plt.plot(radius_kpc, mass_mod)
plt.xlabel('radius (kpc)')
plt.ylabel('surface brightness (solar mass/pc^2)')
plt.title('Surface Brightness of Rings')
plt.savefig(pictureDir+'SB_fit')

# output the result
output=np.transpose(np.vstack((radius_kpc,mass_radius)))
fmt= '%10.3e %10.3e \n'
filename=logDir+'stellarmass.txt'
np.savetxt(filename,output,fmt=fmt,delimiter='',newline='\n')

SB_star=np.transpose(np.vstack((mass_SB,uncertainty_SB)))
fmt= '%10.3e %10.3e \n'
filename=logDir+'SB_star.txt'
np.savetxt(filename,SB_star,fmt=fmt,delimiter='',newline='\n')

## measure gas fraction

#  gas mass fraction.
m_star=np.transpose(np.loadtxt(logDir+'mstar_ring.txt'))
# m_star[0]=m_star[0]/0.6
m_gas=np.transpose(np.loadtxt(logDir+'H2mass_ring.txt'))
fraction=np.divide(m_gas[1],m_star[1]+m_gas[1])
error_gas=np.divide(m_gas[2],m_star[1]+m_gas[1])

def fraction_err(m_star,m_gas):
    error_no=m_gas[2]/m_gas[1]
    error_den=np.sqrt((m_gas[2]**2+m_star[2]**2)/(m_star[1]+m_gas[1])**2)
    error_frac=np.sqrt(error_no**2+error_den**2)
    error=error_frac*m_gas[1]/(m_gas[1]+m_star[1])
    return error

error=fraction_err(m_star,m_gas)    

fig=plt.figure()
plt.errorbar(radius_output,fraction,error)
plt.xlabel('radius (kpc)')
plt.ylabel('molecular gas fraction')
plt.savefig(pictureDir+'fraction.png')

os.chdir(origDir)
