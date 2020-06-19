
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
from matplotlib import ticker
from regions import read_ds9

##plot setting
labelsize=20
nbins=5

############################################################
# directories and files. 

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/ratio/'
maskDir=Dir+'mask/'
regionDir=Dir+'region/'
imageDir=Dir+'image/ratio/1213/contsub_pbcor/'
# measureDir=Dir+'NGC5257/test_measure/1213/'
pictureDir=Dir+'picture/'
scriptDir=Dir+'script/'
image_12CO10='NGC5257_12CO10_combine_smooth_masked'
image_13CO10='NGC5257_13CO10_12m_smooth_masked'
image_ratio=imageDir+'NGC5257_1213_ratio_pbcor'

############################################################
# basic information
galaxy='NGC5257'
rationame='12CO/13CO 1-0'
ratiolabel='1213'
freq12=112
freq13=107

position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)

stat12={}
stat13={}
stat={}
apertures={}
regions=['center','first ring','second ring','third ring','fourth ring']
values=['flux','uncertainty']
ratio=['ratio','uncertainty']
type=['sky','pix']
stat12=dict.fromkeys(regions,{})
stat13=dict.fromkeys(regions,{})
apertures=apertures.fromkeys(regions,{})
stat=dict.fromkeys(regions,{})
for region in regions:
    stat12[region]=dict.fromkeys(values)
    stat13[region]=dict.fromkeys(values)
    apertures[region]=dict.fromkeys(type)
    stat[region]=dict.fromkeys(ratio)

############################################################
# basic settings

testra = 204.97609228
testdec = 0.84611111

linelabel = '$^{12}$CO 2-1'


############################################################
# function

def fits_import(fitsimage, item=0):
    hdr = fits.open(fitsimage)[item].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[item].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

############################################################
# Draw the region 

# position=SkyCoord(dec=50.415*u.arcmin,ra=204.9705*u.degree,frame='icrs')
# apertures['center']['sky']=SkyEllipticalAperture(position,a=2*u.arcsec,b=3*u.arcsec,theta=0*u.degree)
# apertures['center']['pix']=apertures['center']['sky'].to_pixel(wcs=cut.wcs)

# apertures[regions[1]]['sky']=SkyEllipticalAnnulus(position,a_in=2*u.arcsec,a_out=4*u.arcsec,b_out=6*u.arcsec,theta=0*u.degree)
# apertures[regions[1]]['pix']=apertures[regions[1]]['sky'].to_pixel(wcs=cut.wcs)

# apertures[regions[2]]['sky']=SkyEllipticalAnnulus(position,a_in=4*u.arcsec,a_out=6*u.arcsec,b_out=9*u.arcsec,theta=0*u.degree)
# apertures[regions[2]]['pix']=apertures[regions[2]]['sky'].to_pixel(wcs=cut.wcs)

# apertures[regions[3]]['sky']=SkyEllipticalAnnulus(position,a_in=6*u.arcsec,a_out=9*u.arcsec,b_out=13.5*u.arcsec,theta=0*u.degree)
# apertures[regions[3]]['pix']=apertures[regions[3]]['sky'].to_pixel(wcs=cut.wcs)

# apertures[regions[4]]['sky']=SkyEllipticalAnnulus(position,a_in=9*u.arcsec,a_out=15*u.arcsec,b_out=22.5*u.arcsec,theta=0*u.degree)
# apert

## spiral arm
maskimage=maskDir+'spiralarm_mask.fits'
wcs=fits_import(maskimage)[0]
spiral_masked=fits_import(maskimage)[1]
spiral_cut=Cutout2D(data=spiral_masked,position=position,size=size,wcs=wcs).data
wcs_cut=Cutout2D(data=spiral_masked,position=position,size=size,wcs=wcs).wcs

maskimage=maskDir+'anomaly_mask.fits'
anomaly_masked=fits_import(maskimage)[1]
anomaly_cut=Cutout2D(data=anomaly_masked,position=position,size=size,wcs=wcs).data

maskimage=maskDir+'disk_mask.fits'
restdisk_masked=fits_import(maskimage)[1]
restdisk_cut=Cutout2D(data=restdisk_masked,position=position,size=size,wcs=wcs).data

## south continuum source 
file=regionDir+'south.reg'
south_sky=read_ds9(file,errors='warn')[0]
south_pix=south_sky.to_pixel(wcs_cut)


# origDir=os.getcwd()
# os.chdir(measureDir)

fitsimage=image_ratio+'.fits'

# import data
hdr = fits.open(fitsimage)[0].header
wcs = WCS(hdr).celestial
data=fits.open(fitsimage)[0].data[0][0]
data=data*freq13**2/freq12**2

# cut out the region with data
position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)
cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
data_cut=cut.data

####################
# private region

# mask the data. 
data_masked=np.ma.masked_invalid(data_cut)
mask=data_masked.mask
co12_masked=data_masked

####################

fig=plt.figure()
ax=plt.subplot(projection=cut.wcs)
im=ax.imshow(data_cut,cmap='rainbow',origin='lower',vmax=30,vmin=0)
ax.tick_params(labelsize=8, direction='in')
ax.contour(spiral_cut, levels=[0.5])
ax.contour(anomaly_cut, levels=[0.5])
south_pix.plot(color='black')

x, y = (cut.wcs).wcs_world2pix(testra, testdec, 1)
plt.text(x, y, galaxy + ' ' + rationame, fontsize=15)

# apertures['center']['pix'].plot(color='red')
# apertures[regions[1]]['pix'].plot(color='red')
# apertures[regions[2]]['pix'].plot(color='red')
# apertures[regions[3]]['pix'].plot(color='red')
# apertures[regions[4]]['pix'].plot(color='red')
lon = ax.coords[0]
lat = ax.coords[1]
# lon.set_ticklabel_visible(False)
# lat.set_ticklabel_visible(False)
cax=fig.add_axes()
cbar=fig.colorbar(im,cax=cax)
# cbar.set_label(rationame, fontsize=20)
cbar.ax.tick_params(labelsize=labelsize)
tick_locator = ticker.MaxNLocator(nbins=nbins)
cbar.locator = tick_locator
cbar.update_ticks() 
fig.savefig(pictureDir+galaxy+'_'+ratiolabel+'_ratio_paper_aperture.png')
# plt.show()
