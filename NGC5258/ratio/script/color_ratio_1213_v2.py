
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
from reproject import reproject_interp


##plot setting
labelsize=20
nbins=5

############################################################
# directories and files. 

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/ratio/'
maskDir=Dir+'mask/'
regionDir=Dir+'region/'
imageDir=Dir+'image/ratio/1213/contsub_pbcor/'
# measureDir=Dir+'NGC5257/test_measure/1213/'
pictureDir=Dir+'picture/'
scriptDir=Dir+'script/'
image_12CO10='NGC5258_12CO10_combine_smooth_masked'
image_13CO10='NGC5258_13CO10_12m_smooth_masked'
image_ratio=imageDir+'NGC5258_1213_ratio_pbcor'

############################################################
# basic information
galaxy='NGC5258'
rationame='12CO/13CO 1-0'
ratiolabel='1213'
freq12=112
freq13=107

center=SkyCoord(dec=0.8319*u.degree,ra=204.9908*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)

stat12={}
stat13={}
stat={}
apertures={}
regions=['center','ring around center','northarm','southarm']
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

testra = 204.99631309
testdec = 0.83772222

linelabel = '$^{12}$CO 2-1'

poster = 'off'

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

position=SkyCoord(dec=0.8309*u.degree,ra=204.9906*u.degree,frame='icrs')
apertures['center']['sky']=SkyCircularAperture(position,r=3*u.arcsec)
apertures['ring around center']['sky']=SkyCircularAnnulus(position,r_in=3*u.arcsec,r_out=7*u.arcsec)
position=SkyCoord(dec=0.8340*u.degree,ra=204.9935*u.degree,frame='icrs')
apertures['northarm']['sky']=SkyEllipticalAperture(position,a=13*u.arcsec,b=4*u.arcsec,theta=185*u.degree)
position=SkyCoord(dec=0.8283*u.degree,ra=204.9882*u.degree,frame='icrs')
apertures['southarm']['sky']=SkyEllipticalAperture(position,a=9*u.arcsec,b=4*u.arcsec,theta=360*u.degree)

from regions import read_ds9
file=regionDir+'south_SFR.reg'
southsfr_sky=read_ds9(file,errors='warn')[0]

fitsimage=image_ratio+'.fits'


###  import data
hdr = fits.open(fitsimage)[0].header
wcs = WCS(hdr).celestial
data=fits.open(fitsimage)[0].data[0][0]
data=data*freq13**2/freq12**2

# cut out the region with data
position=center
size=u.Quantity((54,42),u.arcsec)
cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
data_cut=cut.data

for key in apertures.keys():
    apertures[key]['pix']=apertures[key]['sky'].to_pixel(wcs=cut.wcs)
southsfr_pix=southsfr_sky.to_pixel(cut.wcs)

# ## spiral arm
# maskimage=maskDir+'spiralarm_mask.fits'
# wcs=fits_import(maskimage)[0]
# spiral_masked=fits_import(maskimage)[1]
# spiral_cut=Cutout2D(data=spiral_masked,position=position,size=size,wcs=wcs).data
# wcs_cut=Cutout2D(data=spiral_masked,position=position,size=size,wcs=wcs).wcs

# maskimage=maskDir+'anomaly_mask.fits'
# anomaly_masked=fits_import(maskimage)[1]
# anomaly_cut=Cutout2D(data=anomaly_masked,position=position,size=size,wcs=wcs).data

# maskimage=maskDir+'disk_mask.fits'
# restdisk_masked=fits_import(maskimage)[1]
# restdisk_cut=Cutout2D(data=restdisk_masked,position=position,size=size,wcs=wcs).data

# ## south continuum source 
# file=regionDir+'south.reg'
# south_sky=read_ds9(file,errors='warn')[0]
# south_pix=south_sky.to_pixel(wcs_cut)


# origDir=os.getcwd()
# os.chdir(measureDir)


####################
# private region

# mask the data. 
data_masked=np.ma.masked_invalid(data_cut)
mask=data_masked.mask
co12_masked=data_masked

####################

## import the 12CO2-1 moment 0 map.
co21_dir = Dir+'image/12CO21/'
fitsimage = co21_dir+'NGC5258_12CO21_combine_pbcor_mom0.fits' 
wcs_co21, data_masked = fits_import(fitsimage)

contour, footprint = reproject_interp((data_masked, wcs_co21), cut.wcs, shape_out=np.shape(data_cut)) 

color = 'maroon'
fig=plt.figure()
ax=plt.subplot(projection=cut.wcs)
im=ax.imshow(data_cut,cmap='rainbow',origin='lower',vmax=30,vmin=0)
ax.tick_params(labelsize=8, direction='in')
if poster == 'off':
	apertures['center']['pix'].plot(color=color, linewidth=2.0)
	apertures['ring around center']['pix'].plot(color=color, linewidth=2.0)
	apertures['northarm']['pix'].plot(color=color, linewidth=2.0)
	apertures['southarm']['pix'].plot(color=color, linewidth=2.0)
southsfr_pix.plot(color=color, linewidth=2.0)

ax.contour(contour, levels=[1.1, 2.2], colors='black', linewidths=0.8)

x, y = (cut.wcs).wcs_world2pix(testra, testdec, 1)
plt.text(x,y, galaxy + ' ' + rationame, fontsize=15)

plt.xlabel('J2000 Right Ascension')
plt.ylabel('J2000 Declination')

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
if poster == 'off':
	fig.savefig(pictureDir+galaxy+'_'+ratiolabel+'_ratio_paper_aperture.png', bbox_inches='tight')
else: 
	fig.savefig(pictureDir+galaxy+'_'+ratiolabel+'_ratio_poster_aperture.png', bbox_inches='tight')


