
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

# ##plot setting
labelsize=20
nbins=5

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/ratio/'
imageDir=Dir+'image/ratio/1213/contsub_pbcor/'
# measureDir=Dir+'NGC5258/test_measure/1213/'
image_12CO10='NGC5258_12CO10_combine_smooth_masked'
image_13CO10='NGC5258_13CO10_12m_smooth_masked'
image_ratio=imageDir+'NGC5258_1213_ratio_pbcor'
scriptDir=Dir+'script/'
pictureDir=Dir+'picture/'
regionDir=Dir+'region/'

stat12={}
stat13={}
stat={}
apertures={}
regions=['northarm','southarm','center','ring around center']
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

# origDir=os.getcwd()
# os.chdir(measureDir)

fitsimage=image_ratio+'.fits'

# import data
hdr = fits.open(fitsimage)[0].header
wcs = WCS(hdr).celestial
data=fits.open(fitsimage)[0].data[0][0]

# cut out the region with data
position=SkyCoord(dec=0.8319*u.degree,ra=204.9908*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)
cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
data_cut=cut.data

# mask the data. 

data_masked=np.ma.masked_invalid(data_cut)
mask=data_masked.mask

# draw the shape of the region.

position=SkyCoord(dec=0.8309*u.degree,ra=204.9906*u.degree,frame='icrs')
apertures['center']['sky']=SkyCircularAperture(position,r=3*u.arcsec)
apertures['center']['pix']=apertures['center']['sky'].to_pixel(wcs=cut.wcs)

apertures['ring around center']['sky']=SkyCircularAnnulus(position,r_in=3*u.arcsec,r_out=7*u.arcsec)
apertures['ring around center']['pix']=apertures['ring around center']['sky'].to_pixel(wcs=cut.wcs)

position=SkyCoord(dec=0.8349*u.degree,ra=204.9930*u.degree,frame='icrs')
apertures['northarm']['sky']=SkyEllipticalAperture(position,a=15*u.arcsec,b=7*u.arcsec,theta=340*u.degree)
apertures['northarm']['pix']=apertures['northarm']['sky'].to_pixel(wcs=cut.wcs)

position=SkyCoord(dec=0.8275*u.degree,ra=204.9884*u.degree,frame='icrs')
apertures['southarm']['sky']=SkyEllipticalAperture(position,a=10*u.arcsec,b=5*u.arcsec,theta=340*u.degree)
apertures['southarm']['pix']=apertures['southarm']['sky'].to_pixel(wcs=cut.wcs)

## import the regions. 

from regions import read_ds9

file=regionDir+'south_SFR.reg'
southsfr_sky=read_ds9(file,errors='warn')[0]
southsfr_pix=southsfr_sky.to_pixel(cut.wcs)

# draw apertures in the image
fig=plt.figure()
ax=plt.subplot(projection=cut.wcs)
im=ax.imshow(data_cut,cmap='rainbow',origin='lower',vmax=30,vmin=0)
apertures['center']['pix'].plot(color='black')
apertures['ring around center']['pix'].plot(color='black')
apertures['northarm']['pix'].plot(color='black')
apertures['southarm']['pix'].plot(color='black')
southsfr_pix.plot(color='black')
lon = ax.coords[0]
lat = ax.coords[1]
# lon.set_major_formatter('hh:mm:ss')
# lon.set_ticklabel_visible(False)
# lat.set_ticklabel_visible(False)
cax=fig.add_axes()
cbar=fig.colorbar(im,cax=cax)
cbar.set_label('12CO/13CO 1-0 ratio', fontsize=20)
cbar.ax.tick_params(labelsize=labelsize)
tick_locator = ticker.MaxNLocator(nbins=nbins)
cbar.locator = tick_locator
cbar.update_ticks() 
plt.show()
fig.savefig(pictureDir+'NGC5258_1213_ratio_paper_aperture.png')


# ############################################################
# # NGC 5257 1213

# Dir='/1/home/heh15/workingspace/Arp240/ratio/'
# imageDir=Dir+'NGC5257/1213/combine/'
# measureDir=Dir+'NGC5257/test_measure/1213/'
# scriptDir=Dir+'script/'
# image_12CO10='NGC5257_12CO10_combine_smooth_masked'
# image_13CO10='NGC5257_13CO10_12m_smooth_masked'
# image_ratio='NGC5257_1213_ratio'
# stat12={}
# stat13={}
# stat={}
# apertures={}
# regions=['center','first ring','second ring','third ring','fourth ring']
# values=['flux','uncertainty']
# ratio=['ratio','uncertainty']
# type=['sky','pix']
# stat12=dict.fromkeys(regions,{})
# stat13=dict.fromkeys(regions,{})
# apertures=apertures.fromkeys(regions,{})
# stat=dict.fromkeys(regions,{})
# for region in regions:
#     stat12[region]=dict.fromkeys(values)
#     stat13[region]=dict.fromkeys(values)
#     apertures[region]=dict.fromkeys(type)
#     stat[region]=dict.fromkeys(ratio)

# origDir=os.getcwd()
# os.chdir(measureDir)

# fitsimage=image_ratio+'.fits'

# # import data
# hdr = fits.open(fitsimage)[0].header
# wcs = WCS(hdr).celestial
# data=fits.open(fitsimage)[0].data[0][0]

# # cut out the region with data
# position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')
# size=u.Quantity((54,42),u.arcsec)
# cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
# data_cut=cut.data

# # mask the data. 

# data_masked=np.ma.masked_invalid(data_cut)
# mask=data_masked.mask
# co12_masked=data_masked

# # draw the apertures in the image
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
# apertures[regions[4]]['pix']=apertures[regions[4]]['sky'].to_pixel(wcs=cut.wcs)

# fig=plt.figure()
# ax=plt.subplot(projection=cut.wcs)
# im=ax.imshow(data_cut,cmap='rainbow',origin='lower',vmax=60,vmin=0)
# apertures['center']['pix'].plot(color='red')
# apertures[regions[1]]['pix'].plot(color='red')
# apertures[regions[2]]['pix'].plot(color='red')
# apertures[regions[3]]['pix'].plot(color='red')
# apertures[regions[4]]['pix'].plot(color='red')
# lon = ax.coords[0]
# lat = ax.coords[1]
# lon.set_ticklabel_visible(False)
# lat.set_ticklabel_visible(False)
# cax=fig.add_axes()
# cbar=fig.colorbar(im,cax=cax)
# cbar.ax.tick_params(labelsize=labelsize)
# tick_locator = ticker.MaxNLocator(nbins=nbins)
# cbar.locator = tick_locator
# cbar.update_ticks() 
# fig.savefig(pictureDir+'NGC5257_1213_ratio_poster')
# plt.show()


############################################################
# NGC 5258 2110
Dir='/1/home/heh15/workingspace/Arp240/NGC5258/ratio/'
pictureDir=Dir+'picture/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/ratio/2110/uvtaper_5rms/'
image_12CO10='NGC5258_12CO10_combine_uvrange_smooth_regrid21_masked'
image_12CO21='NGC5258_12CO21_combine_uvtaper_smooth_masked'
image_ratio=imageDir+'NGC5258_2110_ratio_uvtaper_pbcor'
stat10={}
stat21={}
stat={}
apertures={}
regions=['northarm','southarm','center','ring around center']
values=['flux','uncertainty']
values2=['flux','uncertainty','peak','mean','median','area']
ratio=['ratio','uncertainty']
type=['sky','pix']
stat10=dict.fromkeys(regions,{})
stat21=dict.fromkeys(regions,{})
apertures=apertures.fromkeys(regions,{})
stat=dict.fromkeys(regions,{})
for region in regions:
    stat10[region]=dict.fromkeys(values)
    stat21[region]=dict.fromkeys(values2)
    apertures[region]=dict.fromkeys(type)
    stat[region]=dict.fromkeys(ratio)

# origDir=os.getcwd()
# os.chdir(measureDir)

fitsimage=image_ratio+'.fits'

hdr = fits.open(fitsimage)[0].header
wcs = WCS(hdr).celestial
data=np.squeeze(fits.open(fitsimage)[0].data)
data=data*112.73**2/225.46**2

# cut out the region with data
position=SkyCoord(dec=0.8319*u.degree,ra=204.9908*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)
cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
data_cut=cut.data

# mask the data. 

data_masked=np.ma.masked_invalid(data_cut)
mask=data_masked.mask

position=SkyCoord(dec=0.8309*u.degree,ra=204.9906*u.degree,frame='icrs')
apertures['center']['sky']=SkyCircularAperture(position,r=3*u.arcsec)
apertures['center']['pix']=apertures['center']['sky'].to_pixel(wcs=cut.wcs)

apertures['ring around center']['sky']=SkyCircularAnnulus(position,r_in=3*u.arcsec,r_out=7*u.arcsec)
apertures['ring around center']['pix']=apertures['ring around center']['sky'].to_pixel(wcs=cut.wcs)

position=SkyCoord(dec=0.8349*u.degree,ra=204.9930*u.degree,frame='icrs')
apertures['northarm']['sky']=SkyEllipticalAperture(position,a=15*u.arcsec,b=7*u.arcsec,theta=340*u.degree)
apertures['northarm']['pix']=apertures['northarm']['sky'].to_pixel(wcs=cut.wcs)

position=SkyCoord(dec=0.8275*u.degree,ra=204.9884*u.degree,frame='icrs')
apertures['southarm']['sky']=SkyEllipticalAperture(position,a=10*u.arcsec,b=5*u.arcsec,theta=340*u.degree)
apertures['southarm']['pix']=apertures['southarm']['sky'].to_pixel(wcs=cut.wcs)

# draw apertures in the image
fig=plt.figure()
ax=plt.subplot(projection=cut.wcs)
im=ax.imshow(data_cut,cmap='rainbow',origin='lower',vmin=0.5,vmax=1.1)
apertures['center']['pix'].plot(color='black')
apertures['ring around center']['pix'].plot(color='black')
apertures['northarm']['pix'].plot(color='black')
apertures['southarm']['pix'].plot(color='black')
southsfr_pix.plot(color='black')
lon = ax.coords[0]
lat = ax.coords[1]
# lon.set_ticklabel_visible(False)
# lat.set_ticklabel_visible(False)
cax=fig.add_axes()
cbar=fig.colorbar(im,cax=cax)
cbar.set_label('12CO 2-1/1-0 ratio', fontsize=20)
cbar.ax.tick_params(labelsize=labelsize)
tick_locator = ticker.MaxNLocator(nbins=nbins)
cbar.locator = tick_locator
cbar.update_ticks() 
fig.savefig(pictureDir+'NGC5258_2110_ratio_paper_aperture.png')
plt.show()
plt.show()


# ############################################################
# # NGC 5257 2110

# Dir='/1/home/heh15/workingspace/Arp240/ratio/'
# pictureDir=Dir+'picture/'
# scriptDir=Dir+'script/'
# imageDir=Dir+'NGC5257/2110/uvtaper/'
# measureDir=Dir+'NGC5257/test_measure/2110/'
# image_12CO10='NGC5257_12CO10_combine_uvrange_smooth_regrid21_masked'
# image_12CO21='NGC5257_12CO21_combine_uvtaper_smooth_masked'
# image_ratio='NGC5257_2110_ratio_uvtaper.image'
# stat10={}
# stat21={}
# stat={}
# apertures={}
# regions=['center','first ring','second ring','third ring','forth ring']
# values=['flux','uncertainty']
# values2=['flux','uncertainty','peak','mean','median','area']
# ratio=['ratio','uncertainty']
# type=['sky','pix']
# stat10=dict.fromkeys(regions,{})
# stat21=dict.fromkeys(regions,{})
# apertures=apertures.fromkeys(regions,{})
# stat=dict.fromkeys(regions,{})
# for region in regions:
#     stat10[region]=dict.fromkeys(values)
#     stat21[region]=dict.fromkeys(values2)
#     apertures[region]=dict.fromkeys(type)
#     stat[region]=dict.fromkeys(ratio)

# origDir=os.getcwd()
# os.chdir(measureDir)

# fitsimage=image_ratio+'.fits'

# hdr = fits.open(fitsimage)[0].header
# wcs = WCS(hdr).celestial
# data=fits.open(fitsimage)[0].data[0][0]
# data=data*107.78**2/225.46**2

# position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')
# size=u.Quantity((54,42),u.arcsec)
# cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
# data_cut=cut.data

# # mask the data. 
# data_masked=np.ma.masked_invalid(data_cut)
# mask=data_masked.mask
# co10_masked=data_masked

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
# apertures[regions[4]]['pix']=apertures[regions[4]]['sky'].to_pixel(wcs=cut.wcs)


# fig=plt.figure()
# ax=plt.subplot(projection=cut.wcs)
# im=ax.imshow(data_cut,cmap='rainbow',origin='lower',vmin=0,vmax=3)
# apertures['center']['pix'].plot(color='red')
# apertures[regions[1]]['pix'].plot(color='red')
# apertures[regions[2]]['pix'].plot(color='red')
# apertures[regions[3]]['pix'].plot(color='red')
# apertures[regions[4]]['pix'].plot(color='red')
# lon = ax.coords[0]
# lat = ax.coords[1]
# lon.set_ticklabel_visible(False)
# lat.set_ticklabel_visible(False)
# cax=fig.add_axes()
# cbar=fig.colorbar(im,cax=cax)
# cbar.ax.tick_params(labelsize=labelsize)
# tick_locator = ticker.MaxNLocator(nbins=nbins)
# cbar.locator = tick_locator
# cbar.update_ticks() 
# # fig.savefig(pictureDir+'NGC5257_2110_ratio_poster')
# plt.show()
