'''
Apr. 24th
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

# Comparision with measurement from ds9.
'''
value={}
value['north_arm']={}
value['north_arm']['mean']=4.3363494
value['north_arm']['median']=3.9840026

value['south_arm']={}
value['south_arm']['mean']=3.2792174
value['south_arm']['median']=3.2592201

value['center']={}
value['center']['mean']=3.134277
value['center']['median']=3.1407537

value['center_ring']={}
value['center_ring']['mean']=3.0586879
value['center_ring']['median']=3.0209136

'''


############################################################
# measure the flux ratio of 12CO2-1/1-0

# basic setting
Dir='/1/home/heh15/workingspace/Arp240/ratio/'
pictureDir=Dir+'picture/'
scriptDir=Dir+'script/'
imageDir=Dir+'NGC5258/2110/'
measureDir=Dir+'NGC5258/test_measure/2110/'
image_12CO10='NGC5258_12CO10_combine_uvrange_smooth_regrid21_masked'
image_12CO21='NGC5258_12CO21_combine_uvtaper_smooth_masked'
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

origDir=os.getcwd()
os.chdir(measureDir)

##############################
# calculate the flux and uncertainty of 12CO1-0.

#copytree(imageDir+image_12CO10,image_12CO10)
#exportfits(imagename=image_12CO10,fitsimage=image_12CO+'.fits')
fitsimage=image_12CO10+'.fits'


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


# draw apertures in the image
fig=plt.figure()
ax=plt.subplot(projection=cut.wcs)
ax.imshow(data_cut,cmap='gray',origin='lower')
apertures['center']['pix'].plot(color='red')
apertures['ring around center']['pix'].plot(color='red')
apertures['northarm']['pix'].plot(color='red')
apertures['southarm']['pix'].plot(color='red')
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_major_formatter('hh:mm:ss')
plt.show()


# Calculate the flux and uncertainties
shape=data_cut.shape
rms_12CO10=0.0015
chan_width=10
beam_area=2.186*1.896*1.1331
beam_area_pix=beam_area/0.01
error_value=rms_12CO10*10*sqrt(50/beam_area_pix)
error_12CO10=np.full((shape[0],shape[1]),error_value)
data_cut=data_cut/beam_area_pix

for region in regions:
    stat10[region]['flux']=aperture_photometry(data_cut,apertures=apertures[region]['pix'],error=error_12CO10,mask=mask)['aperture_sum'][0]
    stat10[region]['uncertainty']=aperture_photometry(data_cut,apertures=apertures[region]['pix'],error=error_12CO10,mask=mask)['aperture_sum_err'][0]

##############################
# Calculate the flux and uncertainty of 12CO2-1
fitsimage=image_12CO21+'.fits'

# input the 12CO2-1 data
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


# draw apertures in the image
fig=plt.figure()
ax=plt.subplot(projection=cut.wcs)
ax.imshow(data_cut,cmap='gray',origin='lower')
apertures['center']['pix'].plot(color='red')
apertures['ring around center']['pix'].plot(color='red')
apertures['northarm']['pix'].plot(color='red')
apertures['southarm']['pix'].plot(color='red')
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_major_formatter('hh:mm:ss')
plt.show()
plt.savefig('NGC5258_12CO2-1.png')

os.remove(pictureDir+'NGC5258_12CO2-1.png')
shutil.move('NGC5258_12CO2-1.png',pictureDir)


# Calculate the flux and uncertainties
shape=data_cut.shape
rms_12CO21=0.004
chan_width=10
chans=50
beam_area=2.186*1.896*1.1331
beam_area_pix=beam_area/0.01
error_value=rms_12CO21*chan_width*sqrt(chans/beam_area_pix)
error_12CO21=np.full((shape[0],shape[1]),error_value)
data_cut=data_cut/beam_area_pix

for region in regions:
    stat21[region]['flux']=aperture_photometry(data_cut,apertures=apertures[region]['pix'],error=error_12CO21,mask=mask)['aperture_sum'][0]
    stat21[region]['uncertainty']=aperture_photometry(data_cut,apertures=apertures[region]['pix'],error=error_12CO21,mask=mask)['aperture_sum_err'][0]


# Calculate the peak, mean and median brightness of 12CO2-1

for region in regions: 
    region_mask_tmp=apertures[region]['pix'].to_mask(method='center')[0]
    region_mask_num=region_mask_tmp.to_image(shape=((shape[0],shape[1])))
    region_mask=np.ma.mask_or(region_mask_num==0,mask)
    region_masked=np.ma.masked_where(region_mask,data_cut)
    stat21[region]['peak']=np.ma.max(region_masked)*beam_area_pix
    stat21[region]['mean']=np.ma.mean(region_masked)*beam_area_pix
    stat21[region]['median']=np.ma.median(region_masked)*beam_area_pix
    stat21[region]['area']=region_masked.count()*0.01

##############################
# Calculate luminosity weighted ratio and uncertainty

for region in regions:
    stat[region]['ratio']=stat21[region]['flux']/stat10[region]['flux']*107.78**2/225.46**2
    stat[region]['uncertainty']=stat[region]['ratio']*\
                                math.sqrt((stat10[region]['uncertainty']/stat10[region]['flux'])**2+(stat21[region]['uncertainty']/stat21[region]['flux'])**2)
##############################
# record the value into the file.

# record the ratio
output_ratio=[]
header=ratio

for region in regions:
    tmp=[]
    tmp.append('{0: <20}'.format(region))
    for value in header:
        tmp1=round(stat[region][value],2)
        tmp1=str(tmp1)
        tmp1='{0: <12}'.format(tmp1)
        tmp.append(tmp1)
    output_ratio.append(tmp)

for i in range(len(header)):
    header[i]='{0: <12}'.format(header[i])

header.insert(0,'                    ')
output_ratio.insert(0,header)
 
filename='NGC5258_2110_ratio_quadrature.txt'
with open(filename,'w') as output:
    for i in range(len(output_ratio)):
        for j in range(len(output_ratio[i])):
            if j==(len(output_ratio[i])-1):
                output.write(output_ratio[i][j]+'\n')
            else:
                output.write(output_ratio[i][j])

os.remove(scriptDir+filename)
shutil.move(filename,scriptDir)

# record the brightness
output_12CO21=[]
header=values2
width=16
for region in regions:
    tmp=[]
    tmp.append('{0: <20}'.format(region))
    for value in values2:
        tmp1=round(stat21[region][value],2)
        tmp1=str(tmp1)
        tmp1='{0: <12}'.format(tmp1)
        tmp.append(tmp1)
    output_12CO21.append(tmp)

for i in range(len(header)):
    header[i]='{0: <12}'.format(header[i])

header.insert(0,'                    ')
output_12CO21.insert(0,header)


filename='NGC5258_12CO2-1_brightness.txt'
with open(filename,'w') as output:
    for i in range(len(output_12CO21)):
        for j in range(len(output_12CO21[i])):
            if j==(len(output_12CO21[i])-1):
                output.write(output_12CO21[i][j]+'\n')
            else:
                output.write(output_12CO21[i][j])

os.remove(scriptDir+filename)
shutil.move(filename,scriptDir)

os.chdir(scriptDir)
