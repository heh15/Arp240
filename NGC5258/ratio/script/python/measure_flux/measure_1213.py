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

# comparison with measurement from ds9
'''
value={}
value['north_arm']={}
value['north_arm']['mean']=15.805886
value['north_arm']['median']=14.15788

value['south_arm']={}
value['south_arm']['mean']=16.725381
value['south_arm']['median']=16.100594

value['center']={}
value['center']['mean']=9.72295
value['center']['median']=8.5737753

value['center_ring']={}
value['center_ring']['mean']=12.757023
value['center_ring']['median']=11.060942
'''

############################################################
# measure the flux ratio of 12CO/13CO

# basic setting
Dir='/1/home/heh15/workingspace/Arp240/ratio/'
imageDir=Dir+'NGC5258/1213/'
measureDir=Dir+'NGC5258/test_measure/1213/'
image_12CO10='NGC5258_12CO10_combine_smooth_masked'
image_13CO10='NGC5258_13CO10_12m_smooth_masked'
scriptDir=Dir+'script/'
pictureDir=Dir+'picture/'
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
data_cut.shape
rms_12CO10=0.0016
chan_width=10
beam_area=2.186*1.896*1.1331
beam_area_pix=beam_area/0.09
error_value=rms_12CO10*10*sqrt(50/beam_area_pix)
error_12CO10=np.full((180,140),error_value)
data_cut=data_cut/beam_area_pix

for region in regions:
    stat12[region]['flux']=aperture_photometry(data_cut,apertures=apertures[region]['pix'],error=error_12CO10,mask=mask)['aperture_sum'][0]
    stat12[region]['uncertainty']=aperture_photometry(data_cut,apertures=apertures[region]['pix'],error=error_12CO10,mask=mask)['aperture_sum_err'][0]

##############################
# Calculate the flux and uncertainty of 13CO1-0
fitsimage=image_13CO10+'.fits'

# input the 13CO1-0 data
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

# Calculate the flux and uncertainties
data_cut.shape
rms_13CO10=0.0006
chan_width=10
beam_area=2.186*1.896*1.1331
beam_area_pix=beam_area/0.09
error_value=rms_13CO10*20*sqrt(25/beam_area_pix)
error_13CO10=np.full((320,320),error_value)
data_cut=data_cut/beam_area_pix

for region in regions:
    stat13[region]['flux']=aperture_photometry(data_cut,apertures=apertures[region]['pix'],error=error_13CO10,mask=mask)['aperture_sum'][0]
    stat13[region]['uncertainty']=aperture_photometry(data_cut,apertures=apertures[region]['pix'],error=error_13CO10,mask=mask)['aperture_sum_err'][0]


##############################
# Calculate luminosity weighted ratio and uncertainty

for region in regions:
    stat[region]['ratio']=stat12[region]['flux']/stat13[region]['flux']
    stat[region]['uncertainty']=stat[region]['ratio']*\
                                math.sqrt((stat12[region]['uncertainty']/stat12[region]['flux'])**2+(stat13[region]['uncertainty']/stat13[region]['flux'])**2)

##############################
# record the value into the file

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
 
filename='NGC5258_1213_ratio_quadrature.txt'
with open(filename,'w') as output:
    for i in range(len(output_ratio)):
        for j in range(len(output_ratio[i])):
            if j==(len(output_ratio[i])-1):
                output.write(output_ratio[i][j]+'\n')
            else:
                output.write(output_ratio[i][j])

os.remove(scriptDir+filename)
shutil.move(filename,scriptDir)
