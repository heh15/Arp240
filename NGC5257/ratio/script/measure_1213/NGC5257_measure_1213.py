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


############################################################
# measure the flux ratio of 12CO/13CO

# basic setting
Dir='/1/home/heh15/workingspace/Arp240/ratio/'
imageDir=Dir+'NGC5257/1213/combine/'
measureDir=Dir+'NGC5257/test_measure/1213/'
scriptDir=Dir+'script/'
image_12CO10='NGC5257_12CO10_combine_smooth_masked'
image_13CO10='NGC5257_13CO10_12m_smooth_masked'
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

origDir=os.getcwd()
os.chdir(measureDir)

##############################
# calculate the flux and uncertainty of 12CO1-0.

if os.path.isdir(image_12CO10+'.image.mom0')==False:
    copytree(imageDir+image_12CO10+'.image.mom0',image_12CO10+'.image.mom0')

# in casa 
# exportfits(imagename=image_12CO10+'.image.mom0',fitsimage=image_12CO10+'.fits')
fitsimage=image_12CO10+'.fits'

# import data
hdr = fits.open(fitsimage)[0].header
wcs = WCS(hdr).celestial
data=fits.open(fitsimage)[0].data[0][0]

# cut out the region with data
position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)
cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
data_cut=cut.data

# mask the data. 

data_masked=np.ma.masked_invalid(data_cut)
mask=data_masked.mask
co12_masked=data_masked

# draw the apertures in the image
position=SkyCoord(dec=50.415*u.arcmin,ra=204.9705*u.degree,frame='icrs')
apertures['center']['sky']=SkyEllipticalAperture(position,a=2*u.arcsec,b=3*u.arcsec,theta=0*u.degree)
apertures['center']['pix']=apertures['center']['sky'].to_pixel(wcs=cut.wcs)

apertures[regions[1]]['sky']=SkyEllipticalAnnulus(position,a_in=2*u.arcsec,a_out=4*u.arcsec,b_out=6*u.arcsec,theta=0*u.degree)
apertures[regions[1]]['pix']=apertures[regions[1]]['sky'].to_pixel(wcs=cut.wcs)

apertures[regions[2]]['sky']=SkyEllipticalAnnulus(position,a_in=4*u.arcsec,a_out=6*u.arcsec,b_out=9*u.arcsec,theta=0*u.degree)
apertures[regions[2]]['pix']=apertures[regions[2]]['sky'].to_pixel(wcs=cut.wcs)

apertures[regions[3]]['sky']=SkyEllipticalAnnulus(position,a_in=6*u.arcsec,a_out=9*u.arcsec,b_out=13.5*u.arcsec,theta=0*u.degree)
apertures[regions[3]]['pix']=apertures[regions[3]]['sky'].to_pixel(wcs=cut.wcs)

apertures[regions[4]]['sky']=SkyEllipticalAnnulus(position,a_in=9*u.arcsec,a_out=15*u.arcsec,b_out=22.5*u.arcsec,theta=0*u.degree)
apertures[regions[4]]['pix']=apertures[regions[4]]['sky'].to_pixel(wcs=cut.wcs)

fig=plt.figure()
ax=plt.subplot(projection=cut.wcs)
ax.imshow(data_cut,cmap='gray',origin='lower')
apertures['center']['pix'].plot(color='red')
apertures[regions[1]]['pix'].plot(color='red')
apertures[regions[2]]['pix'].plot(color='red')
apertures[regions[3]]['pix'].plot(color='red')
apertures[regions[4]]['pix'].plot(color='red')
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
error_co12=error_value

for region in regions:
    stat12[region]['flux']=aperture_photometry(data_cut,apertures=apertures[region]['pix'],error=error_12CO10,mask=mask)['aperture_sum'][0]
    stat12[region]['uncertainty']=aperture_photometry(data_cut,apertures=apertures[region]['pix'],error=error_12CO10,mask=mask)['aperture_sum_err'][0]

##############################
# Calculate the flux and uncertainty of 13CO1-0
fitsimage=image_13CO10+'.fits'

if os.path.isdir(image_13CO10+'.image.mom0')==False:
    copytree(imageDir+image_13CO10+'.image.mom0',image_13CO10+'.image.mom0')
# in casa 
# exportfits(imagename=image_13CO10+'.image.mom0',fitsimage=image_13CO10+'.fits')
fitsimage=image_13CO10+'.fits'

# input the 13CO1-0 data
hdr = fits.open(fitsimage)[0].header
wcs = WCS(hdr).celestial
data=fits.open(fitsimage)[0].data[0][0]

# cut out the region with data
position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)
cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
data_cut=cut.data
data_cut=data_cut

# mask the data. 

data_masked=np.ma.masked_invalid(data_cut)
mask=data_masked.mask
co13_masked=data_masked

# Calculate the flux and uncertainties
data_cut.shape
rms_13CO10=0.0006
chan_width=10
beam_area=2.186*1.896*1.1331
beam_area_pix=beam_area/0.09
error_value=rms_13CO10*20*sqrt(25/beam_area_pix)
error_13CO10=np.full((180,140),error_value)
data_cut=data_cut/beam_area_pix
error_co13=error_value

for region in regions:
    stat13[region]['flux']=aperture_photometry(data_cut,apertures=apertures[region]['pix'],error=error_13CO10,mask=mask)['aperture_sum'][0]
    stat13[region]['uncertainty']=aperture_photometry(data_cut,apertures=apertures[region]['pix'],error=error_13CO10,mask=mask)['aperture_sum_err'][0]


##############################
# Calculate luminosity weighted ratio and uncertainty

for region in regions:
    stat[region]['ratio']=stat12[region]['flux']/stat13[region]['flux']
    stat[region]['uncertainty']=stat[region]['ratio']*\
                                math.sqrt((stat12[region]['uncertainty']/stat12[region]['flux'])**2+(stat13[region]['uncertainty']/stat13[region]['flux'])**2+2*0.05**2)

##############################
# measure the total value and record it into the file

regions.append('total')
region='total'
stat13[region]={}
stat13[region]['flux']=np.ma.sum(co13_masked)/beam_area_pix
stat13[region]['uncertainty']=error_co13*math.sqrt(np.ma.count(co13_masked))

stat12[region]={}
error_12CO10_masked=np.ma.masked_where(error_12CO10,mask)
stat12[region]['flux']=np.ma.sum(co12_masked)/beam_area_pix
stat12[region]['uncertainty']=error_co12*math.sqrt(np.ma.count(co12_masked))

stat[region]={}
stat[region]['ratio']=stat12[region]['flux']/stat13[region]['flux']
stat[region]['uncertainty']=stat[region]['ratio']*\
                                math.sqrt((stat12[region]['uncertainty']/stat12[region]['flux'])**2+(stat13[region]['uncrtainty']/stat13[region]['flux'])**2)

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
 
filename='NGC5257_1213_ratio_quadrature.txt'
with open(filename,'w') as output:
    for i in range(len(output_ratio)):
        for j in range(len(output_ratio[i])):
            if j==(len(output_ratio[i])-1):
                output.write(output_ratio[i][j]+'\n')
            else:
                output.write(output_ratio[i][j])

os.remove(scriptDir+filename)
shutil.move(filename,scriptDir)



