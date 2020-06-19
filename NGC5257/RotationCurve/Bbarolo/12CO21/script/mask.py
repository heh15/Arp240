# import cube
from astropy.utils import data
from spectral_cube import SpectralCube, Projection
from radio_beam import Beam
from astropy import units as u
import time
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from regions import read_ds9

from matplotlib import rcParams
rcParams['mathtext.default']='regular'

############################################################
# basic information

position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')

############################################################
# function

def cut_2d(data,position,size,wcs):
    cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

def Regmask3d(data,region_pix,lowchan,highchan):
    region_mask=region_pix.to_mask()
    shape=np.shape(data)
    mask=region_mask.to_image(shape=((shape[1],shape[2])))
    mask3d=np.zeros((shape[0],shape[1],shape[2]))
    mask3d[lowchan:highchan]=mask
    maskTF=mask3d==1

    data_masked=np.copy(data)
    data_masked[maskTF]='nan'

    return data_masked

############################################################
# main program

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/RotationCurve/Bbarolo/12CO21/'
imageDir=Dir+'image/'
picDir=Dir+'picture/'
regionDir=Dir+'region/'

fitsfile=imageDir+'NGC5257_12CO21_pbcor_cube.fits'
kcube=SpectralCube.read(fitsfile)
wcs=WCS(kcube.hdulist[0].header).celestial
data=kcube.hdulist[0].data

### south west region
region=read_ds9(regionDir+'masks.reg')[0]
region_pix=region.to_pixel(wcs)
region_mask=region_pix.to_mask()
lowchan=30;highchan=51
data_masked=Regmask3d(data,region_pix,lowchan,highchan)
data=np.copy(data_masked)

### south continuum source
region=read_ds9(regionDir+'masks.reg')[1]
region_pix=region.to_pixel(wcs)
region_mask=region_pix.to_mask()
lowchan=29;highchan=49
data_masked=Regmask3d(data,region_pix,lowchan,highchan)
data=np.copy(data_masked)

### north spiral arms 
region=read_ds9(regionDir+'masks.reg')[2]
region_pix=region.to_pixel(wcs)
region_mask=region_pix.to_mask()
lowchan=32;highchan=62
data_masked=Regmask3d(data,region_pix,lowchan,highchan)
data=np.copy(data_masked)

fig=plt.figure()
ax=plt.subplot('111',projection=wcs)
plt.imshow(data_masked[40], origin='lower')

outputfits=imageDir+'NGC5257_12CO21_pbcor_cube_masked.fits'
hdu=fits.PrimaryHDU(data_masked)
hdu.header=kcube.hdulist[0].header
hdu.writeto(outputfits,overwrite=True)

