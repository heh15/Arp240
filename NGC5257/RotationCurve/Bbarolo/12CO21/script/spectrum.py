import cube
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
# directory

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/RotationCurve/'
baroloDir=Dir+'Bbarolo/12CO21/'
scriptDir=baroloDir+'script/'
imageDir=baroloDir+'image/'
picDir=Dir+'picture/'
regionDir=baroloDir+'region/'

############################################################
# basic settings

position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')
############################################################
# function

# mask 3d cube with a region file in specified channels (June 28th)
def Regmask3d(data,region_pix,lowchan,highchan):
    region_masks=region_pix.to_mask()
    if type(region_masks)==list:
        region_mask=region_masks[0]
    else:
        region_mask=region_masks
    shape=np.shape(data)
    mask=region_mask.to_image(shape=((shape[1],shape[2])))
    mask3d=np.zeros((shape[0],shape[1],shape[2]))
    mask3d[lowchan:highchan]=mask
    maskTF=mask3d==1

    data_masked=np.copy(data)
    data_masked[maskTF]='nan'

    return data_masked, maskTF

def fits_import(fitsimage, item=0):
    hdr = fits.open(fitsimage)[item].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[item].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def cut_2d(data,position,size,wcs):
    cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

############################################################
# extract the spectrum of the anomaly region. 

fitsimage=imageDir+'NGC5257_12CO21_pbcor_cube.fits'
imagecube=SpectralCube.read(fitsimage)
imcube=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)
wcs=WCS(imcube.hdu.header).celestial
data=imcube.hdu.data

file=regionDir+'anomaly.reg'
region_sky= read_ds9(file,errors='warn')[0]
region_pix= region_sky.to_pixel(wcs)

lowchan=0; highchan=69
mask=Regmask3d(data, region_pix, lowchan, highchan)[1]

regioncube=imcube.with_mask(mask)
spectrum=regioncube.mean(axis=(1,2))
# spectrum.quicklook()

fig=plt.figure()
ax=plt.subplot(111)
plt.plot(spectrum.spectral_axis, spectrum.value, drawstyle='steps-mid')
ax.set_xlabel('km/s', fontsize=20)
ax.set_ylabel('Jy/beam', fontsize=20)
ax.tick_params(labelsize=20)
fig.tight_layout()
plt.savefig(picDir+'region_spectral.png')

# test=regioncube.moment(order=0)
# test.quicklook()

fitsimage=imageDir+'NGC5257_12CO21_pbcor_cube_mom1.fits'
wcs=fits_import(fitsimage)[0]
mom1=fits_import(fitsimage)[1]

fig=plt.figure()

size=u.Quantity((54,42),u.arcsec)
wcs_cut=cut_2d(mom1,position,size,wcs)[0]
mom1_cut=cut_2d(mom1,position,size,wcs)[1]
region_pix_cut=region_sky.to_pixel(wcs_cut)

ax=plt.subplot(111, projection=wcs_cut)
im=plt.imshow(mom1_cut, origin='lower', cmap='rainbow')
region_pix_cut.plot(color='black')
cbar=plt.colorbar(im)
cbar.ax.tick_params(labelsize=20)
cbar.set_label('$Jy km \dot s^{-1} beam^{-1}$', fontsize=20)
fig.tight_layout()
plt.savefig(picDir+'region_map.png', bbox_inches='tight',pad_inches=0.2)
