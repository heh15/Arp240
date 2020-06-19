import sys
sys.path.append('/home/heh15/workingspace/script')
from cmapNorms import CASApower

import time
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy.wcs import WCS
from photutils import SkyEllipticalAperture
from photutils import SkyEllipticalAnnulus
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D

from matplotlib import rcParams
rcParams['mathtext.default']='regular'

############################################################
#directory

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/RotationCurve/'
baroloDir=Dir+'Barolo/12CO21/'
scriptDir=baroloDir+'script/'
imageDir=baroloDir+'image/'

############################################################
# basic setting. 

incl=58;cosi=0.55
pa=100
ra=15*(13*u.degree+39*u.arcmin+52.921*u.arcsec)
dec=50*u.arcmin+24.332*u.arcsec
center=SkyCoord(dec=dec,ra=ra,frame='icrs')
size=u.Quantity((36,36),u.arcsec)

############################################################
# main function

def fits_import(fitsimage, item=0):
    hdr = fits.open(fitsimage)[item].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[item].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def cut_2d(data_masked,position,size,wcs):
    cut=Cutout2D(data=data_masked,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

############################################################
# main program

fitsimage='../test7/run2/maps/NGC_5257_2mom.fits'
wcs=fits_import(fitsimage)[0]
data=fits_import(fitsimage)[1]
wcs_cut=cut_2d(data,center, size, wcs)[0]
data_cut=cut_2d(data,center, size, wcs)[1]

fitsimage='../test7/run2/maps/NGC_5257_local_2mom.fits'
model=fits_import(fitsimage)[1]
model_cut=cut_2d(model,center, size, wcs)[1]

# draw rings
# a_in=4*u.arcsec
# a_out=8*u.arcsec
# b_out=a_out*0.5
# theta=204*u.deg

# ring_sky=SkyEllipticalAnnulus(positions=center, a_in=a_in, a_out=a_out, b_out=b_out, theta=theta)
# ring_pix=ring_sky.to_pixel(wcs_cut)

# phi=204; xmax=360; xmin=0; ymax=360; ymin=0;xcen=180;ycen=180
# x = np.arange(0,xmax-xmin,0.1) 
# y = np.tan(np.radians(phi-90))*(x-xcen)+ycen 
vmin=10;vmax=70

fig=plt.figure()
ax=plt.subplot(121)
ax.imshow(data_cut, origin='lower', cmap='hsv', vmin=vmin, vmax=vmax)
# ring_pix.plot()
ax.set_xticks([])
ax.set_yticks([])
# ax.plot(x,y,color='#808080',linewidth=2)
ax.set_ylim(0,360)
ax2=plt.subplot(122)
im=ax2.imshow(model_cut, origin='lower', cmap='hsv', vmin=vmin, vmax=vmax)
ax2.set_xticks([])
ax2.set_yticks([])
# ax2.plot(x,y,color='#808080',linewidth=2)
ax2.set_ylim(0,360)
plt.subplots_adjust(wspace=0, hspace=0)
cbar_ax=fig.add_axes([0.92, 0.25, 0.02, 0.5])
cbar=plt.colorbar(im, cax=cbar_ax)
plt.savefig('../picture/NGC5257_barolo_mom2_80km.png', bbox_inches='tight', pad_inches=0.5)

 
