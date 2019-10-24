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
import pandas as pd
from astropy.utils import data
from spectral_cube import SpectralCube, Projection
from radio_beam import Beam
from astropy import units as u
from regions import read_ds9

############################################################
# directory

Dir='/home/heh15/workingspace/Arp240/NGC5257/SFR/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
regionDir=Dir+'region/'
picDir=Dir+'picture/'

############################################################
# basic settings

center=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')

beam_area=1.019*0.522*1.1331
beam_area_pix=beam_area/0.01
regions=['center','hinge','south','rest','all']
values=['FWMH']
linewidth=pd.DataFrame(index=regions,columns=values)
rms=0.003*u.Jy/u.beam

D=99
ratio=0.77

############################################################
# function

def fits_import(fitsimage, item=0):
    hdr = fits.open(fitsimage)[item].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[item].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

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

    return data_masked

def cut_2d(data,position,size,wcs):
    cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

############################################################
# main program
filename=logDir+'NGC5258_cluster.txt'
info=pd.read_csv(filename, header=None, sep=r"\s*", skiprows=5)
coordinates=info.iloc[:,1:3]

positions=list();apertures=list()
for i in range(coordinates.shape[0]):
    position=SkyCoord(ra=coordinates[1][i]*u.degree,dec=coordinates[2][i]*u.degree, frame='fk5')
    position_icrs=position.transform_to('icrs')
    positions.append(position_icrs)
    aperture=SkyCircularAperture(positions=position_icrs, r=0.5*u.arcsec)
    apertures.append(aperture)

# cube
fitsimage=imageDir+'NGC5257_12CO21_pbcor_cube.fits'
imcube=SpectralCube.read(fitsimage)
Imcube=imcube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)
wcs=WCS(Imcube.hdu.header).celestial
data=Imcube.hdu.data

# moment 0 map
fitsimage=imageDir+'NGC5257_12CO21_pbcor_cube_smooth_co21_mom0.fits'
mom0_wcs=fits_import(fitsimage)[0]
mom0=fits_import(fitsimage)[1]

size=u.Quantity((54,42),u.arcsec)
wcs_cut=cut_2d(mom0,center,size,wcs)[0]
mom0_cut=cut_2d(mom0,center,size,wcs)[1]

# 33GHz map
fitsimage=imageDir+'NGC5257_33GHz_pbcor_smooth_co21.fits'
cont_wcs=fits_import(fitsimage)[0]
cont=fits_import(fitsimage)[1]

size=u.Quantity((54,42),u.arcsec)
cont_wcs_cut=cut_2d(cont,center,size,cont_wcs)[0]
cont_cut=cut_2d(cont,center,size,cont_wcs)[1]

# spitzer 3.6um image. 
fitsimage=imageDir+'spitzer_regrid_36um.fits'
ir_wcs=fits_import(fitsimage)[0]
ir=fits_import(fitsimage)[1]

### south continuum source
ra=15*(13*u.degree+39*u.arcmin+52.939*u.arcsec)
dec=50*u.arcmin+12.473*u.arcsec
position=SkyCoord(ra=ra, dec=dec,frame='icrs') 

south_sky=SkyCircularAperture(positions=position, r=1.1*u.arcsec)
south_pix=south_sky.to_pixel(wcs)
south_pix_cut=south_sky.to_pixel(wcs_cut)


highchan=np.shape(data)[0]-1;lowchan=0
south_masked=Regmask3d(data,south_pix,lowchan,highchan)
south_mask=np.isnan(south_masked)

south=Imcube.with_mask(south_mask)
spectrum=south.mean(axis=(1,2))
# spectrum.quicklook()

specarray=spectrum.hdu.data
velocity=np.linspace(-300, 390, 70)

vel_mean=np.nansum(specarray*velocity)/np.nansum(specarray)
dispersion=math.sqrt(np.nansum(specarray*(velocity-vel_mean)**2)/np.nansum(specarray))  

## Leroy et al. 2018
l=1.1*D*4.85
mass_vir=892*l*dispersion**2


## 12CO 2-1 flux. 
flux=aperture_photometry(mom0, apertures=south_pix)['aperture_sum'][0]
flux=flux/beam_area*0.1**2
test=flux
flux_CO10=flux/4/ratio

mass_gas=1.05e4*(0.5/2)*flux_CO10*D**2*1.4


#### show the coordinates of the objects. 
fig=plt.figure()
ax=plt.subplot(projection=cont_wcs_cut)
plt.imshow(cont_cut,origin='lower', cmap='gist_ncar_r')
apertures_pix=list()
for aperture in apertures:
    aperture_pix=aperture.to_pixel(cont_wcs_cut)
    apertures_pix.append(aperture_pix)
    aperture_pix.plot(color='black', linewidth=2.0)
# ra=204.9692242*u.degree; dec=0.837503427*u.degree
# position=SkyCoord(ra=ra, dec=dec)
# test_sky=SkyCircularAperture(positions=position, r=0.5*u.arcsec)
# test_pix=test_sky.to_pixel(cont_wcs)
# test_pix.plot(color='red')
plt.savefig(picDir+'starcluster.png',bbox_inches='tight',pad_inches=0.2)

# calculate the offset between identified position and source peak. 
peak_position=position
linden_position=SkyCoord(ra=204.9703117*u.degree, dec=0.836690288*u.degree)

sep=position.separation(linden_position).arcsecond


##### the total south region ####

file=regionDir+'south_co21_gas.reg'
south_sky=read_ds9(file)[0]
south_pix2=south_sky.to_pixel(wcs)
south_pix2_cut=south_sky.to_pixel(wcs_cut)
highchan=np.shape(data)[0]-1;lowchan=0
south_masked=Regmask3d(data,south_pix2,lowchan,highchan)
south_mask=np.isnan(south_masked)

# 12CO 10 flux
flux=np.nansum(south_masked)
flux=flux/beam_area*0.1**2
flux_CO10=flux/4/ratio
mass_gas2=1.05e4*(0.5/2)*flux_CO10*D**2*1.4

south=Imcube.with_mask(south_mask)
spectrum=south.mean(axis=(1,2))
# spectrum.quicklook()
specarray=spectrum.hdu.data
velocity=np.linspace(-300, 390, 70)
vel_mean=np.nansum(specarray*velocity)/np.nansum(specarray)
dispersion=math.sqrt(np.nansum(specarray*(velocity-vel_mean)**2)/np.nansum(specarray))  


l2_tmp=np.sqrt(4.0*1.5)
l2=l2_tmp*D*4.85
mass_vir2=892*l*dispersion**2

fig=plt.figure()
ax=plt.subplot(projection=wcs_cut)
plt.imshow(mom0_cut, origin='lower', cmap='gist_ncar_r')
south_pix_cut.plot(color='red')
south_pix2_cut.plot(color='blue')
plt.savefig(picDir+'cluster_co21_mass.png',bbox_inches='tight',pad_inches=0.2)
