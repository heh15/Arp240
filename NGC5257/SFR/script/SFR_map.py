'''
Oct. 25th, 2018

make the star formation rate map

'''
import numpy as np
import os,glob
import scipy.ndimage as sni
import sys
import re
import itertools
from shutil import copytree
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import make_lupton_rgb
# import aplpy
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.wcs import WCS
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

############################################################
Dir='/home/heh15/workingspace/Arp240/NGC5257/SFR/'
workDir=Dir+'SFR_map/'
logDir=Dir+'log/'
picDir=Dir+'picture/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
regionDir=Dir+'region/'
os.chdir(workDir)

############################################################
# basic setting (for both galaxies. )

ra= 204*u.degree+58*u.arcmin+50*u.arcsec
dec=50*u.arcmin+7*u.arcsec
center_point= SkyCoord(dec=dec,ra=ra,frame='icrs')
size=u.Quantity((125,125),u.arcsec)

# common pixelsize while regriding 
pix_size=1.6

ratio=1.2

############################################################
# function

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

def sfr_24um(f_mjsr):
    flux=f_mjsr*10**(-17)/(4.25*10**10)*2.45**2
    L_nu=4*math.pi*(100*10**6*3.086*10**18)**2*flux
    L=L_nu*1.2*10**13
    SFR_tot=10**(-42.69)*L
    
    return SFR_tot

def sfr_70um(f_her):
    flux_erg=f_her*10**(-23)
    L_nu=4*math.pi*(100*10**6*3.086*10**18)**2*flux_erg
    L=4.283*10**12*L_nu
    SFR_tot=10**(-43.23)*L

    return SFR_tot

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask(method='center')
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)

array_round=np.vectorize(round_sig)

def flux_aperture_get(data_cut,aperture):
    flux=aperture_photometry(data_cut,apertures=aperture)['aperture_sum'][0]

    return flux

############################################################
# main program

wavelength=['24um','70um','95GHz','33GHz']
value=['flux','uncertainty','SFR']
df=pd.DataFrame(index=wavelength,columns=value)

regions=['center', 'south', 'southarm']
wavelength=['24','70']
Coordinates=pd.DataFrame(index=regions,columns=wavelength)

# ## add the header to the image. 
# filename='herschel_70um.fits'
# hdul=fits.open('herschel_70um.fits')
# hdr=hdul[1].header

# hdr['BMAJ']=0.0015166666666666666
# hdr['BMIN']=0.0016000000000000000
# hdr['BPA']=0.00

# hdul.writeto(filename.strip('.fits')+'_header.fits')
# hdul.close()

# filename='galex_FUV.fits'
# hdul=fits.open(filename)
# hdr=hdul[0].header

# hdr['BMAJ']=4.3/3600
# hdr['BMIN']=4.3/3600
# hdr['BPA']=0.00

# hdul.writeto(filename.strip('.fits')+'_header.fits')
# hdul.close()

# filename='spitzer_24um.fits'
# hdul=fits.open(filename)
# hdr=hdul[0].header

# hdr['BMAJ']=6.0/3600
# hdr['BMIN']=6.0/3600
# hdr['BPA']=0.00

# hdul.writeto(filename.strip('.fits')+'_header.fits')
# hdul.close()

## execute the imsmooth.py

## add the data together. 

filename=imageDir+'herschel_70um.fits' # pixel size 1.6 arcesec
wcs_her=fits_import(filename, item=1)[0]
data_70um=fits_import(filename, item=1)[1]/1.6**2*0.3**2
intensity_70um=data_70um/0.3**2*(3600*180/math.pi)**2/10**6

herschel_wcs_cut=cut_2d(data_70um, center_point, size, wcs_her)[0]
herschel_cut=cut_2d(data_70um, center_point,size, wcs_her)[1]

filename=imageDir+'spitzer_24um_regrid70.fits' # regrid to herschel pixel size.
wcs_spi=fits_import(filename)[0]
# hdul=fits.open(filename)
# hdr=hdul[0].header
# wcs=WCS(hdr)
# data_24um=hdul[0].data
data_24um=fits_import(filename)[1]
data_24um=data_24um-48
# hdul.close()
Spitzer_wcs_cut=cut_2d(data_24um, center_point, size, wcs_spi)[0]
Spitzer_cut=cut_2d(data_24um, center_point,size, wcs_spi)[1]

# filename='galex_FUV_header_smooth_regrid.fits'
# # hdul=fits.open(filename)
# # hdr=hdul[0].header
# # wcs=WCS(hdr)
# # data_FUV=hdul[0].data
# # hdul.close()

# data_FUV=fits_import(filename)[1]
# data_FUV_tmp=data_FUV*1.4*10**(-15)*1528**2/(3*10**18)
# data_FUV_tmp=data_FUV_tmp*10**17/(1.5**2)*4.25*10**10

# # regrid the image and add the spitzer and FUV map together
# SFR=8.1*10**(-2)*data_FUV_tmp+3.2*10**-3*data_24um
# SFR_ob=3.2*10**-3*data_24um
# SFR_uob=8.1*10**(-2)*data_FUV_tmp

# ratio=intensity_70um/data_24um

## add the aperture
ra=204*u.degree+58*u.arcmin+15*u.arcsec
dec=50*u.arcmin+13*u.arcsec
south=SkyCoord(ra=ra, dec=dec, frame='icrs')
Coordinates['70']['south']=south
south_her=SkyCircularAperture(positions=Coordinates['70']['south'],r=3*1.24*u.arcsec).to_pixel(herschel_wcs_cut)

ra=204*u.degree+58*u.arcmin+13*u.arcsec
dec=50*u.arcmin+13*u.arcsec
south=SkyCoord(ra=ra, dec=dec, frame='icrs')
Coordinates['24']['south']=south
south_spi=SkyCircularAperture(positions=Coordinates['24']['south'],r=3*1.24*u.arcsec).to_pixel(herschel_wcs_cut)

ra=204*u.degree+58*u.arcmin+15*u.arcsec
dec=50*u.arcmin+25*u.arcsec
Coordinates['70']['center']=SkyCoord(ra=ra,dec=dec,frame='icrs')
center_her=SkyCircularAperture(positions=Coordinates['70']['center'],r=3*1.24*u.arcsec).to_pixel(herschel_wcs_cut)

ra=204*u.degree+58*u.arcmin+14*u.arcsec
dec=50*u.arcmin+23*u.arcsec
south=SkyCoord(ra=ra,dec=dec,frame='icrs')
Coordinates['24']['center']=SkyCoord(ra=ra,dec=dec,frame='icrs')
center_spi=SkyCircularAperture(positions=Coordinates['24']['center'],r=3*1.24*u.arcsec).to_pixel(herschel_wcs_cut)

ra=(13*u.degree+39*u.arcmin+57.14*u.arcsec)*15 ;dec=49*u.arcmin+44.309*u.arcsec
position=SkyCoord(dec=dec,ra=ra,frame='icrs')
southarm=SkyCircularAperture(positions=position,r=3*2*u.arcsec).to_pixel(herschel_wcs_cut)


## add the region file 
from regions import read_ds9
file=regionDir+'NGC5257_arm.reg'
arm_sky=read_ds9(file)[0]
armpix_spi=arm_sky.to_pixel(Spitzer_wcs_cut)
armpix_her=arm_sky.to_pixel(herschel_wcs_cut)

### plot the figure
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig=plt.figure()
ax=plt.subplot('111',projection=Spitzer_wcs_cut)
ax.tick_params(direction='in')
ax.title.set_text('Spitzer 24um')
im=ax.imshow(Spitzer_cut,origin='lower', cmap='gist_ncar_r')
center_spi.plot(color='black')
south_spi.plot(color='black')
southarm.plot(color='black')
armpix_spi.plot(color='black')
plt.savefig(picDir+'spitzer_24um.png',bbox_inches='tight',pad_inches=0.2)

fig=plt.figure()
ax1=plt.subplot('111',projection=herschel_wcs_cut)
ax1.title.set_text('Herschel 70um')
ax1.tick_params(direction='in')
im=ax1.imshow(herschel_cut,origin='lower', cmap='gist_ncar_r')
center_her.plot(color='black')
south_her.plot(color='black')
southarm.plot(color='black')
armpix_her.plot(color='black')
plt.savefig(picDir+'herschel_70um.png',bbox_inches='tight',pad_inches=0.2)

# fig=plt.figure()
# ax2=plt.subplot('111',projection=wcs_her)
# ax2.title.set_text('70um/24um ratio')
# im=ax2.imshow(ratio,origin='lower',vmin=0)
# # divider= make_axes_locatable(ax2)
# # cax=divider.append_axes('right',size='10%',pad=0.2)
# cbar=plt.colorbar(im,fraction=0.06,pad=0.1)
# plt.savefig(picDir+'70um_over_24um.png')
 
# fig=plt.figure()
# ax=plt.subplot('111',projection=wcs)
# im=ax.imshow(SFR_uob,origin='lower',vmax=0.02)
# plt.colorbar(im)

counts=33
flux=counts*1.4*10**(-15)*1528**2/(3*10**18)*0.68*10**(23)
flux_erg=flux/10**23
L_uv=4*math.pi*(100*10**6*3.086*10**18)**2*flux_erg
test=0.68*10**(-28)*L_uv
print(test)

pixel_area=2.45*2.45
pixel_sr=pixel_area/(60**2*180/math.pi)**2

f_mjsr=9503
flux=f_mjsr*10**(-17)/(4.25*10**10)*2.45**2
flux=f_mjsr*pixel_sr*10**6
#flux=1.34*10**(-23)
L_nu=4*math.pi*(100*10**6*3.086*10**18)**2*flux
L=L_nu*1.2*10**13
SFR_spitzer=10**(-42.69)*L

# center_sky=SkyCircularAperture(position,a=a*u.arcsec)
# center_pix=center_sky.to_pixel(wcs=wcs_cut)
# center_mask=Apmask_convert(center_pix,data_cut)


L_uvnu=L_uv*1.96*10**15+3.89*L
SFR=10**(-43.35)*L_uvnu


# import 12CO3-2 data

CO_data=fits_import('NGC5257co32_all.map40r.mom0.fits')[1]
wcs_co=fits_import('NGC5257co32_all.map40r.mom0.fits')[0]
 
# fig=plt.figure()
# ax=plt.subplot(projection=wcs)
# ax.imshow(data_70um,origin='lower')
# ax.contour(CO_data,transform=ax.get_transform(wcs_co))

# fig=plt.figure()
# ax=plt.subplot('121',projection=wcs_spi)
# ax.title.set_text('spitzer 24um')
# im=ax.imshow(SFR_ob,origin='lower',vmax=0.5,vmin=0)
# ax2=plt.subplot('122',projection=wcs_spi)
# ax2.imshow(SFR_uob,origin='lower',vmax=0.5,vmin=0)
# ax2.title.set_text('galex FUV')
# plt.draw()
# p0 = ax.get_position().get_points().flatten()
# p1 = ax2.get_position().get_points().flatten()
# ax_cbar = fig.add_axes([p0[0],0.05, p1[2]-p0[0], 0.05])
# ax_cbar.set_label('$M_{\solar}/(kpc^2) $'  )
# ax.imshow(data_70um,origin='lower')
# levels=[5.0,10.0,15.0,20.0]
# ax.contour(CO_data,transform=ax.get_transform(wcs_co),levels=levels)
# plt.colorbar(im, cax=ax_cbar, orientation='horizontal')
# plt.savefig(picDir+'SFR_map.png')

# fig=plt.figure()
# ax2=plt.subplot('111',projection=wcs_spi)
# im=ax2.imshow(SFR_uob,origin='lower')
# ax2.title.set_text('galex FUV')
# plt.colorbar(im)
# plt.savefig(picDir+'SFR_uv.png')


# continuum image 33GHz
fitsimage=imageDir+'ngc5257_Ka_c_r0.5_ms.fits'
wcs_33GHz=fits_import('ngc5257_Ka_c_r0.5_ms.fits')[0]
data_33GHz=fits_import('ngc5257_Ka_c_r0.5_ms.fits')[1]

L_33=(6.4*10**-4)*4*math.pi*(99*3.086*10**24)**2*10**(-23)
south_tot=10**-27*(2.18*33**(-0.1)+15.1*33**(-0.7))**(-1)*L_33
south_syn=6.64*10**(-29)*33**0.7*L_33
south_fre=4.6*10**(-28)*33**(0.1)*L_33

fitsimage='NGC5257_33GHz_smooth.fits'
wcs_33smo=fits_import('NGC5257_33GHz_smooth.fits')[0]
data_33smo=fits_import('NGC5257_33GHz_smooth.fits')[1]


# test the position of fitted beam for 33GHz data for smoothed beam
position=SkyCoord(dec=0.0146051786921*u.rad,ra=-2.7057736283*u.rad,frame='fk5')
south_sky=SkyEllipticalAperture(position,a=3*ratio*u.arcsec,b=3*ratio*u.arcsec,theta=0.0*u.degree)
south_pix=south_sky.to_pixel(wcs=wcs_33smo)
# fig=plt.figure()
# ax=plt.subplot(projection=wcs_33smo)
# ax.imshow(data_33smo,origin='lower')
# south_pix.plot()
# # plt.savefig(picDir+'33GHz_south.png')
south_33=south_pix


beamarea=1.1331*5.76*5.46
beamarea_pix=beamarea/(0.12**2)

data_cut=data_33smo
aperture=south_33
flux=flux_aperture_get(data_33smo,aperture)
flux_33=flux/beamarea_pix
df['flux']['33GHz']=flux_33

L_33=(flux_33)*4*math.pi*(99*3.086*10**24)**2*10**(-23)
south_tot=10**-27*(2.18*33**(-0.1)+15.1*33**(-0.7))**(-1)*L_33
south_syn=6.64*10**(-29)*33**0.7*L_33
south_fre=4.6*10**(-28)*33**(0.1)*L_33
df['SFR']['33GHz']=south_tot

peak=8.1e-4
rms=5.17e-5
error_33=rms/peak*south_tot
df['uncertainty']['33GHz']=flux_33*rms/peak

south=Apmask_convert(south_33,data_33smo)
test=np.ma.sum(south)
# # test the position of fitted beam for 33GHz data. 
# position=SkyCoord(dec=0.0146051786921*u.rad,ra=-2.7057736283*u.rad,frame='fk5')
# south_sky=SkyEllipticalAperture(position,a=1.13708546612*u.arcsec,b=1.01120920045*u.arcsec,theta=15.604037779*u.degree)
# south_pix=south_sky.to_pixel(wcs=wcs_33GHz)
# fig=plt.figure()
# ax=plt.subplot(projection=wcs)
# ax.imshow(data_33GHz,origin='lower')
# south_pix.plot()
# plt.savefig(picDir+'33GHz_south.png')


# # test the position in the fitted aperture in herschel 70um data
position=SkyCoord(dec=0.014607*u.rad,ra=-2.705771*u.rad,frame='fk5')
south_sky=SkyEllipticalAperture(position,a=3*ratio*u.arcsec,b=3*ratio*u.arcsec,theta=123.58*u.degree)
south_pix=south_sky.to_pixel(wcs=wcs_her)
# fig=plt.figure()
# ax=plt.subplot(projection=wcs_her)
# ax.imshow(data_70um,origin='lower')
# south_pix.plot()
# plt.savefig(picDir+'70um_south.png')
south_70=south_pix

# # test the position in the fitted aperture in herschel 24um data
position=SkyCoord(dec=0.014607*u.rad,ra=-2.705778*u.rad,frame='fk5')
south_sky=SkyEllipticalAperture(position,a=3*ratio*u.arcsec,b=3*ratio*u.arcsec,theta=123.58*u.degree)
south_pix=south_sky.to_pixel(wcs=wcs_spi)
# fig=plt.figure()
# ax=plt.subplot(projection=wcs_spi)
# ax.imshow(data_24um,origin='lower')
# south_pix.plot(color='red')
# plt.savefig(picDir+'24um_south.png')
south_24=south_pix


south=Apmask_convert(south_24,data_24um)
flux_south_24=np.ma.sum(south)/2.45**2*pix_size**2
flux_24=flux_south_24*10**6/(4.25*10**10)*2.45**2
df['flux']['24um']=flux_24

# the flux from manually drawn circle
# flux_south_24=1.297*10**3
SFR_south_spi=sfr_24um(flux_south_24)
df['SFR']['24um']=SFR_south_spi
peak=149.2-48
rms=48.75-48
error_24=rms/peak*SFR_south_spi
df['uncertainty']['24um']=flux_24*rms/peak

south=Apmask_convert(south_70,data_70um)
flux_south_70=np.ma.sum(south)

# flux_south_70=0.7554
SFR_south_her=sfr_70um(flux_south_70)
peak=6.06e-2
rms=3.33e-4
error_70=rms/peak*SFR_south_her
df['flux']['70um']=flux_south_70
df['uncertainty']['70um']=flux_south_70*rms/peak
df['SFR']['70um']=SFR_south_her

fitsfile='NGC5257_95GHz_smooth.fits'
data_95GHz=fits_import(fitsfile)[1]
wcs_95GHz=fits_import(fitsfile)[0]


position=SkyCoord(dec=0.014607*u.rad,ra=-2.705775*u.rad,frame='fk5')
south_sky=SkyEllipticalAperture(position,a=3*ratio*u.arcsec,b=3*ratio*u.arcsec,theta=123.58*u.degree)
south_pix=south_sky.to_pixel(wcs=wcs_95GHz)
# fig=plt.figure()
# ax=plt.subplot(projection=wcs_her)
# ax.imshow(data_95GHz,origin='lower')
# south_pix.plot()
# plt.title('95GHz')
# plt.savefig(picDir+'70um_south.png')
south_95GHz=south_pix

south=Apmask_convert(south_95GHz,data_95GHz)
flux_temp=np.ma.sum(south)
beamarea=1.1331*5.76*5.46
beamarea_pix=beamarea/(0.3**2)
flux_south_95=flux_temp/beamarea_pix

df['flux']['95GHz']=flux_south_95

L_95=(flux_south_95)*4*math.pi*(99*3.086*10**24)**2*10**(-23)
south_tot=10**-27*(2.18*95**(-0.1)+15.1*95**(-0.7))**(-1)*L_95
peak=6.36e-4
rms=6.6e-5
error_95=rms/peak*south_tot 
df['uncertainty']['95GHz']=flux_south_95*rms/peak
df['SFR']['95GHz']=south_tot

array=array_round(df.values)
output=pd.DataFrame(array,index=df.index,columns=df.columns)
with open(logDir+'south_phot.txt','w') as out:
    out.write('# from NGC 5257, SFR, SFR_map.py \n')
    out.write('\n')
    output.to_string(out)

