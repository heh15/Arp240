'''
Dec. 5th,2019. Hao He

Based on the equation from Kruijssen et al. 2014 equtation 9)
'''

import time
import matplotlib.pyplot as plt
import numpy as np
import math
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
from matplotlib import rcParams
import glob
rcParams['mathtext.default']='regular'

from astropy.wcs.utils import pixel_to_skycoord
from astropy.coordinates import match_coordinates_sky
import pandas as pd
from regions import read_ds9
from reproject import reproject_interp

############################################################
# directory
Dir='/1/home/heh15/workingspace/Arp240/NGC5257/RotationCurve/'
CO10Dir=Dir+'12CO10/'
CO21Dir=Dir+'12CO21/'
picDir=Dir+'picture/'
logDir=Dir+'log/'
imageDir=Dir+'image/'
mapDir=Dir+'map/'
regionDir=Dir+'region/'

regions=['center', 'south', 'arm', 'west']
quantities=['Mstar']

regionfiles=dict.fromkeys(regions)
regionfiles['center']=regionDir+'center_co2-1.reg'
regionfiles['south']=regionDir+'south_co2-1.reg'
regionfiles['arm']=regionDir+'arm_co2-1.reg'
regionfiles['west']=regionDir+'west_co2-1.reg'

regionobjects=dict.fromkeys(regions)

regionproperties=pd.DataFrame(columns=quantities, index=regions)

############################################################
# basic settings

galaxy='NGC5257'

incl=0.42
PA=95
ra=15*(13*u.degree+39*u.arcmin+52.922*u.arcsec)
dec=50*u.arcmin+24.296*u.arcsec
center=SkyCoord(dec=dec,ra=ra,frame='icrs')
freq=230.54
D=99
beammaj=1.004
beammin=0.556
ratio=0.85
alpha=0.8
rms=3e-3
rms_mom0=3e-3*10*np.sqrt(50)

sr_arcsec=(180/math.pi*60**2)**2
arcsec_pc=480

scale_length=4.2
A=290;B=8*0.48 # v_rot=A*(1-exp(-r/B))


############################################################
# function
def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def cut_2d(data_masked,position,size,wcs=None):
    cut=Cutout2D(data=data_masked,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

def aperture_ring(radius_in,radius_out,wcs):
    a_in=radius_in
    a_out=radius_out
    b_out=a_out*incl
    ring_sky=SkyEllipticalAnnulus(center,a_in=a_in*u.arcsec,a_out=a_out*u.arcsec,b_out=b_out*u.arcsec,theta=PA*u.degree)
    ring_pix=ring_sky.to_pixel(wcs=wcs)
    return ring_pix

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask(method='center')
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

'''calculating the SFR based on radio continuum.'''
def sfr_radio(flux,freq,d):
    '''
    flux in Jy
    freq in Hz
    d in Mpc
    '''
    f_erg=flux*1e-23
    L=4*math.pi*(d*1e6*3.086e18)**2*f_erg
    SFR_tot=1e-27*(2.18*freq**-0.1+15.1*freq**-0.7)**(-1)*L

    return SFR_tot

def ssfr_radio(intensity,beammajor,beamminor,freq,d):
    '''
    flux in Jy
    beammajor and beamminor in arcsec
    freq in Hz
    d in Mpc
    sSFR in msun/kpc^2/yr
    '''
    SFR_tot=sfr_radio(intensity,freq,d)
    sSFR=SFR_tot/(1.1331*beammajor*beamminor)/(4.85*d/1000)**2

    return sSFR

def star_mass(flux_36,flux_45):
    mass=math.pow(10,5.65)*np.power(flux_36,2.85)*np.power(flux_45,-1.85)*(D/0.05)**2*0.7
    return mass

def cut_2d(data_masked,position,size,wcs):
    cut=Cutout2D(data=data_masked,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

def reproj_binning(data, wcs, bin_num):
    
    map_in_shape=np.shape(data)
    nx_in, ny_in=map_in_shape
    nx_out=math.trunc(nx_in/bin_num);ny_out=math.trunc(ny_in/bin_num)
    xs,ys=np.meshgrid(np.arange(nx_out), np.arange(ny_out))
    wcs_out=wcs.deepcopy()
    wcs_out.wcs.crpix =[math.trunc(nx_out/2), math.trunc(ny_out/2)]
    wcs_out.wcs.cdelt=wcs.wcs.cdelt*bin_num
    wcs_out.wcs.ctype = ['RA---SIN', 'DEC--SIN']
    coords_out=pixel_to_skycoord(xs, ys, wcs_out)
    coords_out_flat=coords_out.flatten()
    pixel_labels_out = np.arange(xs.size)
    data_binned=np.zeros((nx_out, ny_out)).flatten()
    map_out_shape=(nx_out, ny_out)
    
    xs_in, ys_in = np.meshgrid(np.arange(nx_in), np.arange(ny_in))
    coords_in = pixel_to_skycoord(xs_in, ys_in, wcs)
    pixel_map_arr = np.full((nx_in, ny_in), np.nan).flatten()

    i_in=0
    npix_in = coords_in.flatten().size
    dra, ddec = np.zeros(npix_in), np.zeros(npix_in)
    i_out, d2d, d3d = match_coordinates_sky(coords_in.flatten(), coords_out_flat)
    dra, ddec = (coords_in.flatten()).spherical_offsets_to(
        coords_out_flat[i_out])
    dra = dra.arcsec
    ddec = ddec.arcsec

    good = (-0.5001 <= dra) & (dra < 0.5001) & (-0.5001 <= ddec) & (ddec < 0.5001)
    pixel_map_arr[good]=pixel_labels_out[i_out[good]]
    data_labeled=np.stack((data.flatten(),pixel_map_arr), axis=1)
    nan_index=np.where(np.isnan(data_labeled[:,1]))
    data_labeled=np.delete(data_labeled, nan_index,axis=0)
    data_labeled=data_labeled[np.argsort(data_labeled[:,1])]
    data_group=np.split(data_labeled[:,0], np.cumsum(np.unique(data_labeled[:,1], return_counts=True)[1])[:-1])
    for i in pixel_labels_out:
        data_binned[i]=np.nanmean(data_group[i])

    data_binned=data_binned.reshape(map_out_shape)
        
    return wcs_out, data_binned
    
def Regmask_convert(aperture,data_cut):
    apmask=aperture.to_mask()
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

############################################################
bin_num=10

###  import the velocity disperstion data

fitsimage=imageDir+'NGC5257_12CO21_pbcor_cube_mom2.fits'
wcs=fits_import(fitsimage)[0]
data_masked=fits_import(fitsimage)[1]
mom2=data_masked.data
threshold=mom2<10
mom2[threshold]='nan'
wcs_out,mom2_binned=reproj_binning(mom2, wcs, bin_num)

## import the surface density data
fitsimage=imageDir+'NGC5257_12CO21_pbcor_cube_mom0.fits'
data_masked=fits_import(fitsimage)[1]
mom0=data_masked.data

## moment 0 maps. 
mom0_binned=reproj_binning(mom0, wcs, bin_num)[1]
sds=mom0_binned/(0.0109*beammaj*beammin*(freq/115.27)**2)*alpha/ratio

fitsfile=imageDir+'mass_map_regrid.fits'
sd_star=fits_import(fitsfile)[1]
sdsstar_std=sd_star*2e30/(3.1e16)**2
sdstar_binned=reproj_binning(sd_star, wcs, bin_num)[1]

length=scale_length
Star_height=length/7.3*1000
star_den=sdstar_binned/Star_height/2
mask=star_den<0.01
star_den[mask]='nan'

radius=8/0.48
ring=aperture_ring(radius, (radius+1), wcs_out)
ring_mask=Apmask_convert(ring, star_den)
rho_compare=np.ma.median(ring_mask)

radius_eff=5
ring2=aperture_ring(radius_eff, (radius_eff+1), wcs_out)
ring2_mask=Apmask_convert(ring2, sdstar_binned)
sd_star_compare=np.ma.mean(ring2_mask)

## gas disk

sds_std=sds*2e30/(3.1e16)**2
dispersions_std=mom2_binned*1000
G=6.67e-11

Gas_height_std=dispersions_std**2/(math.pi*G*sds_std)
Gas_height=Gas_height_std/(3.1e16)
gas_den_masked=(sds/(2*Gas_height)+np.sqrt(sds**2/(4*Gas_height**2)+4*sds*star_den/(2*Gas_height)))/2
gas_den=gas_den_masked.data

ring_gas_mask=Apmask_convert(ring, gas_den)
rho_gas_compare=np.nanmean(ring_gas_mask)


## get the total fraction of the gas component in the system
rho_tot=gas_den+star_den
fraction=gas_den/(gas_den+star_den)
size=u.Quantity((54,42),u.arcsec)
fraction_cut=cut_2d(fraction,center,size,wcs_out)[1]
wcs_cut=cut_2d(fraction,center,size,wcs_out)[0]

for key in regionfiles.keys():
    regionobject=read_ds9(regionfiles[key])[0]
    regionobjects[key]=regionobject

fitsimage = imageDir + 'NGC5257_12CO21_combine_noise40_mom0.fits'
wcs, contour_masked =  fits_import(fitsimage)
contour = contour_masked.data
contour_cut, footprint = reproject_interp((contour_masked, wcs), wcs_cut, shape_out=np.shape(fraction_cut))

fig=plt.figure()
ax=plt.subplot('111', projection=wcs_cut)
ax.tick_params(direction='in', labelsize=8)
# ring.plot(color='red')
ax.contour(contour_cut, colors='orange', linewidths=0.8, levels=[1.1, 2.2])
for key in regionobjects.keys():
    regionpix = regionobjects[key].to_pixel(wcs_cut)
    regionpix.plot(color='red', linewidth=2.0)
plt.imshow(fraction_cut, origin='lower',vmax=0.5)
plt.xlabel('J2000 Right Ascension')
plt.ylabel('J2000 Declination')
plt.colorbar()
plt.title(r'$\rho_{gas}/(\rho_{gas}+\rho_{star})$')
plt.savefig(picDir+galaxy+'_gas_vol_fraction.png')

### calculate the median fraction for each region. 

for key in regionobjects.keys():
    region_pix=regionobjects[key].to_pixel(wcs_cut)
    region_mask=Regmask_convert(region_pix, fraction_cut)
    regionproperties['Mstar'][key]=np.ma.median(region_mask)

filename=logDir+galaxy+'_aperture_Mstar.csv'
regionproperties.to_csv(filename)

## save the fraction to the fitsfile. 
outfits=mapDir+'gas_vol_fraction.fits'
hdu=fits.PrimaryHDU(fraction)
header=wcs_out.to_header()
header.insert(0, 'SIMPLE')
header['SIMPLE']=True
hdu.header=header
hdu.writeto(outfits, overwrite=True)

## save the fraction to the fitsfile. 
outfits=mapDir+'rho_vol.fits'
hdu=fits.PrimaryHDU(rho_tot)
header=wcs_out.to_header()
header.insert(0, 'SIMPLE')
header['SIMPLE']=True
hdu.header=header
hdu.writeto(outfits, overwrite=True)

