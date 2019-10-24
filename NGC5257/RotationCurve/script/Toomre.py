'''
Mar. 27th,2019. Hao He

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


############################################################
# directory
Dir='/1/home/heh15/workingspace/Arp240/NGC5257/RotationCurve/'
CO10Dir=Dir+'12CO10/'
CO21Dir=Dir+'12CO21/'
picDir=Dir+'picture/'
logDir=Dir+'log/'
imageDir=Dir+'image/'
mapDir=Dir+'map/'

############################################################
# basic settings
incl=0.42
PA=111.68
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

def Apmask_convert(aperture,data_masked):
    data_cut=data_masked.data
    data_mask=data_masked.mask
    apmask=aperture.to_mask(method='center')[0]
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    mask_tmp=mask==0
    ap_mask=np.ma.mask_or(mask_tmp,data_mask)
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

############################################################
bin=10

###  import the velocity disperstion data

fitsimage=imageDir+'NGC5257_12CO21_pbcor_cube_mom2.fits'
wcs=fits_import(fitsimage)[0]
data_masked=fits_import(fitsimage)[1]
mom2=data_masked.data
threshold=mom2<10
mom2[threshold]='nan'

## import the surface density data
fitsimage=imageDir+'NGC5257_12CO21_pbcor_cube_mom0.fits'
data_masked=fits_import(fitsimage)[1]
mom0=data_masked.data

mom0_binned=np.nanmean(np.nanmean(mom0.reshape(int(960/bin),bin,int(960/bin),bin),axis=-1),axis=1)
# threshold=mom0<5*rms_mom0
# mom0[threshold]='nan'

sds=mom0/(0.0109*beammaj*beammin*(freq/115.27)**2)*alpha/ratio


## calculate the epicycle of the receding side of the regions. 
# data2=np.loadtxt(logDir+'rotationcurve_rec.txt')
# data2=data2.transpose()
# size=24
# v_rot=data2[1][0:size]
# R_arcsec=data2[0][0:size]
# R=data2[0][0:size]*0.48
# dlnv=np.ediff1d(np.log(data2[1]));dlnv=dlnv[0:size]
# dlnr=np.ediff1d(np.log(data2[0]));dlnr=dlnr[0:size]

step=1
R_arcsec=np.linspace(0,30,int(30/step)+1);R_arcsec=R_arcsec[1:]
R_kpc=R_arcsec*0.48
A=290;B=8*0.48
v_rot=A*(1-np.exp(-R_kpc/B))
dvdR=A/B*np.exp(-R_kpc/B)
Gradln=R_kpc/v_rot*dvdR

v_rot_std=v_rot*1000
G=6.67e-11


### Calculate the azimuthal angle and radius of each pixel

## create array of coordinates 
size2=960;center_pix=[480,482]
x=range(size2)
y=range(size2)
index1=np.meshgrid(x,y)[0]
index2=np.meshgrid(x,y)[1]
vectors=np.stack([index1,index2],axis=2)
vectors[:,:,0]=vectors[:,:,0]-center_pix[0]
vectors[:,:,1]=vectors[:,:,1]-center_pix[1]
dists=np.sqrt(vectors[:,:,0]**2+vectors[:,:,1]**2)

## vector with certain PA
revangle=math.radians((360-PA-90))
majaxis=[math.cos(math.radians(PA+90)),math.sin(math.radians(PA+90))]

## transform the image to new coordinates.
transform=[[math.cos(revangle),-math.sin(revangle)],[math.sin(revangle),math.cos(revangle)]]
# transform=np.transpose(transform)
vec_transform=np.matmul(vectors[:,:],transform)
vec_transform[:,:,0]=vec_transform[:,:,0]/dists
vec_transform[:,:,1]=vec_transform[:,:,1]/dists
distinction=vec_transform[:,:,1]
mask=distinction<0

## the angle between pixel and major axis
cosine=np.dot(vectors[:,:], majaxis)/dists
angle=arccos(cosine)
angle=angle*180/math.pi
angle_binned=np.nanmean(np.nanmean(angle.reshape(int(960/bin),bin,int(960/bin),bin),axis=-1),axis=1)
# angle[mask]=2*math.pi-angle[mask]

# mask_array=np.zeros((960,960),dtype=bool)
# mask_array[0:40,0:40]=True
# cosine_masked=np.ma.masked_where(~mask_array,cosine)

# fig=plt.figure()
# im=plt.imshow(cosine,origin='lower')
# cbar=plt.colorbar(im)

### Calculate the epicycle array

## return the value
dists=dists*0.1*0.48
R=np.sqrt(dists**2*(cosine**2+(1-cosine**2)/incl**2))
R_binned=np.nanmean(np.nanmean(R.reshape(int(960/bin),bin,int(960/bin),bin),axis=-1),axis=1)

R_std=R*1000*3.1*10**16

epsi_tmp=math.sqrt(2)*v_rot_std*(1+Gradln)**0.5

# store the epicycle into the fits file

# draw the mask
epsi_array_tmp=np.zeros((np.shape(mom0)[0],np.shape(mom0)[1]))
epsi_array_masked=np.ma.masked_where(np.ma.nomask,epsi_array_tmp)

a=R_arcsec[0]
b=R_arcsec[0]*incl
center_sky=SkyEllipticalAperture(center,a=a*u.arcsec,b=b*u.arcsec,theta=PA*u.degree)
center_pix=center_sky.to_pixel(wcs=wcs)
center_mask=Apmask_convert(center_pix,epsi_array_masked)

rings=dict.fromkeys((range(len(R_arcsec)-1)))
rings_mask=dict.fromkeys((range(len(R_arcsec)-1)))

for i in range(len(rings)):
    rings[i]=aperture_ring(R_arcsec[i],R_arcsec[i+1],wcs)
    rings_mask[i]=Apmask_convert(rings[i],epsi_array_masked)

# Check the aperture 
# fig=plt.figure()
# plt.imshow(data_masked,origin='lower')
# center_pix.plot(color='blue')
# for i in range(len(data2[0])-1):
#     rings[i].plot(color='red')

# fig=plt.figure()
# plt.imshow(rings_mask[3],origin='lower')

epsi_array_tmp[~center_mask.mask]=epsi_tmp[0]

for i in range(len(rings)):
    epsi_array_tmp[~rings_mask[i].mask]=epsi_tmp[i+1]

mask=epsi_array_tmp==0
epsi_array_tmp[mask]='nan'

epsi_array=epsi_array_tmp/R_std

# radius_array=np.zeros((np.shape(mom1)[0],np.shape(mom1)[1]))
# radius_array[~center_mask.mask]=R[0]

# for i in range(size-1):
#     radius_array[~rings_mask[i].mask]=R[i+1]


# outputfits=mapDir+'NGC5257_12CO21_epsi.fits'
# hdul=fits.open(fitsfile)
# hdul[0].data=epsi_array
# hdul.writeto(outputfits)
# hdul.close()

### Calculate the Toomre factor

sds_std=sds*2e30/(3.1e16)**2
dispersions_std=mom2*1000
G=6.67e-11

Q=dispersions_std*epsi_array/(math.pi*G*sds_std)

## Correct the factor with stellar mass
fitsfile=imageDir+'mass_map_regrid.fits'
sd_star=fits_import(fitsfile)[1]
sdsstar_std=sd_star*2e30/(3.1e16)**2

# calculate the stellar velocity dispersion. 
length=scale_length;length_std=length*1000*3.1e16
vdstar_std=np.sqrt(2*math.pi*G*length_std/7.3)*sdsstar_std**0.5
vdstar=vdstar_std/1000
vdstar=vdstar.data
mask=np.ma.masked_invalid(mom2).mask
vdstar[mask]='nan'

Q_tot=Q*(1+sd_star/sds*mom2/vdstar)**(-1)
Qtot_binned=np.nanmean(np.nanmean(Q_tot.reshape(int(960/bin),bin,int(960/bin),bin),axis=-1),axis=1)

SFR=fits.getdata(imageDir+'NGC5257_33GHz_pbcor_regrid.fits')[0][0]
SFR_binned=np.nanmean(np.nanmean(SFR.reshape(int(960/bin),bin,int(960/bin),bin),axis=-1),axis=1)
beammajor=0.7; beamminor=0.58; freq=33
SFR_phy_binned=ssfr_radio(SFR_binned, beammajor, beamminor, freq, D)
wcs_sfr=WCS(fits.getheader(imageDir+'NGC5257_33GHz_pbcor_regrid.fits'))
wcs_sfr=wcs_sfr.celestial
levels=[4*1.3e-5]

fig=plt.figure()
sc=plt.scatter(R_binned,Qtot_binned,c=mom0_binned,marker='.', cmap='brg_r')
cbar=plt.colorbar(sc)
cbar.set_label('12CO 2-1 intensity $Jy km s^{-1} beam^{-1}$', fontsize=20)
plt.xlabel('Radius (kpc)', fontsize=20)
plt.ylabel('Toomre factor Q', fontsize=20)
fig.tight_layout()
plt.savefig(picDir+'NGC5257_Toomre_scatter.png')
           
# position=SkyCoord(dec=49.8583*u.arcmin,ra=204.9903*u.degree,frame='icrs')
# size=u.Quantity((54,42),u.arcsec)
# cut=Cutout2D(data=Q_tot,position=position,size=size,wcs=wcs)
# Qtot_cut=cut.data
import matplotlib.colors as colors

fig=plt.figure()

position=center; size=u.Quantity((36,36),u.arcsec)
Qtot_cut_wcs=cut_2d(Q_tot, position, size, wcs)[0]
Qtot_cut=cut_2d(Q_tot, position, size, wcs)[1]
# Qtot_cutbin=np.nanmean(np.nanmean(Qtot_cut.reshape(int(540/bin),bin,int(420/bin),bin),axis=-1),axis=1)
SFR_cut=cut_2d(SFR, position, size, wcs)[1]

ax=plt.subplot(projection=Qtot_cut_wcs)
im=ax.imshow(Qtot_cut,origin='lower', vmax=np.nanmax(Qtot_binned))
# rings[11].plot()
# rings[13].plot()
cbar=plt.colorbar(im)
cbar.set_label('$ Q_{tot}$',fontsize=24)
ax.contour(SFR_cut,levels=levels,colors=['red'])
plt.savefig(picDir+'NGC5257_Toomre_map.png')

# fig=plt.figure()
# plt.subplot(projection=wcs)
# im=plt.imshow(R,origin='lower')
# cbar=plt.colorbar(im)


### plot the trend of the correction factor ###
correction=(1+sd_star/sds*mom2/vdstar)**(-1)
correction_binned=np.nanmean(np.nanmean(correction.reshape(int(960/bin),bin,int(960/bin),bin),axis=-1),axis=1)
Q_binned=np.nanmean(np.nanmean(Q.reshape(int(960/bin),bin,int(960/bin),bin),axis=-1),axis=1)
epsiarray_binned=np.nanmean(np.nanmean(epsi_array.reshape(int(960/bin),bin,int(960/bin),bin),axis=-1),axis=1)

fig=plt.figure()
plt.scatter(R_binned, correction_binned, marker='.')


### save the angular velocity, beta and Q ###

Qfile=mapDir+'NGC5257_Toomre_map.fits'
hdu=fits.PrimaryHDU(Q_tot.data)
hdu.writeto(Qfile, overwrite=True)


## save the beta.  

# draw the mask
Gradln_array=np.zeros((np.shape(mom0)[0],np.shape(mom0)[1]))
Gradln_masked=np.ma.masked_where(np.ma.nomask,Gradln_array)

a=R_arcsec[0]
b=R_arcsec[0]*incl
center_sky=SkyEllipticalAperture(center,a=a*u.arcsec,b=b*u.arcsec,theta=PA*u.degree)
center_pix=center_sky.to_pixel(wcs=wcs)
center_mask=Apmask_convert(center_pix,Gradln_masked)

rings=dict.fromkeys((range(len(R_arcsec)-1)))
rings_mask=dict.fromkeys((range(len(R_arcsec)-1)))

for i in range(len(rings)):
    rings[i]=aperture_ring(R_arcsec[i],R_arcsec[i+1],wcs)
    rings_mask[i]=Apmask_convert(rings[i],Gradln_masked)

# fig=plt.figure()
# plt.imshow(rings_mask[3],origin='lower')

Gradln_array[~center_mask.mask]=Gradln[0]

for i in range(len(rings)):
    Gradln_array[~rings_mask[i].mask]=Gradln[i+1]

mask=Gradln_array==0
Gradln_array[mask]='nan'

betafile=mapDir+'NGC5257_beta_map.fits'
hdu=fits.PrimaryHDU(Gradln_array.data)
hdu.writeto(betafile, overwrite=True)

## save the Omega. 

# draw the mask
vrot_array=np.zeros((np.shape(mom0)[0],np.shape(mom0)[1]))
vrot_masked=np.ma.masked_where(np.ma.nomask,vrot_array)

a=R_arcsec[0]
b=R_arcsec[0]*incl
center_sky=SkyEllipticalAperture(center,a=a*u.arcsec,b=b*u.arcsec,theta=PA*u.degree)
center_pix=center_sky.to_pixel(wcs=wcs)
center_mask=Apmask_convert(center_pix,vrot_masked)

rings=dict.fromkeys((range(len(R_arcsec)-1)))
rings_mask=dict.fromkeys((range(len(R_arcsec)-1)))

for i in range(len(rings)):
    rings[i]=aperture_ring(R_arcsec[i],R_arcsec[i+1],wcs)
    rings_mask[i]=Apmask_convert(rings[i],vrot_masked)

# fig=plt.figure()
# plt.imshow(rings_mask[3],origin='lower')

vrot_array[~center_mask.mask]=v_rot[0]

for i in range(len(rings)):
    vrot_array[~rings_mask[i].mask]=v_rot[i+1]

mask=vrot_array==0
vrot_array[mask]='nan'

Omega_array=vrot_array*1000/R_std*(3600*24*365*10**6)

Omegafile=mapDir+'NGC5257_Omega_map.fits'
hdu=fits.PrimaryHDU(Omega_array.data)
hdu.writeto(Omegafile, overwrite=True)