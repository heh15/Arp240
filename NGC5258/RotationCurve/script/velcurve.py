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
Dir='/1/home/heh15/workingspace/Arp240/NGC5258/RotationCurve/'
baroloDir=Dir+'barolo/'
diskfitDir=Dir+'diskfit/'
CO10Dir=Dir+'12CO10/'
CO21Dir=Dir+'12CO21/'
picDir=Dir+'Picture/'
logDir=Dir+'log/'
imageDir=Dir+'Image/'
fitsDir=Dir+'Fits/'

############################################################
# basic settings
incl=0.57
PA=213.3
ra=204.9904
dec=49.8473
center=SkyCoord(dec=dec*u.arcmin,ra=ra*u.degree,frame='icrs')
freq=230.54
D=99
beammaj=1.004
beammin=0.556
ratio=0.77
alpha=0.86

############################################################ 
# function 

# read the diskfile output
def readvels(infile):
    radius=[];velocity=[];error=[]
    line=infile.readline()
    words=line.split()
    marker=words[0]
    while (marker != "Fitted"):
        line=infile.readline()
        words=line.split()
        if words==[]:
            continue
        marker=words[0]
    line=infile.readline()
    line=infile.readline()
    line=infile.readline()
    words=line.split()
    while(words != []):
        radius.append(words[0])
        velocity.append(words[2])
        error.append(words[3])
        line=infile.readline()
        words=line.split()
    radius=0.3*np.array([float(l) for l in radius])
    velocity=np.array([float(l) for l in velocity])
    error=np.array([float(l) for l in error])
    return radius, velocity,error

def aperture_ring(radius_in,radius_out,wcs):
    a_in=radius_in
    a_out=radius_out
    b_out=a_out*incl
    ring_sky=SkyEllipticalAnnulus(center,a_in=a_in*u.arcsec,a_out=a_out*u.arcsec,b_out=b_out*u.arcsec,theta=PA*u.degree)
    ring_pix=ring_sky.to_pixel(wcs=wcs)
    return ring_pix

def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

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

############################################################
# main program

i=0
radius=[];velocity=[]
filename=logDir+'NGC5258_12CO10_vel_radcut.out'
with open (filename,'r') as Input:
    radius,velocity,error=readvels(Input)
radius_arcsec=radius
vel=velocity
error_CO10=error

i=0
radius=[];velocity=[]
# filename=logDir+'NGC5258_12CO21_vel_radcut.out'
filename=diskfitDir+'test/NGC5258_12CO21_vel.out'
with open (filename,'r') as Input:
    radius,velocity,error=readvels(Input)

radius_2d_CO21=radius/3
vel_CO21=velocity
error_CO21=error
radius21=radius_2d_CO21

# read the data from barolo fitting. 
filename=baroloDir+'12CO21/test3/run3/ringlog1.txt'
with open (filename,'r') as Input:
    header=Input.readline()
    header=header.split()
    data=np.loadtxt(Input,skiprows=0)
data=np.transpose(data)
radius_3dB=data[1]
vel_3dB=data[2]
error_upp_3dB=data[13]
error_low_3dB=-data[14]
# yerr_3dB=error_low_3dB


# read the data from the literature
data=np.loadtxt(logDir+'rotationcurve_app.txt')
data=data.transpose()

data2=np.loadtxt(logDir+'rotationcurve_rec.txt')
data2=data2.transpose()


fig=plt.figure()
CO10_2d=plt.errorbar(radius_arcsec,vel,error_CO10,color='red',marker='o',linestyle='None',label='DiskFit 12CO10')
CO10_3d=plt.errorbar(radius_3dB,vel_3dB,[error_low_3dB, error_upp_3dB], color='green',marker='o',label='barolo fitting')
CO21_2d=plt.errorbar(radius_2d_CO21,vel_CO21,error_CO21,color='blue',marker='o',linestyle='None',label='DiskFit 12CO21')
Ha_rec=plt.plot(data[0],data[1],color='black',marker='h',linestyle='none',label='Ha receding side')
Ha_app=plt.errorbar(data2[0],data2[1],color='black',marker='*',linestyle='none',label='Ha approaching side')
plt.legend([CO10_2d,CO21_2d],['CO10 Diskfit','CO21 Diskfit'],loc='upper left')
plt.legend(loc='lower right')
plt.xlabel('radius(arcsec)')
plt.ylabel('velocity(km/s)')
plt.show()
plt.ylim(0,400)
plt.savefig(picDir+'NGC5258_rotationalcurve.png')


# # calculate the tommere factor of the receding side. 
# size=14
# v_rot=data[1][0:size]
# R=data[0][0:size]*0.48
# dlnv=np.ediff1d(np.log(data[1]));dlnv=dlnv[0:size]
# dlnr=np.ediff1d(np.log(data[0]));dlnr=dlnr[0:size]
# epsi=math.sqrt(2)*v_rot/R*(1+dlnv/dlnr)**0.5

# # calculate the velocity dispersion of for the aperture rings
# fitsimage=imageDir+'NGC5258_12CO21_combine_noise45_mom2.fits'
# wcs=fits_import(fitsimage)[0]
# data_masked=fits_import(fitsimage)[1]
# mom2=data_masked.data

# a=data[0][0]
# b=data[0][0]*incl
# center_sky=SkyEllipticalAperture(center,a=a*u.arcsec,b=b*u.arcsec,theta=PA*u.degree)
# center_pix=center_sky.to_pixel(wcs=wcs)
# center_mask=Apmask_convert(center_pix,data_masked)

# dispersions=np.empty((0,0))
# dispersions=np.append(dispersions,np.ma.median(center_mask))

# rings=dict.fromkeys((range(len(data)-1)))
# rings_mask=dict.fromkeys((range(len(data)-1)))

# for i in range(len(data[0])-1):
#     rings[i]=aperture_ring(data[0][i],data[0][i+1],wcs)
#     rings_mask[i]=Apmask_convert(rings[i],data_masked)

# for i in range(len(data[0])-1):
#     dispersion=np.ma.median(rings_mask[i])
#     dispersions=np.append(dispersions,dispersion)

# dispersions=dispersions[0:size]

# # fig=plt.figure()
# # ax=plt.subplot(projection=wcs)
# # im=ax.imshow(mom2,cmap='rainbow',origin='lower',vmax=60)
# # rings[5].plot(color='red')
# # rings[9].plot(color='red')
# # cbar=fig.colorbar(im)
# # lon = ax.coords[0]
# # lat = ax.coords[1]
# # lon.set_major_formatter('hh:mm:ss')

# # calculate the surface density of different rings
# fitsimage=imageDir+'NGC5258_12CO21_combine_noise45_mom0.fits'
# data_masked=fits_import(fitsimage)[1]
# mom0=data_masked.data

# a=data[0][0]
# b=data[0][0]*incl
# center_sky=SkyEllipticalAperture(center,a=a*u.arcsec,b=b*u.arcsec,theta=PA*u.degree)
# center_pix=center_sky.to_pixel(wcs=wcs)
# center_mask=Apmask_convert(center_pix,data_masked)

# # areas=math.pi*np.ediff1d(data[0]**2)
# # areas=np.insert(areas,0,math.pi*data[0][0]**2)
# # areas=areas*480**2

# for i in range(len(data[0])-1):
#     rings[i]=aperture_ring(data[0][i],data[0][i+1],wcs)
#     rings_mask[i]=Apmask_convert(rings[i],data_masked)

# sds=np.empty((0,0))
# sd=np.ma.mean(center_mask)
# sds=np.append(sds,np.ma.mean(center_mask))

# for i in range(len(data[0])-1):
#     sd=np.ma.mean(rings_mask[i])
#     sds=np.append(sds,sd)

# sds=sds/(0.0109*beammaj*beammin*(freq/115.27)**2)*alpha/ratio
# sds=sds[0:size]

# fig=plt.figure()
# ax=plt.subplot(projection=wcs)
# im=ax.imshow(mom0,cmap='rainbow',origin='lower',vmax=6)
# rings[5].plot(color='red')
# rings[9].plot(color='red')
# cbar=fig.colorbar(im)
# lon = ax.coords[0]
# lat = ax.coords[1]
# lon.set_major_formatter('hh:mm:ss')
# plt.savefig(picDir+'Toorme_ring.png')

# # calculate the Toomere factor
# dispersions_std=dispersions*1000
# sds_std=sds*2e30/(3.1*10**16)**2
# R_std=R*1000*3.1*10**16
# v_rot_std=v_rot*1000
# G=6.67e-11

# epsi=math.sqrt(2)*v_rot_std/R_std*(1+dlnv/dlnr)**0.5
# Q=dispersions_std*epsi/(math.pi*G*sds_std)


# # plt.figure()
# # plt.scatter(R,Q)
# # plt.ylabel('Toomre factor Q')
# # plt.xlabel('radius (kpc)')
# # plt.ylim(0,4)
# # plt.savefig(picDir+'Toomre.png')

# # calculate the average stellar density
# fitsfile=fitsDir+'mass_map_regrid.fits'
# sd_star=fits_import(fitsfile)[1]
# center_mask=Apmask_convert(center_pix,sd_star)

# for i in range(len(data[0])-1):
#     rings[i]=aperture_ring(data[0][i],data[0][i+1],wcs)
#     rings_mask[i]=Apmask_convert(rings[i],sd_star)

# sdstar_rad=np.empty((0,0))
# sdstar_rad=np.append(sdstar_rad,np.ma.mean(center_mask))

# for i in range(len(R)-1):
#     sdstar_rad=np.append(sdstar_rad,np.ma.mean(rings_mask[i]))
    
# ############################################################
# # record the results

# output=np.transpose(np.vstack((R,Q,sds,sdstar_rad)))
# filename=logDir+'Toomre_avg.txt'
# np.savetxt(filename,output)
