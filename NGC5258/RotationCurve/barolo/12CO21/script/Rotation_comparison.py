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

from matplotlib import rcParams
rcParams['mathtext.default']='regular'

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/RotationCurve/'
baroloDir=Dir+'barolo/12CO21/'
scriptDir=baroloDir+'script/'
imageDir=baroloDir+'image/'

############################################################
# basic setting. 

incl=58;cosi=0.55
pa=100
ra=(13*u.degree+39*u.arcmin+52.922*u.arcsec)*15
dec=50*u.arcmin+24.4*u.arcsec
center=SkyCoord(ra=ra,dec=dec,frame='icrs')

############################################################
# main function

def aperture_ring(radius_arcsec,wcs):
    a_in=radius_arcsec-1.5
    a_out=radius_arcsec
    b_out=a_out*incl
    ring_sky=SkyEllipticalAnnulus(position,a_in=a_in*u.arcsec,a_out=a_out*u.arcsec,b_out=b_out*u.arcsec,theta=PA*u.degree)
    ring_pix=ring_sky.to_pixel(wcs=wcs)
    return ring_pix

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask(method='center')[0]
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

############################################################
# main program

# with innital velocity to be 100 km/s
filename=baroloDir+'test1/run1/ringlog1.txt'
data1=pd.read_csv(filename,header=0,sep=r"\s*")

filename=baroloDir+'test1/run2/ringlog1.txt'
data2=pd.read_csv(filename,header=0,sep=r"\s*")

filename=baroloDir+'test1/run3/ringlog1.txt'
data3=pd.read_csv(filename,header=0,sep=r"\s*")


# with initial velocity to be 100 km/s, the inclination angle and position angle is not fixed.  

filename=baroloDir+'test2/run1/ringlog1.txt'
data1=pd.read_csv(filename,header=0,sep=r"\s*")

filename=baroloDir+'test2/run2/ringlog1.txt'
data2=pd.read_csv(filename,header=0,sep=r"\s*")

filename=baroloDir+'test2/run3/ringlog1.txt'
data3=pd.read_csv(filename,header=0,sep=r"\s*")

fig=plt.figure()
plt.errorbar(data1['RAD(arcs)'],data1['VROT(km/s)'],[-data1['E_VROT1'], data1['E_VROT2']], label='run 1')
plt.errorbar(data2['RAD(arcs)'],data2['VROT(km/s)'],[-data2['E_VROT1'], data2['E_VROT2']], label='run 2')
plt.errorbar(data3['RAD(arcs)'],data3['VROT(km/s)'],[-data3['E_VROT1'], data3['E_VROT2']], label='run 3')

plt.ylim(0,300)
plt.legend()
plt.xlabel('Radius (arcsec)')
plt.ylabel('Velocity (km/s)')
plt.savefig('../picture/NGC5258_barolo_nocut.png')

# cut the radius tobe 12 arcsec. 
filename=baroloDir+'test3/run1/ringlog1.txt'
data1=pd.read_csv(filename,header=0,sep=r"\s*")

filename=baroloDir+'test3/run2/ringlog1.txt'
data2=pd.read_csv(filename,header=0,sep=r"\s*")

filename=baroloDir+'test3/run3/ringlog1.txt'
data3=pd.read_csv(filename,header=0,sep=r"\s*")

fig=plt.figure()
plt.errorbar(data1['RAD(arcs)'],data1['VROT(km/s)'],[-data1['E_VROT1'], data1['E_VROT2']], label='run 1')
plt.errorbar(data2['RAD(arcs)'],data2['VROT(km/s)'],[-data2['E_VROT1'], data2['E_VROT2']], label='run 2')
plt.errorbar(data3['RAD(arcs)'],data3['VROT(km/s)'],[-data3['E_VROT1'], data3['E_VROT2']], label='run 3')

plt.ylim(0,300)
plt.legend()
plt.xlabel('radius (arcsec)')
plt.ylabel('velocity (km/s)')
plt.savefig('../picture/NGC5258_barolo_cut.png')


####### Comparison #######

# imcube=SpectralCube.read(imageDir+'NGC5257_12CO21_pbcor_cube_masked.fits')
# Imcube=imcube.with_spectral_unit(u.km / u.s, velocity_convention='radio')

# # #### compare with moment 2 map.
# seperation=data3['RAD(arcs)'][1]-data3['RAD(arcs)'][0]
# rings=dict.fromkeys(range(data3.shape[0]))
# rings_mask=dict.fromkeys(range(data3.shape[0]))
# size=data3.shape[0]

# mom2=Imcube.linewidth_sigma()
# wcs=WCS(mom2.hdu.header)
# radius_arc=data3['RAD(arcs)']+seperation/2

# center_sky= SkyEllipticalAperture(positions=center,a=radius_arc[0]*u.arcsec,b=radius_arc[0]*cosi*u.arcsec,theta=pa*u.degree)
# rings[0]=center_sky.to_pixel(wcs=wcs)
# rings_mask[0]=Apmask_convert(rings[0],np.array(mom2))

# def aperture_rings(a_in,a_out,wcs,cosi,pa):
#     ring_sky=SkyEllipticalAnnulus(positions=center,a_in=a_in*u.arcsec,a_out=a_out*u.arcsec,b_out=a_out*cosi*u.arcsec,theta=pa*u.degree)
#     ring_pix=ring_sky.to_pixel(wcs=wcs)
#     return ring_pix

# for i in range(1,size):
#     rings[i]=aperture_rings(radius_arc[i-1],radius_arc[i],wcs,cosi,pa)
#     rings_mask[i]=Apmask_convert(rings[i],np.array(mom2))

# meandisp=np.empty((0,0))

# for i in range(size):
#     meandisp=np.append(meandisp,np.nanmean(rings_mask[i]))

# # fig=plt.figure()
# # plt.plot(data3['RAD(arcs)'],data3['DISP(km/s)'])
# # plt.plot(data3['#RAD(Kpc)'],meandisp,label='moment2')
# plt.legend()

# mom0=imcube.moment(order=0)
# wcs=WCS(mom0.hdu.header)

# fig=plt.figure()
# ax=plt.subplot(projection=wcs)
# plt.imshow(np.array(mom0),origin='lower')

# incl=56.317;cosi=0.55
# pa=100;
# ra=(13*u.degree+39*u.arcmin+52.922*u.arcsec)*15
# dec=50*u.arcmin+24.4*u.arcsec
# position=SkyCoord(ra=ra,dec=dec,frame='icrs')
# R=8.5*u.arcsec
# ring=SkyEllipticalAperture(position,a=R, b=R*cosi,theta=pa*u.degree)
# ring_pix=ring.to_pixel(wcs=wcs)

# ring_pix.plot(color='red')


# #### Check the moment 2 map made by barolo

# def fits_import(fitsimage):
#     hdr = fits.open(fitsimage)[0].header
#     wcs = WCS(hdr).celestial
#     data=fits.open(fitsimage)[0].data
#     data=np.squeeze(data)
#     data_masked=np.ma.masked_invalid(data)

#     return wcs, data_masked

# fitsimage=baroloDir+'test3/run3/maps/NGC_5257_local_2mom.fits'
# wcs=fits_import(fitsimage)[0]
# mom2_model=fits_import(fitsimage)[1].data

# pas=data3['P.A.(deg)']
# cosis=np.cos(data3['INC(deg)']*math.pi/180)

# center_sky= SkyEllipticalAperture(positions=center,a=radius_arc[0]*u.arcsec,b=radius_arc[0]*cosis[0]*u.arcsec,theta=pas[0]*u.arcsec)
# rings[0]=center_sky.to_pixel(wcs=wcs)
# rings_mask[0]=Apmask_convert(rings[0],mom2_model)

# for i in range(1,size):
#     rings[i]=aperture_rings(radius_arc[i-1],radius_arc[i],wcs,cosis[i],pas[i])
#     rings_mask[i]=Apmask_convert(rings[i],mom2_model)

# meandisp=np.empty((0,0))

# for i in range(size):
#     meandisp=np.append(meandisp,np.nanmean(rings_mask[i]))

# fig=plt.figure()
# plt.imshow(mom2_model,origin='lower')
# rings[2].plot()


# fig=plt.figure()
# plt.plot(data3['RAD(arcs)'],data3['DISP(km/s)'])
# plt.plot(data3['RAD(arcs)'],meandisp)


## import the moment 1 map, 12 arcsec ring. 
a_in=12-0.5
a_out=12+0.5

ring_12=aperture_rings(radius_arc[i-1],radius_arc[i],wcs,cosis[i],pas[i])
