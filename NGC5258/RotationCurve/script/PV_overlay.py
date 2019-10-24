from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import make_lupton_rgb
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from matplotlib.patches import Arrow
import math
import glob


############################################################
# directory
Dir='/home/heh15/workingspace/Arp240/NGC5258/RotationCurve/'
imageDir=Dir+'Image/'
diskfitDir=Dir+'diskfit/'
logDir=Dir+'log/'
picDir=Dir+'Picture/'

############################################################
# function
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
    radius=np.array([float(l) for l in radius])
    velocity=np.array([float(l) for l in velocity])
    error=np.array([float(l) for l in error])
    return radius, velocity,error

############################################################
# main program

# import the PV fits data
filename=imageDir+'NGC5258_12CO21_pv.fits'
hdul=fits.open(filename)
hdr=hdul[0].header
wcs=WCS(hdr)
data=hdul[0].data
hdul.close()

data=np.array(data)
data=np.swapaxes(data,0,2)
data=np.swapaxes(data,0,1)
data=data[:,:,0]
wcs.wcs.crpix[1]=30
wcs.wcs.crval[1]=0
temp=wcs.wcs.cdelt[1]
wcs.wcs.cdelt[1]=temp/1000
wcs=wcs.dropaxis(2)

# import the DiskFit result
filename=logDir+'NGC5258_12CO10_vel_radcut.out'
with open (filename,'r') as Input:
    radius,velocity,error=readvels(Input)
radius=0.3*radius;radius_n=-radius
radius=wcs.wcs.crpix[0]+radius/wcs.wcs.cdelt[0]
radius_n=wcs.wcs.crpix[0]+radius_n/wcs.wcs.cdelt[0]
vel_proj=velocity*math.sin(math.radians(55.0));vel_projn=-vel_proj
vel_proj=wcs.wcs.crpix[1]-vel_proj/wcs.wcs.cdelt[1]
vel_projn=wcs.wcs.crpix[1]-vel_projn/wcs.wcs.cdelt[1]
error_proj=error*math.sin(math.radians(55.0))
error_proj=error_proj/wcs.wcs.cdelt[1]
error_projn=error_proj


filename=logDir+'NGC5258_12CO21_vel_radcut.out'
with open (filename,'r') as Input:
    radius21,velocity,error21=readvels(Input)
radius21=0.1*radius21;radius21_n=-radius21

# radius in the image wcs system
radius21=wcs.wcs.crpix[0]+radius21/wcs.wcs.cdelt[0]
radius21_n=wcs.wcs.crpix[0]+radius21_n/wcs.wcs.cdelt[0]

# velocity projected on the image plane
vel21_proj=velocity*math.sin(math.radians(55.0));vel21_projn=-vel21_proj
# velocity in the image wcs system
vel21_proj=wcs.wcs.crpix[1]-vel21_proj/wcs.wcs.cdelt[1]
vel21_projn=wcs.wcs.crpix[1]-vel21_projn/wcs.wcs.cdelt[1]
error21_proj=error21*math.sin(math.radians(55.0))
error21_proj=error21_proj/wcs.wcs.cdelt[1]
error21_projn=error21_proj

# Draw the figure. 
fig=plt.figure()
ax=plt.subplot(111,projection=wcs)
ax.imshow(data,cmap='rainbow',origin='lower',aspect='auto')
ax.errorbar(radius,vel_proj,error_proj,color='red',label='12CO10 diskfit')
ax.errorbar(radius_n,vel_projn,error_proj,color='red')
ax.errorbar(radius21,vel21_proj,error21_proj,color='black',label='12CO21 diskfit')
ax.errorbar(radius21_n,vel21_projn,error21_projn,color='black')
ax.set_xlabel('offset (arcsec)')
ax.set_ylabel('line of sight velocity (km/s)')
ax.legend()
fig.savefig(picDir+'NGC5258_12CO21_pv.png')

