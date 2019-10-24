
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
    apmask=aperture.to_mask(mode='center')
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    mask_tmp=mask==0
    ap_mask=np.ma.mask_or(mask_tmp,data_mask)
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

# importfits(fitsimage='mass_map.fits',imagename='mass_map.image')

# ia.open('mass_map.image')
# ia.adddegaxes(spectral=True,stokes='Q',outfile='test.image')
# ia.close()

# image='mass_map_stokes.image'
# region='source.crtf'
# output='test.mask'
# makemask(inpimage=image,inpmask=region,mode='copy',output=output)

from regions import read_ds9

fitsimage='mass_map.fits'
hdr = fits.open(fitsimage)[0].header
wcs=fits_import(fitsimage)[0]
mass=fits_import(fitsimage)[1]
mass_data=mass.data


file='source.reg'
sources=read_ds9(file,errors='warn')
# source_sky=read_ds9(file,errors='warn')[1]
# source_pix=source_sky.to_pixel(wcs)
# source_masked=Apmask_convert(source_pix,mass)
# source_mask=~source_masked.mask

for source in sources:
    source_pix=source.to_pixel(wcs)
    source_masked=Apmask_convert(source_pix,mass)
    source_mask=~source_masked.mask
    mass_data[source_mask]='nan'
    

# fig=plt.figure()
# plt.subplot(projection=wcs)
# plt.imshow(source_masked,origin='lower')
# source_pix.plot(color='red')


fig=plt.figure()
plt.subplot(projection=wcs)
plt.imshow(mass_data,origin='lower')

outfits='mass_map_masked.fits'
hdu=fits.PrimaryHDU(mass_data)
hdu.header=hdr
hdu.writeto(outfits)
