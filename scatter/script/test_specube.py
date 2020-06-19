import cube
from astropy.utils import data
from spectral_cube import SpectralCube, Projection
from radio_beam import Beam
from astropy import units as u


Dir='/1/home/heh15/workingspace/Arp240/scatter/'
imageDir=Dir+'image/'

# import the image into the spectral cube

name=imageDir+'NGC5257/NGC5257_12CO21_combine_sinbeam_cube.fits'
imagecube=SpectralCube.read(name)

# newbeam=Beam(1.02*u.arcsec,0.55*u.arcsec,-64.5*u.deg)
# outcube=cube.convolve_cube(imagecube,newbeam)

rmscube=cube.calc_noise_in_cube(imagecube)
outcube=cube.find_signal_in_cube(imagecube,rmscube)
