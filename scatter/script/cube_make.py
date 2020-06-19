import cube
from spectral_cube import SpectralCube, Projection
from astropy import units as u

Dir = '/home/heh15/workingspace/Arp240/scatter/'
imageDir = Dir + 'image/'

# NGC 5257
name=imageDir+'NGC5257/12CO21/NGC5257_12CO21_combine_smooth_pbcor.fits'
imagecube=SpectralCube.read(name)
Imcube=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)


## create rms cube
rmscube=cube.calc_noise_in_cube(Imcube)

# mask the the low value of rmscube. 
mask=rmscube<3.0e-3*u.Jy/u.beam
lowrms=rmscube.with_mask(~mask)
newrms=lowrms.with_fill_value(3.0e-3)

## find the signal of the cube. 
outcube=cube.find_signal_in_cube(Imcube,newrms,snr_hi=5)
kcube=outcube.to(u.K)
kcube.write(imageDir+'NGC5257/NGC5257_kcube_pbcor.fits')

# NGC 5258

name=imageDir+'NGC5258/NGC5258_12CO21_combine_smooth_pbcor.fits'
imagecube=SpectralCube.read(name)
Imcube=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)

## create rms cube
rmscube=cube.calc_noise_in_cube(Imcube)

# mask the the low value of rmscube. 
mask=rmscube<3.0e-3*u.Jy/u.beam
lowrms=rmscube.with_mask(~mask)
newrms=lowrms.with_fill_value(3.0e-3)

## find the signal of the cube. 
outcube=cube.find_signal_in_cube(Imcube,newrms,snr_hi=5)
kcube=outcube.to(u.K)
kcube.write(imageDir+'NGC5258/NGC5258_kcube_pbcor.fits')

