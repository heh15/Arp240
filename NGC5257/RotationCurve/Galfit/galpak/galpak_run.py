import galpak
from galpak import GalPaK3D
from astropy.io import fits

my_psf=galpak.ImagePointSpreadFunction('NGC5257_12CO10_combine_psf_2D.fits')
gk = GalPaK3D('NGC5257_12CO10_combine_image.fits', seeing=1.0, instrument=galpak.ALMA())
gk.instrument.psf=my_psf
gk.run_mcmc(max_iterations=500)
gk.save('my_galpak_run')



ALMA_b7 = galpak.Instrument(psf=galpak.GaussianPointSpreadFunction(fwhm=0.749,pa=5.826,ba=float(0.665/0.749)),lsf=galpak.GaussianLineSpreadFunction(fwhm=1.0))
