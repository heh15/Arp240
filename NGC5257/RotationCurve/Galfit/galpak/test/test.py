import numpy
import scipy
import astropy
import matplotlib
import galpak
import asciitable
from astropy.io import fits
from galpak import run
import time


# restoring beam is 2.004",1.605",-70.245 deg
ALMA_b7 = galpak.Instrument(psf=galpak.GaussianPointSpreadFunction(fwhm=2.004,pa=-70.245,ba=float(1.605/2.004)),lsf=galpak.GaussianLineSpreadFunction(fwhm=1.0))

'''
min_bounds = galpak.GalaxyParameters(radius=0.5,inclination=0.0, velocity_dispersion=1e-5, maximum_velocity=-1600.0,turnover_radius=1e-5)
max_bounds = galpak.GalaxyParameters(radius=10.0,inclination=90.0, velocity_dispersion=300.0, maximum_velocity=1600.0,turnover_radius=5.0)
'''
initial_params = galpak.GalaxyParameters(x=159,y=161,z=67)#,flux=149.0,radius=3.08,inclination=85.28,pa=-41.92,turnover_radius=0.012,maximum_velocity=350.0, velocity_dispersion=160.0)


t_start = time.time()

# redshift 0.0225
mydisk = galpak.DiskModel(flux_profile='gaussian',redshift=0.0225)
NGC5257 = run('../NGC5257_12CO10_combine_image.fits', instrument=ALMA_b7,max_iterations=int(1500), model=mydisk)

t_end = time.time()
t_tot = t_end-t_start
#tell me how long the run took
print 'run took: ' + str(int(t_tot/60.0)) + ' minutes'
