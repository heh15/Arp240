'''
Mar. 22nd, 2018
casa 5.1.1

'''
import glob, os
import re
from distutils.dir_util import copy_tree

imagename='NGC5257_12CO10_noise42_mom0.fits'
imsmooth(imagename=imagename,
         kernel='gauss',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         targetres=True,
         outfile='NGC5257_12CO10_noise42_smooth.mom0')
