'''
Feb. 28th, 2020

'''

import numpy as np
import os
import scipy.ndimage as sni
import sys

############################################################
# directories
vis = '/home/heh15/Data/Arp240/' \
      'arp240-110GHz.ms'
cleanDir = '/home/heh15/workingspace/Arp240/NGC5257/'\
             '13CO10/findlines/'
preName = cleanDir + 'NGC5257_13CO'

############################################################
# basic settings 

field = '0'
phasecenter='J2000 13h39m52.922 0d50m24.1'
mode = 'velocity'
width = '40km/s' 
cell='0.3arcsec'  
imsize = [320,320]
weighting = 'briggs'
robust = 0.5
imagermode = 'mosaic'
cyclefactor = 1.0  # default value
# stuff related to auto-masking
rmsCRTF = '/1/home/heh15/workingspace/Arp240/13CO10/' \
          'NGC5257/emission_free_54.crtf'
rmsChan = '1~5'

# Additional parameter for tclean
specmode='cube'
outframe='BARY'
deconvolver='hogbom' 
gridder='mosaic' 
pblimit=0.2
spws=['1', '2', '3']

############################################################
# main program

delmod(vis=vis)

for spw in spws:
    tclean(vis=vis,
          imagename=preName+'_spw'+spw,
          phasecenter=phasecenter,
          field=field,
          specmode='cube',
          outframe='BARY',
          width=width,
          cell=cell,
          imsize=imsize,
          weighting=weighting,
          robust=robust,
          deconvolver='hogbom',
          gridder='mosaic',
          niter=None,
          threshold='0.0Jy',
          cyclefactor=cyclefactor,
          pblimit=0.2,
           interactive=False, 
           spw=spw)
