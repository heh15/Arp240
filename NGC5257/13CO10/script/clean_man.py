'''
Mar 7th
'''

import numpy as np
import os
import scipy.ndimage as sni
import sys

vis = '/home/heh15/practice/Arp240/' \
      'workingspace/CO13_10/NGC5257_CO13_12m_contsub.ms'
cleanDir = '/home/heh15/practice/Arp240/'\
             'workingspace/CO13_10/NGC5257/test_man/'
preName = cleanDir + 'NGC5257_CO13_10_12m_man'

dirtyname = '/home/heh15/practice/Arp240/' \
           '/workingspace/CO13_10/NGC5257/' \
           '/dirtyimage/NGC5257_CO13_10_12m_mthresh_20km_dirty.image/'
field = '0'
#phasecenter = 0
mode = 'velocity'
restfreq='107.78GHz'
width = '20km/s' 
nchan = 35 
start = '-300km/s' 
cell='0.3arcsec'  
imsize = [320,320]
weighting = 'briggs'
robust = 0.5
imagermode = 'mosaic'
cyclefactor = 1.0  # default value
# stuff related to auto-masking
rmsCRTF = '/1/home/heh15/practice/Arp240/workingspace/CO13_10/' \
          'NGC5257/emission_free.crtf'
rmsChan = '1'

# Additional parameter for tclean
specmode='cube'
outframe='BARY'
deconvolver='hogbom' 
gridder='mosaic' 
pblimit=0.2


# CLEAN output names
myImage = preName + '.image'
myFlux = preName + '.flux'
myMask = preName + '.mask'
myResidual = preName + '.residual'

rms = imstat(imagename=dirtyname,
             region=rmsCRTF, chans=rmsChan)['rms'][0]
thresh=2*rms

delmod(vis=vis)

tclean(vis=vis,
       imagename=preName,
       field=field,
       antenna=antenna,
       specmode=specmode,
       outframe=outframe,
       restfreq=restfreq,
       width=width,
       nchan=nchan,
       start=start,
       deconvolver=deconvolver,
       cell=cell,
       imsize=imsize,
       weighting=weighting,
       robust=robust,
       gridder=gridder,
       niter=100000,
       threshold=str(thresh)+'Jy/beam',
       restoringbeam=restoringbeam,
       cyclefactor=cyclefactor,
       pblimit=pblimit,
       interactive=True)
