import numpy as np
import os
import scipy.ndimage as sni
import sys

vis = '/1/home/heh15/practice/Arp240/'\
      'workingspace/CO12_21/calibrated/'\
      'NGC5257_12m.ms'
cleanDir = '/home/heh15/practice/Arp240/'\
             'workingspace/CO12_21/continuum/'
preName = cleanDir + 'NGC5257_12CO21_cont_12m'
contspw='0:0~200,0:450~800,1,2:0~200,2:400~800,3:0~700'

tclean(vis=vis,
      imagename=preName,
      specmode='mfs',
      spw=contspw,
      outframe='BARY',
      cell='0.3arcsec',
      imsize=[320,320],
      weighting='briggs',
      robust=0.5,
      deconvolver='hogbom',
      gridder='mosaic',
      niter=0,
      cyclefactor=1.0,
      pblimit=0.2,
      interactive=False,
      chanchunks=-1,
      threshold='0.16mJy/beam')
