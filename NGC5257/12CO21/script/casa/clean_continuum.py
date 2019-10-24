import numpy as np
import os
import scipy.ndimage as sni
import sys

vis = '/1/home/heh15/workingspace/Arp240/NGC5257/'\
      '12CO21/calibrated/NGC5257_combine.ms'
cleanDir = '/home/heh15/workingspace/Arp240/NGC5257/'\
             '12CO21/continuum/'
preName = cleanDir + 'NGC5257_12CO21_cont_combine'
contspw='0:0~200,0:450~800,1,2,3,'\
        '4:0~200,4:450~800,5,6,7,'\
        '8:0~200,8:450~800,9,10,11,'\
        '12:0~200,12:450~800,13,14,15'

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
      threshold='0.06mJy/beam')
