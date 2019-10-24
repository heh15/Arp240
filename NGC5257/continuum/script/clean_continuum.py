import numpy as np
import os
import scipy.ndimage as sni
import sys

vis = '/home/heh15/practice/Arp240/' \
      'workingspace/CO12_10/calibrated/arp240_12m_CalibratedData.ms/'
cleanDir = '/home/heh15/practice/Arp240/'\
             'workingspace/CO12_10/continuum/'
maskDir = '/home/heh15/practice/Arp240/'\
          'workingspace/continuum/'
preName = cleanDir + 'NGC5257_12CO10_cont_12m_spw0'

tclean(vis=vis,
      imagename=preName,
      field='0',
      specmode='mfs',
      spw='0:0~400;700~900',
      outframe='BARY',
      cell='0.3arcsec',
      imsize=[320,320],
      weighting='briggs',
      robust=0.5,
      deconvolver='hogbom',
      gridder='mosaic',
      niter=100000,
      cyclefactor=1.0,
      pblimit=0.2,
      interactive=True,
      chanchunks=-1,
      threshold='0.1mJy/beam',
      mask=maskDir+'NGC5257_13CO.mask')
