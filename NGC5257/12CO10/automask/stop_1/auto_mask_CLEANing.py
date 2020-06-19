""" Jan. 15,2018
Tool to non-interactively image and CLEAN
visibilities with iteratively improved and
automatically generated CLEAN masks. Basic process is
    1. a dirty map is made
    2. you give this script a CASA region text file
       (CRTF) from the viewer to estimate the RMS in
       that dirty map
    3. you specify a stopping criteria as a multiple
       of that RMS
    4. masks are made where pixels are above that
       stopping criteria as well as being within the
       specified minimum primary beam region and not
       being smaller than a specified multiple of the
       synthesized beam area
    5. an initial CLEAN threshold is set to a
       fraction of the dirty map peak
    6. the map is cumulatively CLEANed with the
       latest mask and current CLEAN threshold
    7. the CLEAN threshold is halved
    8. repeat steps 6 and 7 until the CLEAN threshold
       is less than the specified stopping criteria
It must be run within CASA as it uses multiple tasks.

There are two ways of running this script.
    1. Edit the parameters under the 'initialize
       CLEAN parameters' and 'stuff related to
       auto-masking' comments in this file and then
       run
         $ casa --nologger --nogui -c auto_mask_CLEANing.py
       in a terminal.
    2. Start an instance of CASA (DO NOT include the
       --nologger and --nogui startup options), enter
       each variable definition for the user defined
       parameters at the CASA prompt and then run
       this script with execfile
         CASA <1>: execfile('auto_mask_CLEANing.py')
       The user defined parameters saved in this file
       will be completely ignored when running it
       like this.

This was originally written using CASA 4.6.0-REL
(r36590).
"""
"""
The stop is set as 1.5.
This file is used to deal with 12m, 7m and combine single field data.
This script is run in casa 4.7.1
Add the step to make the moment 0 map
"""
# if makemask path bug is ever fixed can get rid of
# directory changing back and forth along with
# os.path.basename stuff in makemask call

import numpy as np
import os
import scipy.ndimage as sni
import sys


# if running from outside of CASA CLI
if 'casapy.py' not in sys.argv[-1]:
    # intialize CLEAN parameters
    vis = '/home/heh15/practice/Arp240/' \
          'workingspace/automask/arp240_combine_CO.ms/'
    cleanDir = '/home/heh15/practice/Arp240/' \
                 'workingspace/automask/stop_1/combine/'
    preName = cleanDir + 'NGC_5257_CO10_combine_auto1'
    field = '0'
    antenna = ''
    phasecenter = None
    mode = 'velocity'
    restfreq='112.73GHz'
    width = '10km/s' 
    nchan = 120 
    start = '-600km/s' 
    psfmode = 'clark'
    cell='0.5arcsec'  
    imsize = [320,320]
    weighting = 'briggs'
    robust = 0.5
    imagermode = 'mosaic'
    uvrange=''
    uvtaper = False
    outertaper = []
    finalBeam = []
    cyclefactor = 1.5
    # stuff related to auto-masking
    stop = 1.5
    beamMin = 0.5
    rmsCRTF = '/1/home/heh15/practice/Arp240/workingspace/automask/' \
              'emission_free.crtf'
    rmsChan = '21'
    minpb=0.2
    
# CLEAN output names
myImage = preName + '.image'
myFlux = preName + '.flux'
myMask = preName + '.mask'
myResidual = preName + '.residual'

# change directory to workaround makemask failing
# when specifying paths
origDir = os.getcwd()
os.chdir(cleanDir)

# restore the calibrated data
delmod(vis=vis)

# make dirty image to build first mask from
clean(vis=vis,
      imagename=preName,
      field=field,
      antenna=antenna,
      phasecenter=phasecenter,
      mode=mode,
      restfreq=restfreq,
      width=width,
      nchan=nchan,
      start=start,
      psfmode=psfmode,
      cell=cell,
      imsize=imsize,
      weighting=weighting,
      robust=robust,
      imagermode=imagermode,
      uvrange=uvrange,
      uvtaper=uvtaper,
      outertaper=outertaper,
      niter=0,
      threshold='0.0Jy',
      cyclefactor=cyclefactor,
      minpb=minpb,
      interactive=False)

# tell user to go make a CRTF to estimate the RMS in
# the dirty map
if rmsCRTF == '':
    done = False
    rmsCRTF = raw_input('Go draw and save a CASA'
                        +'viewer region for '
                        +'estimating the RMS in the '
                        +'dirty image \n'+myImage
                        +'.\nEnter the full path to'
                        +'that CRTF and press '
                        +'return.\n')
    while not done:
        if not os.path.exists(rmsCRTF):
            rmsCRTF = raw_input(rmsCRTF+' does not'
                                +'exist, try again.'
                                +'\n')
        else:
            print rmsCRTF, 'successfully found.', \
                  'Continuing...'
        done = True
rms = imstat(imagename=myImage,
             region=rmsCRTF, chans=rmsChan)['rms'][0]

# grab some details from the dirty image
major = imhead(imagename=myImage, mode='get',
               hdkey='beammajor')['value']
minor = imhead(imagename=myImage, mode='get',
               hdkey='beamminor')['value']
pixelSize = float(cell.split('arcsec')[0])
num = major*minor*np.pi
den = 4.0*np.log(2.0)
beamArea = (num/den)/(pixelSize**2)
peak = imstat(imagename=myImage)['max'][0]
thresh = peak/4.0

# loop over CLEAN making masks automatically until
# we hit our threshold
n = 0
while (thresh >= stop*rms):
    # make masks based on threshold
    autoMask = preName + '_autoMask' + str(n)
    # get values for masking as numpy arrays
    ia.open(myFlux)
    pbVals = ia.getchunk()
    ia.close()
    ia.open(myResidual)
    resVals = ia.getchunk()
    cs = ia.coordsys()
    ia.close()
    if n == 0:
        mVals = np.zeros(resVals.shape)
    else:
        ia.open(myMask)
        mVals = ia.getchunk()
        ia.close()
    # create masking numpy array
    nextMask = np.zeros(mVals.shape)
    nextMask[resVals > thresh] = 1.0
    nextMask[mVals > 0] = 1.0
    nextMask[pbVals < minpb] = 0.0
    # remove masking regions that are too small
    labeled, j = sni.label(nextMask)
    myHistogram = sni.measurements.histogram(labeled,
                                             0, j+1,
                                             j+1)
    object_slices = sni.find_objects(labeled)
    threshold = beamArea*beamMin
    for i in range(j):
        if myHistogram[i+1] < threshold:
            nextMask[object_slices[i]] = 0.0
    # save masking array as CASA, floats-only, image
    im = ia.newimagefromshape(outfile=autoMask,
                          shape=list(nextMask.shape),
                              csys=cs.torecord())
    im.close()
    ia.open(autoMask)
    ia.putchunk(nextMask)
    ia.close()

    # clean with automatically generated mask
    if (thresh/2.0 < stop*rms
            and thresh < 1.05*stop*rms):
        restoringbeam = finalBeam
    else:
        restoringbeam = []
    os.system('rm -rf '+myMask)
    clean(vis=vis,
          imagename=preName,
          field=field,
          antenna=antenna,
          phasecenter=phasecenter,
          mode=mode,
          restfreq=restfreq,
          width=width,
          nchan=nchan,
          start=start,
          psfmode=psfmode,
          cell=cell,
          imsize=imsize,
          weighting=weighting,
          robust=robust,
          imagermode=imagermode,
          uvrange=uvrange,
          uvtaper=uvtaper,
          outertaper=outertaper,
          mask=autoMask,
          niter=10000,
          threshold=str(thresh)+'Jy/beam',
          restoringbeam=restoringbeam,
          cyclefactor=cyclefactor,
          minpb=minpb,
          interactive=False)
    thresh /= 2.0

    # if a little more than stopping threshold, run
    # with thresh=stop*rms
    if 2.0*thresh == stop*rms: break
    if (thresh < stop*rms
            and 2.0*thresh > 1.05*stop*rms):
        thresh = stop*rms
    n += 1

# make moment map
imagename=preName+'.image'
outname=imagename+'.mom0'
immoments(imagename=imagename,moments=0,outfile=outname)

# go back to where we started
os.chdir(origDir)
