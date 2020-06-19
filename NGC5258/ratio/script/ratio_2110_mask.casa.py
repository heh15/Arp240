'''
Mar 26th. 2018
casa 5.1.1
'''

import numpy as np
import os,glob
import scipy.ndimage as sni
import sys
import re
import itertools
from shutil import copyfile
from shutil import copytree

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/ratio/'
imageDir=Dir+'image/'
workDir=Dir+'image/ratio/2110/uvtaper_mask/'
rmsCRTF=Dir+'region/emission_free.crtf'
rmsChan='1~10'
regionDir=Dir+'region/'

imageName={}
imageName['12CO10']=glob.glob(imageDir+'12CO10/NGC5258_12CO10_combine_contsub_uvrange.image')\
         [0].split('/')[-1].replace('.image','')
imageName['12CO21']=glob.glob(imageDir+'12CO21/NGC5258_12CO21_combine_uvtaper.image')\
         [0].split('/')[-1].replace('.image','')
Mask='whole.mask'

rationame='NGC5258_2110_ratio_uvtaper'
rms2={}
ratio=[]
niter=5
CO12_10_chan='10~60'

origDir=os.getcwd()
os.chdir(workDir)
rmtables('*')

############################################################
# smooth the image

imagenames=[imageDir+'12CO10/'+imageName['12CO10']+'.image',imageDir+'12CO21/'+imageName['12CO21']+'.image']
for imagename in imagenames:
    imagename_tmp=re.split('\_|\.',imagename)
    CO=[string for string in imagename_tmp\
        if re.match('.*CO.*',string)][1]
    imsmooth(imagename=imagename,
             kernel='gauss',
             major='2.186arcsec',
             minor='1.896arcsec',
             pa='-87.314deg',
             targetres=True,
             outfile=imageName[CO]+'_smooth.image')        
imageName['12CO10']=imageName['12CO10']+'_smooth'
imageName['12CO21']=imageName['12CO21']+'_smooth'


# measure the rms of the image.
imagenames=[imageName['12CO10']+'.image',imageName['12CO21']+'.image']
for imagename in imagenames:
    imagename_tmp=re.split('\_|\.',imagename)
    CO=[string for string in imagename_tmp\
        if re.match('.*CO.*',string)][0]
    rms2[CO]= imstat(imagename=imagename,
             region=rmsCRTF, chans=rmsChan)['rms'][0]

regionfile=regionDir+'whole.crtf'
makemask(inpimage=imageName['12CO10']+'.image', inpmask=regionfile, mode='copy', output='whole.mask')

# # for 12CO10, we would regrid the image. 
# imregrid(imagename=imageName['12CO10']+'.image',
#          template=imageName['12CO21']+'.image',
#          output=imageName['12CO10']+'_regrid21.image')
# imregrid(imagename=imageDir+imageMask['12CO10'],
#          template=imageName['12CO21']+'.image',
#          output=imageMask['12CO10'].replace('.mask','')+'_regrid21.mask')

# copytree(imageDir+imageMask['12CO21'],imageMask['12CO21'])
# imageName['12CO10']=imageName['12CO10']+'_regrid21'
# imageMask['12CO10']=imageMask['12CO10'].replace('.mask','')+'_regrid21.mask'
    
# make the moment map for regrid image. 

immoments(imagename=imageName['12CO10']+'.image',moments=0,\
          includepix=[2*rms2['12CO10'],100],mask=Mask,
              chans='10~60',outfile=imageName['12CO10']+'.image'+'.mom0')
immoments(imagename=imageName['12CO21']+'.image',moments=0,
          outfile=imageName['12CO21']+'_tmp0.image'+'.mom0',
          chans='10~60',
          stretch=True,
          includepix=[2*rms2['12CO21'],100],mask=Mask)

peak21=imstat(imagename=imageName['12CO21']+'_tmp0.image'+'.mom0')['max'][0]
peak10=imstat(imagename=imageName['12CO10']+'.image.mom0')['max'][0]
ratio.append(peak21/peak10)

rmsRatio=rms2['12CO21']/rms2['12CO10']

for i in range(niter):
    immoments(imagename=imageName['12CO21']+'.image',moments=0,\
              includepix=[2*ratio[i]*rms2['12CO10'],100],\
                  chans=CO12_10_chan,mask=Mask,
              outfile=imageName['12CO21']+'_tmp'+str(i+1)+'.image'+'.mom0')
    immath(imagename=[imageName['12CO21']+'_tmp'+str(i+1)+'.image'+'.mom0',
                      imageName['12CO10']+'.image'+'.mom0'],
           expr='IM0/IM1',
           outfile=rationame+'_tmp'+str(i+1)+'.image')
    ratio.append(imstat(imagename=rationame+'_tmp'+str(i+1)+'.image')\
                 ['mean'][0])

# put the contour of imageName['12CO21'] overlaid on the image of 12CO10
immoments(imagename=imageName['12CO21']+'.image',moments=0,
          outfile=imageName['12CO21']+'.image'+'.mom0',
          chans='10~60',
          stretch=True,
          includepix=[2*ratio[niter-1]*rms2['12CO10'],100],mask=Mask)
immath(imagename=[imageName['12CO21']+'.image'+'.mom0',imageName['12CO10']+'.image'+'.mom0'],
       expr='IM0/IM1',
       outfile=rationame+'.image')

# make a mask out of ratio map and apply it back to the moment map.
immath(imagename=[rationame+'.image',rationame+'.image'],
       outfile=rationame+'_mask.image',expr='(IM0/IM1)')
immath(imagename=[rationame+'_mask.image',
                  imageName['12CO10']+'.image.mom0'],
       outfile=imageName['12CO10']+'_masked.image.mom0',
       expr='(IM0*IM1)')
immath(imagename=[rationame+'_mask.image',
                  imageName['12CO21']+'.image.mom0'],
       outfile=imageName['12CO21']+'_masked.image.mom0',
       expr='(IM0*IM1)')
imageName['12CO10']=imageName['12CO10']+'_masked'
imageName['12CO21']=imageName['12CO21']+'_masked'

impbcor(imagename= imageName['12CO10']+'.image.mom0', pbimage=imageDir+'12CO10/NGC5258_12CO10_combine_contsub_uvrange_mom0.pb', outfile=imageName['12CO10']+'_pbcor.image.mom0')
impbcor(imagename= imageName['12CO21']+'.image.mom0', pbimage=imageDir+'12CO21/NGC5258_12CO21_combine_uvtaper_mom0.pb/', outfile=imageName['12CO21']+'_pbcor.image.mom0')

imageName['12CO10']=imageName['12CO10']+'_pbcor'
imageName['12CO21']=imageName['12CO21']+'_pbcor'

exportfits(imagename=imageName['12CO10']+'.image.mom0',fitsimage=imageName['12CO10']+'_mom0.fits', overwrite=True)
exportfits(imagename=imageName['12CO21']+'.image.mom0',fitsimage=imageName['12CO21']+'_mom0.fits', overwrite=True)

immath(imagename=[imageName['12CO21']+'.image.mom0', imageName['12CO10']+'.image.mom0'], expr='IM0/IM1', outfile=rationame+'_pbcor.image')

rationame=rationame+'_pbcor'
exportfits(imagename=rationame+'.image', fitsimage=rationame+'.fits', overwrite=True)

flux_12CO10=imstat(imageName['12CO10']+'.image.mom0')['flux'][0]
flux_12CO21=imstat(imageName['12CO21']+'.image.mom0')['flux'][0]
flux_ratio=flux_12CO21/flux_12CO10*112.78**2/225.46**2

os.chdir(origDir)
