'''
Mar 26th. 2018
casa 5.1.1

June 25th, 2019
Use the image made with casa 5.4
'''

import numpy as np
import os,glob
import scipy.ndimage as sni
import sys
import re
import itertools

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/ratio/'
imageDir=Dir+'startImage/'
workDir=Dir+'2110/uvtaper/'
ImageDir=Dir+'image'
rmsCRTF=Dir+'region/emission_free.crtf'
rmsChan='60~69'

imageName={}
imageName['12CO10']=glob.glob(imageDir+'*12CO10*uv*.image')\
         [0].split('/')[-1].replace('.image','')
imageName['12CO21']=glob.glob(imageDir+'*12CO21*uv*.image')\
         [0].split('/')[-1].replace('.image','')

rationame='NGC5258_2110_ratio_uvtaper'
rms2={}
ratio=[]
niter=5
CO12_10_chan='10~60'

origDir=os.getcwd()
os.chdir(workDir)

############################################################
# smooth the image

imagenames=glob.glob(imageDir+'*12CO*uv*.image')
for imagename in imagenames:
    imagename_tmp=re.split('\_|\.',imagename)
    CO=[string for string in imagename_tmp\
        if re.match('.*CO.*',string)][0]
    imsmooth(imagename=imagename,
             kernel='gauss',
             major='2.186arcsec',
             minor='1.896arcsec',
             pa='-87.314deg',
             targetres=True,
             outfile=imageName[CO]+'_smooth.image')        
imageName['12CO10']=imageName['12CO10']+'_smooth'
imageName['12CO21']=imageName['12CO21']+'_smooth'

# # for 12CO10, we would regrid the image. 
# imregrid(imagename=imageName['12CO10']+'.image',
#          template=imageName['12CO21']+'.image',
#          output=imageName['12CO10']+'_regrid21.image')
# imageName['12CO10']=imageName['12CO10']+'_regrid21'

# measure the rms of the image.
imagenames=[imageName['12CO10']+'.image',imageName['12CO21']+'.image']
for imagename in imagenames:
    imagename_tmp=re.split('\_|\.',imagename)
    CO=[string for string in imagename_tmp\
        if re.match('.*CO.*',string)][0]
    rms2[CO]= imstat(imagename=imagename,
             region=rmsCRTF, chans=rmsChan)['rms'][0]
    
# make the moment map for regrid image. 

immoments(imagename=imageName['12CO10']+'.image',moments=0,\
          includepix=[5*rms2['12CO10'],100],
              chans='10~60',outfile=imageName['12CO10']+'.image'+'.mom0')
immoments(imagename=imageName['12CO21']+'.image',moments=0,
          outfile=imageName['12CO21']+'_tmp0.image'+'.mom0',
          chans='10~60',
          stretch=True,
          includepix=[5*rms2['12CO21'],100])

peak21=imstat(imagename=imageName['12CO21']+'_tmp0.image'+'.mom0')['max'][0]
peak10=imstat(imagename=imageName['12CO10']+'.image.mom0')['max'][0]
ratio.append(peak21/peak10)

rmsRatio=rms2['12CO21']/rms2['12CO10']

for i in range(niter):
    immoments(imagename=imageName['12CO21']+'.image',moments=0,\
              includepix=[5*ratio[i]*rms2['12CO10'],100],\
                  chans=CO12_10_chan,
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
          includepix=[5*ratio[niter-1]*rms2['12CO10'],100])
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

flux_12CO10=imstat(imageName['12CO10']+'.image.mom0')['flux'][0]
flux_12CO21=imstat(imageName['12CO21']+'.image.mom0')['flux'][0]
flux_ratio=flux_12CO21/flux_12CO10*107.78**2/225.46**2

# os.chdir(origDir)
