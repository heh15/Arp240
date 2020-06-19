'''
calculate the ratio of 2110 based on 13CO 1-0 threshold

rms_12CO21=5.4e-3
rms_12CO10=6.3e-4

'''

import numpy as np
import os,glob
import scipy.ndimage as sni
import sys
import re
import itertools

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/ratio/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
workDir=imageDir+'ratio/2113/'
ratioCRTF=Dir+'region/majorpart.crtf'
rmsCRTF=Dir+'region/emission_free.crtf'


############################################################
# basic setting

imageName={}
imageName['12CO21']=glob.glob(imageDir+'12CO21/NGC5257_12CO21_combine_contsub_uvtaper.image')\
         [0].split('/')[-1].replace('.image','')
imageName['13CO10']=glob.glob(imageDir+'13CO10/*13CO10*_uvrange.image')\
         [0].split('/')[-1].replace('.image','')

rationame='NGC5257_2113_ratio'
rms1={}
ratios=[]
niter= 20
CO12_10_chan='10~60'

origDir=os.getcwd()
os.chdir(workDir)
rmtables('*')

############################################################
# smooth the image
imagenames=[imageDir+'12CO21/'+imageName['12CO21']+'.image',imageDir+'13CO10/'+imageName['13CO10']+'.image']
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
imageName['12CO21']=imageName['12CO21']+'_smooth'
imageName['13CO10']=imageName['13CO10']+'_smooth'


# measure the rms of the image.
imagenames=[imageName['12CO21']+'.image',imageName['13CO10']+'.image']
for imagename in imagenames:
    imagename_tmp=re.split('\_|\.',imagename)
    CO=[string for string in imagename_tmp\
        if re.match('.*CO.*',string)][0]
    rms1[CO]= imstat(imagename=imagename,
             region=rmsCRTF)['rms'][0]


# after measuring the rms, make the moment map for each image with cutoff 2*rms or 3*rms.
immoments(imagename=imageName['12CO21']+'.image',moments=0,
          includepix=[3*rms1['12CO21'],100],
              chans=CO12_10_chan,
          outfile=imageName['12CO21']+'_tmp0.image'+'.mom0')
immoments(imagename=imageName['13CO10']+'.image',moments=0,
          includepix=[2*rms1['13CO10'],100],
          chans='5~25',outfile=imageName['13CO10']+'.image'+'.mom0')
immath(imagename=[imageName['12CO21']+'_tmp0.image'+'.mom0',
                  imageName['13CO10']+'.image'+'.mom0'],
       expr='IM0/IM1',
       outfile=rationame+'_tmp0.image')
# ratios.append(imstat(imagename=rationame+'_tmp0.image',
#              region=ratioCRTF)['mean'][0])

peak12=imstat(imagename=imageName['12CO21']+'_tmp0.image'+'.mom0',region=ratioCRTF)['max'][0]
peak13=imstat(imagename=imageName['13CO10']+'.image.mom0',region=ratioCRTF)['max'][0]
peakRatio=peak12/peak13

ratios.append(35.0)


for i in range(niter):
    immoments(imagename=imageName['12CO21']+'.image',moments=0,\
              includepix=[2*ratios[i]*rms1['13CO10'],100],\
                  chans=CO12_10_chan,
              outfile=imageName['12CO21']+'_tmp'+str(i+1)+'.image'+'.mom0')
    immath(imagename=[imageName['12CO21']+'_tmp'+str(i+1)+'.image'+'.mom0',
                      imageName['13CO10']+'.image'+'.mom0'],
           expr='IM0/IM1',
           outfile=rationame+'_tmp'+str(i+1)+'.image')
    flux1=imstat(imagename=imageName['12CO21']+'_tmp'+str(i+1)+'.image'+'.mom0', region=ratioCRTF)['flux'][0]
    flux2=imstat(imageName['13CO10']+'.image'+'.mom0', region=ratioCRTF)['flux'][0]
    # ratios.append(imstat(imagename=rationame+'_tmp'+str(i+1)+'.image',
    #               region=ratioCRTF)['mean'][0])
    ratios.append(flux1/flux2)

# If around 10, then rename this ratio map as final image.
if os.path.isdir(rationame+'.image')==False:
    os.rename(rationame+'_tmp'+str(niter)+'.image',rationame+'.image')
    if os.path.isdir(imageName['12CO21']+'.image.mom0')==False:
        os.rename(imageName['12CO21']+'_tmp'+str(niter)+'.image'+'.mom0',
                  imageName['12CO21']+'.image.mom0')

immath(imagename=[rationame+'.image',rationame+'.image'],
       outfile=rationame+'_mask.image',expr='(IM0/IM1)')
immath(imagename=[rationame+'_mask.image',
                  imageName['13CO10']+'.image.mom0'],
       outfile=imageName['13CO10']+'_masked.image.mom0',
       expr='(IM0*IM1)')
immath(imagename=[rationame+'_mask.image',
                  imageName['12CO21']+'.image.mom0'],
       outfile=imageName['12CO21']+'_masked.image.mom0',
       expr='(IM0*IM1)')
imageName['12CO21']=imageName['12CO21']+'_masked'
imageName['13CO10']=imageName['13CO10']+'_masked'

impbcor(imagename= imageName['12CO21']+'.image.mom0', pbimage=imageDir+'12CO21/NGC5257_12CO21_combine_uvtaper_mom0.pb', outfile=imageName['12CO21']+'_pbcor.image.mom0')
impbcor(imagename= imageName['13CO10']+'.image.mom0', pbimage=imageDir+'13CO10/NGC5257_13CO10_12m_uvrange_mom0.pb', outfile=imageName['13CO10']+'_pbcor.image.mom0')

imageName['12CO21']=imageName['12CO21']+'_pbcor'
imageName['13CO10']=imageName['13CO10']+'_pbcor'

immath(imagename=[imageName['12CO21']+'.image.mom0', imageName['13CO10']+'.image.mom0'], expr='IM0/IM1', outfile=rationame+'_pbcor.image')

flux_12CO21=imstat(imageName['12CO21']+'.image.mom0')['flux'][0]
flux_13CO10=imstat(imageName['13CO10']+'.image.mom0')['flux'][0]
flux_ratio=flux_12CO21/flux_13CO10

exportfits(imagename=rationame+'.image',fitsimage=rationame+'.fits')

os.chdir(origDir)

