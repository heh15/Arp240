import os
from subprocess import call
import re

############################################################
# basic settings 

Dir= '/1/home/heh15/workingspace/Arp240/NGC5257/ratio/'
workDir= Dir+'countchannel/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'

threshold_13= 2*0.00062
threshold_12= 2*0.0031
#threhold_12=2*14.2*0.00062

############################################################
# functions

def mask_shape(image,region,output):
    makemask(inpimage=image,inpmask=region,mode='copy',output=output)

def mask_threshold(image,lower,upper,output):
    ia.open(image)
    ia.calcmask('"'+image+'">'+str(lower)+' && '+'"'+image+'" <'+str(upper),name='mask0')
    ia.close()
    makemask(inpimage=image,inpmask=image+':mask0',mode='copy',output=output)
    makemask(mode='delete',inpmask=image+':mask0')

def mask_invert(mask,output):
    ia.open(image)
    ia.calcmask('"'+image+'"< 0.5')
    ia.close()
    makemask(inpimage=image,inpmask=image+':mask0',mode='copy',output=output)
    makemask(mode='delete',inpmask=mask+':mask0')

def mask_subtract(mask1,mask2,output):
    immath(imagename=[mask1,mask2],expr='IM0-IM1',outfile='temp.mask')
    ia.open('temp.mask')
    ia.calcmask('"temp.mask">0.5',name='mask0')
    ia.close()
    makemask(mode='copy', inpimage='temp.mask', inpmask='temp.mask:mask0', output=output)
    rmtables('temp.mask')

def mask_and(mask1,mask2,output):
    immath(imagename=[mask1,mask2],expr='IM0+IM1',outfile='temp.mask')
    ia.open('temp.mask')
    ia.calcmask('"temp.mask">1.5',name='mask0')
    ia.close()
    makemask(mode='copy', inpimage='temp.mask', inpmask='temp.mask:mask0', output=output)
    rmtables('temp.mask')

def count_channel(imagename, outDir, threshold, chans, axes=3):
    os.mkdir('temp')
    os.chdir('temp')
    name=imagename.split('/')[-1]
    outfile=re.sub('\..*$', '.mask', name)
    mask_threshold(imagename,threshold,100,outfile)
    countmask=outfile.replace('.mask', '.nchan')
    imcollapse(imagename=outfile,
           outfile=countmask,
           function='sum',
           axes=axes,   # 3 for ALMA, 2 for SMA
           chans=chans,
           overwrite=True)
    fitsimage=countmask.replace('.nchan','_nchan.fits')
    exportfits(imagename=countmask,fitsimage=fitsimage, overwrite=True)
    call(['cp', fitsimage, outDir])
    os.chdir('..')
    call(['rm','-rf','temp'])
    

############################################################
# main program

# go to the work directory
os.mkdir('temp')
os.chdir('temp')

# in the log/basic.txt, check the threshold used for makeing ratio map. 

# image= 'NGC5257_12CO10_combine_smooth.image'
# outfile='NGC5257_12CO10_combine_smooth.mask'
# mask_threshold(image,threshold_12,100,outfile)

# countmask= 'NGC5257_12CO10_combine_smooth_nchan.mask'
# imcollapse(imagename= outfile,
#            outfile=countmask,
#            function='sum',
#            axes=3,   # 3 for ALMA
#            chans='30~90',
#            overwrite=True)

# fitsimage=countmask.replace('.mask','.fits')
# exportfits(imagename=countmask,fitsimage=fitsimage)

image=imageDir+'13CO10/NGC5257_13CO10_12m_contsub_smooth.image/'
outfile='NGC5257_13CO10_12m_smooth.mask'
mask_threshold(image,threshold_13,100,outfile)

countmask='NGC5257_13CO10_12m_smooth_nchan.mask'
imcollapse(imagename=outfile,
           outfile=countmask,
           function='sum',
           axes=3,   # 3 for ALMA
           chans='5~30',
           overwrite=True)

fitsimage=countmask.replace('.mask','.fits')
exportfits(imagename=countmask,fitsimage=fitsimage)

call(['mv',fitsimage,imageDir+'13CO10/'])
os.chdir('..')
call(['rm','-rf','temp'])


### smooth the 
image=imageDir+'13CO10/NGC5257_13CO10_12m_contsub_smooth.image/'
outfile='NGC5257_13CO10_12m_smooth.mask'
mask_threshold(image,threshold_13,100,outfile)

countmask='NGC5257_13CO10_12m_smooth_nchan.mask'
imcollapse(imagename=outfile,
           outfile=countmask,
           function='sum',
           axes=3,   # 3 for ALMA
           chans='5~30',
           overwrite=True)

fitsimage=countmask.replace('.mask','.fits')
exportfits(imagename=countmask,fitsimage=fitsimage)

call(['mv',fitsimage,imageDir+'13CO10/'])
os.chdir('..')
call(['rm','-rf','temp'])


# count the channel number for smoothed 13CO 1-0 cube. 
imagename=imageDir+'13CO10/NGC5257_13CO10_12m_uvrange_smooth_co32.image'
outDir=imageDir+'13CO10/'
threshold=2*7.3e-4
chans='5~30'

count_channel(imagename, outDir, threshold, chans)


# count the channel number for 12CO 3-2 cube. 
imagename=imageDir+'12CO32/NGC5257co32_all_map40r.image'
outDir=imageDir+'12CO32/'
threshold=2*0.09
chans='10~20'

count_channel(imagename, outDir, threshold, chans, axes=2)
