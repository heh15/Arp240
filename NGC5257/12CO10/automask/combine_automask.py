x1############################################################
# Import the packages
import numpy as np
import os
import scipy.ndimage
import sys

############################################################
# Image parameter

vis='arp240_combine_CO.ms'
prename='NGC_5257_CO10_combine'
myimage=prename+'.image'
myflux=prename+'.flux'
mymask=prename+'.mask'
myresidual=prename+'.residual'
imsize=[320,320]
cell='0.5arcsec'
minpb=0.2
restfreq='112.73GHz'
outframe='bary'
spw='0~1'
width='10km/s'
start='-600km/s'
nchan=120
robust=0.5
phasecenter=0
scales=[0]
smallscalebias=0.6

stop=3.
pixelmin=0.5

############################################################
# Make inital dirty image
os.system('rm -rf '+prename+'.* ' +prename+'_*')
clean(vis=vis,imagename=prename,
      imagermode='mosaic',ftmachine='mosaic',minpb=minpb,
      imsize=imsize,cell=cell,spw=spw,
      weighting='briggs',robust=robust,phasecenter=phasecenter,
      mode='velocity',width=width,start=start,nchan=nchan,      
      restfreq=restfreq,outframe=outframe,veltype='radio',
      mask='',
      niter=0,interactive=False)

# Determine the beam area in pixels for later removal of very small mask regions
major=imhead(imagename=myimage,mode='get',hdkey='beammajor')['value']
minor=imhead(imagename=myimage,mode='get',hdkey='beamminor')['value']
pixelsize=float(cell.split('arcsec')[0])
beamarea=(major*minor*pi/(4*log(2)))/(pixelsize**2)
print 'beamarea in pixels =', beamarea

############################################################
# Find properties of dirty image.
myimage=prename+'.image'
bigstat=imstat(imagename=myimage)
peak= bigstat['max'][0]
print 'peak (Jy/beam) in cube = '+str(peak)
### Sets threshold of first loop, try 2-4. Subsequent loops are set thresh/2.
thresh = peak / 4.

### If True: find the rms in two line-free channels; If False:  Set rms by hand in else statement.
if True:  
    chanstat=imstat(imagename=myimage,chans='4')
    rms1= chanstat['rms'][0]
    chanstat=imstat(imagename=myimage,chans='66')
    rms2= chanstat['rms'][0]
    rms=0.5*(rms1+rms2)        
else:
    rms=0.011


print 'rms (Jy/beam) in a channel = '+str(rms)

############################################################
# Automasking Loop
os.system('rm -rf ' + prename +'_threshmask*')
os.system('rm -rf ' + prename +'_fullmask*')
os.system('rm -rf ' + prename +'.image*')
n=-1
while (thresh >= stop*rms):   
    n=n+1
    print 'clean threshold this loop is', thresh
    threshmask = prename+'_threshmask' +str(n)
    maskim = prename+'_fullmask' +str(n)
    immath(imagename = [myresidual],
           outfile = threshmask,
           expr = 'iif(IM0 > '+str(thresh) +',1.0,0.0)',
           mask=myflux+'>'+str(minpb))
    if (n==0):
        os.system('cp -r '+threshmask+' '+maskim+'.pb')
        print 'This is the first loop'
    else:
        makemask(mode='copy',inpimage=myimage,
                 inpmask=[threshmask,mymask],
                 output=maskim)
        imsubimage(imagename=maskim, mask=myflux+'>'+str(minpb),
                   outfile=maskim+'.pb')     
    print 'Combined mask ' +maskim+' generated.'

    # Remove small masks
    os.system('cp -r '+maskim+'.pb ' +maskim+'.pb.min')
    maskfile=maskim+'.pb.min'
    ia.open(maskfile)
    mask=ia.getchunk()           
    labeled,j=scipy.ndimage.label(mask)                     
    myhistogram = scipy.ndimage.measurements.histogram(labeled,0,j+1,j+1)
    object_slices = scipy.ndimage.find_objects(labeled)
    threshold=beamarea*pixelmin
    for i in range(j):
        if myhistogram[i+1]<threshold:
            mask[object_slices[i]] = 0


    ia.putchunk(mask)
    ia.done()
    print 'Small masks removed and ' +maskim +'.pb.min generated.'

    os.system('rm -rf '+mymask+'')
    clean(vis=vis,imagename=prename,
          imagermode='mosaic',ftmachine='mosaic',minpb=minpb,
          imsize=imsize,cell=cell,spw=spw,
          weighting='briggs',robust=robust,phasecenter=phasecenter,
          mode='velocity',width=width,start=start,nchan=nchan,      
          restfreq=restfreq,outframe=outframe,veltype='radio',
          mask = maskim+'.pb.min',
          multiscale=scales,smallscalebias=smallscalebias,
          interactive = False,
          niter = 10000,
          threshold = str(thresh) +'Jy/beam')


    if thresh==stop*rms: break
    thresh = thresh/2.
    # Run a final time with stop*rms if more than a little above
    # stop*rms. Also make a back-up of next to last image
    if (thresh < stop*rms and thresh*2.>1.05*stop*rms):
        thresh=stop*rms  
        os.system('cp -r '+myimage+' '+myimage+str(n))
##############################################
# Apply a primary beam correction

import glob

myimages = glob.glob("*.image")

rmtables('*.pbcor')
for image in myimages:
    pbimage = image.rsplit('.',1)[0]+'.flux'
    outfile = image.rsplit('.',1)[0]+'.pbcor'
    impbcor(imagename=image, pbimage=pbimage, outfile = outfile)

##############################################
# Create Diagnostic PNGs


os.system("rm -rf *.png")
mycontimages = glob.glob("*.image")
for cimage in mycontimages:
    max=imstat(cimage)['max'][0]
    min=-0.1*max
    outimage = cimage+'.png'
    os.system('rm -rf '+outimage)
    imview(raster={'file':cimage},out=outimage)


# this will have to be run for each sourcename
sourcename='NGC_5257' # insert source here, if it isn't already set
mylineimages = glob.glob(sourcename+"*.image")
for limage in mylineimages:
    rms=imstat(limage,chans='1')['rms'][0]
    mom0=limage+'.mom0'
    os.system("rm -rf "+mom0)
    immoments(limage,moments=[0],outfile=mom0)
    os.system("rm -rf "+mom0+".png")
    imview(raster={'file':mom0},out=mom0+'.png')
