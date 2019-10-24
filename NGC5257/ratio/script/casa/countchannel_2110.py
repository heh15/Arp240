

############################################################
# basic settings 

Dir= '/1/home/heh15/workingspace/Arp240/ratio/'
workDir= Dir+'NGC5257/2110/uvtaper/'
scriptDir=Dir+'script/'
threshold_10= 5*0.001497
threshold_21= 5*3.42*0.001497

############################################################
# functions

def mask_shape(image,region,output):
    makemask(inpimage=image,inpmask=region,mode='copy',output=output)

def mask_threshold(image,lower,upper,output):
    ia.open(image)
    ia.calcmask(image+'>'+str(lower)+' && '+image+'<'+str(upper),name='mask0')
    ia.close()
    makemask(inpimage=image,inpmask=image+':mask0',mode='copy',output=output,overwrite=True)
    makemask(mode='delete',inpmask=image+':mask0')

def mask_invert(mask,output):
    ia.open(image)
    ia.calcmask(image+'< 0.5')
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

############################################################
# main program

# go to the work directory
os.chdir(workDir)

# in the log/basic.txt, check the threshold used for makeing ratio map. 

image= 'NGC5257_12CO10_combine_uvrange_smooth_regrid21.image'
outfile=image.replace('.image', '.mask')
mask_threshold(image,threshold_10,100,outfile)

countmask= outfile.replace('.mask','_nchan.mask')
imcollapse(imagename= outfile,
           outfile=countmask,
           function='sum',
           axes=3,   # 3 for ALMA
           chans='10~60',
           overwrite=True)

fitsimage=countmask.replace('.mask','.fits')
exportfits(imagename=countmask,fitsimage=fitsimage,overwrite=True)

image='NGC5257_12CO21_combine_uvtaper_smooth.image'
outfile=image.replace('.image', '.mask')
mask_threshold(image,threshold_21,100,outfile)

countmask=outfile.replace('.mask','_nchan.mask')
imcollapse(imagename= outfile,
           outfile=countmask,
           function='sum',
           axes=3,   # 3 for ALMA
           chans='10~60',
           overwrite=True)

fitsimage=countmask.replace('.mask','.fits')
exportfits(imagename=countmask,fitsimage=fitsimage,overwrite=True)

os.chdir(scriptDir)
