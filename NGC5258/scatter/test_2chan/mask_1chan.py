'''
Oct. 23rd, 2018

mask the pixel with only one channel selected

Dir=/home/heh15/workingspace/Arp240/scatter/test_2chan/
'''


############################################################
# functions

# with single beam
def mask_shape(image,region,output):
    makemask(inpimage=image,inpmask=region,mode='copy',output=output)

def mask_threshold(image,lower,upper,output):
    ia.open(image)
    ia.calcmask(image+'>'+str(lower)+' && '+image+'<'+str(upper),name='mask0')
    ia.close()
    makemask(inpimage=image,inpmask=image+':mask0',mode='copy',output=output)
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

# with multiple beams

def masks_threshold(image,lower,upper,output):
    ia.open(image)
    ia.calcmask(image+'>'+str(lower)+' && '+image+'<'+str(upper),name='mask0')
    ia.close()
    immath(imagename=[image,image],
           expr='IM0/IM1',
           outfile=output)
    makemask(mode='delete',inpmask=image+':mask0')

############################################################
# main program

image= 'NGC5257_12CO21_combine_rebin.image'
outfile='NGC5257_12CO21_combine_mask_rebin.mask'
rms=0.0031
threshold=4*rms
masks_threshold(image,threshold,100,outfile)

countmask='NGC5257_12CO21_combine_rebin_nchan.mask'
imcollapse(imagename= outfile,
           outfile=countmask,
           function='sum',
           axes=3,   # 3 for ALMA
           chans='10~60',
           overwrite=True)

# exportfits(imagename='NGC5257_12CO21_combine_mask_4sig.mask',fitsimage='NGC5257_12CO21_combine_mask_4sig.fits')
