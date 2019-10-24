'''
make the mask based on the 33GHz SFR continuum map. 

1) The shape is approximately 5rms boundary. 
'''
import os
from subprocess import call

############################################################
# directory
Dir='/1/home/heh15/workingspace/Arp240/NGC5258/ratio/'
regionDir=Dir+'region/'
imageDir=Dir+'image/'
maskDir=Dir+'mask/'
crtDir=os.getcwd()

############################################################
# function
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

def mask_or(mask1,mask2,output):
    immath(imagename=[mask1,mask2],expr='IM0+IM1',outfile='temp.mask')
    ia.open('temp.mask')
    ia.calcmask('"temp.mask">0.5',name='mask0')
    ia.close()
    makemask(mode='copy', inpimage='temp.mask', inpmask='temp.mask:mask0', output=output)
    rmtables('temp.mask')

############################################################

image=imageDir+'33GHz/NGC5258_33GHz_smooth_regrid.image/'
region=regionDir+'south_SFR.crtf'
output=maskDir+'south_SFR.mask'
mask_shape(image,region,output)

exportfits(imagename=output,fitsimage=output.replace('.mask','_mask.fits'))

image=imageDir+'33GHz/NGC5258_33GHz_smooth_regrid.image/'
region=regionDir+'whole.crtf'
output=maskDir+'whole.mask'
mask_shape(image,region,output)

exportfits(imagename=output,fitsimage=output.replace('.mask','_mask.fits'))

mask1=maskDir+'whole.mask'
mask2=maskDir+'south_SFR.mask'
output=maskDir+'rest_SFR.mask'
mask_subtract(mask1,mask2,output)

exportfits(imagename=output,fitsimage=output.replace('.mask','_mask.fits'))

# # merge the 3 masks to the one. 
# mask1=maskDir+'hinge_SFR.mask'
# mask2=maskDir+'center_SFR.mask'
# output='temp1.mask'
# mask_or(mask1,mask2,output)

# mask1=maskDir+'south_SFR.mask'
# mask2='temp1.mask'
# output='temp2.mask'
# mask_or(mask1,mask2,output)

# # subtract the merged mask from the whole. 
# region=regionDir+'whole.crtf'
# output=maskDir+'whole_init.mask'
# mask_shape(image,region,output)

# mask1=maskDir+'whole_init.mask'
# mask2='temp2.mask'
# output=maskDir+'rest.mask'
# mask_subtract(mask1,mask2,output)

os.chdir('..')
call(['rm','-rf','temp'])
