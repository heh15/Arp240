'''
Aug 29th, 2018

'''
from subprocess import call
import glob


Dir='/1/home/heh15/workingspace/Arp240/NGC5257/ratio/'
imageDir=Dir+'image/'
scriptDir=Dir+'script/'
workDir=Dir+'mask/'
regionDir=Dir+'region/'

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


############################################################
# main program

os.chdir(workDir)

immoments(imagename=imageDir+'NGC5257_12CO10_combine_contsub_smooth.image',chans='10~60',outfile='NGC5257_12CO10_combine_contsub_smooth.mom0')
immoments(imagename=imageDir+'NGC5257_13CO10_12m_contsub_smooth.image',chans='5~30',outfile='NGC5257_13CO10_combine_contsub_smooth.mom0')
immoments(imagename=imageDir+'NGC5257_12CO21_combine_contsub_uvtaper_smooth.image',chans='10~60',outfile='NGC5257_12CO21_combine_contsub_uvtaper_smooth.mom0')

Im12CO10='NGC5257_12CO10_combine_contsub_smooth.mom0'
Im13CO10='NGC5257_13CO10_combine_contsub_smooth.mom0'

# copy the image into the current directory
# subprocess.call(['cp','-r',imageDir+Im12CO10,'.'])

ratioDir=Dir+'1213/contsub/'
Im12CO10=ratioDir+'NGC5257_12CO10_combine_contsub_smooth_masked.image.mom0'

image=Im12CO10
lower=2.0; upper=5.0
output='spiralarm_temp.mask'
mask_threshold(image,lower,upper,output)

image=Im12CO10
lower=1.0; upper=100.0
output='lowercut.mask'
mask_threshold(image,lower,upper,output)
image=Im12CO10
region=regionDir+'anomaly.crtf'
output='anomaly_tmp.mask'
mask_shape(image,region,output)
mask1='lowercut.mask'
mask2='anomaly_tmp.mask'
output='anomaly.mask'
mask_and(mask1,mask2,output)
rmtables('lowercut.mask')

image=Im12CO10
region=regionDir+'nonarm.crtf'
output='nonarm.mask'
mask_shape(image,region,output)


mask1='spiralarm_temp.mask'
mask2='anomaly.mask'
output='spiralarm_temp1.mask'
mask_subtract(mask1,mask2,output)

mask1='spiralarm_temp1.mask'
mask2='nonarm.mask'
output='spiralarm.mask'
mask_subtract(mask1,mask2,output)

# disk without spiral arm
image=Im12CO10
lower=0.5; upper=2.0
output='disk.mask'
mask_threshold(image,lower,upper,output)

# central region
image=Im12CO10
lower=5.0;upper=100
output='center.mask'
mask_threshold(image,lower,upper,output)

# anomalous region shown in the moment 1 map. 
mom1_12CO10=imageDir+'NGC5257_12CO10_combine_contsub.mom1'
region=regionDir+'anomaly_kin.crtf'
output='anomaly_kin_int.mask'
mask_shape(mom1_12CO10,region,output)

lower=0;upper=300
output='away.mask'
mask_threshold(mom1_12CO10,lower,upper,output)

mask1='anomaly_kin_int.mask'
mask2='away.mask'
output='anomaly_kin.mask'
mask_and(mask1,mask2,output)


# convert them to 12CO21 image. 

# imregrid(imagename='anomaly.mask',template='NGC5257_12CO21_combine_uvtaper_smooth_masked.image.mom0',output='anomaly_12CO21.mask')

# imregrid(imagename='nonarm.mask',template='NGC5257_12CO21_combine_uvtaper_smooth_masked.image.mom0',output='nonarm_12CO21.mask')

# imregrid(imagename='spiralarm.mask',template='NGC5257_12CO21_combine_uvtaper_smooth_masked.image.mom0',output='spiralarm_12CO21.mask')

# imregrid(imagename='restdisk.mask',template='NGC5257_12CO21_combine_uvtaper_smooth_masked.image.mom0',output='restdisk_12CO21.mask')

# imregrid(imagename='center.mask',template='NGC5257_12CO21_combine_uvtaper_smooth_masked.image.mom0',output='center_12CO21.mask')

# move the file to the directory
lists=glob.glob('*.mask')
for mask in lists:
    exportfits(imagename=mask,fitsimage=mask.replace('.mask','_mask.fits'))

