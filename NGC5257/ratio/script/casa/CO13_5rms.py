'''
Mar 13th, 2019
Create the file to give the 5rms cut for the 13CO10 image.
'''

from subprocess import call
import glob


Dir='/1/home/heh15/workingspace/Arp240/NGC5257/ratio/'
imageDir=Dir+'image/'
scriptDir=Dir+'script/'
workDir=Dir+'mask/'
regionDir=Dir+'region/'
maskDir=Dir+'mask/'

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

image=imageDir+'13CO10/NGC5257_13CO10_12m_contsub_smooth.mom0'
lower=2*0.00062*math.sqrt(25)*20
upper=100

makemask(mode='delete',inpmask=image+':mask0')
ia.open(image)
ia.calcmask('"'+image+'">'+str(lower)+' && '+'"'+image+'" <'+str(upper),name='mask0')
ia.close()

makemask(mode='delete',inpmask=image2+':mask0')
image2=imageDir+'13CO10/NGC5257_13CO10_12m_contsub_smooth_pbcor.mom0'
ia.open(image2)
ia.maskhandler('copy',[image+':mask0','mask10'])
ia.close()

# image=imageDir+'13CO10/NGC5257_13CO10_12m_contsub_smooth.mom0'
# lower=2*0.00062*math.sqrt(25)*20
# upper=100
# output=maskDir+'13CO_2rms_mom0.mask'
# mask_threshold(image,lower,upper,output)
