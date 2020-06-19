from subprocess import call
import glob,os
import numpy as np
import math

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/ratio/'
imageDir=Dir+'image/'
scriptDir=Dir+'script/'
workDir=Dir+'anomaly_kin/'
regionDir=Dir+'region/'

beamarea=2.186*1.896*1.1331
beamarea_pix=beamarea/0.09

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

# 12CO10
makemask(inpmask=regionDir+'anomaly_kin.mask',inpimage=imageDir+'NGC5257_12CO10_combine_contsub_smooth.image/',output='anomaly_kin_chans.mask',mode='expand')
immoments(imagename=imageDir+'NGC5257_12CO10_combine_contsub_smooth.image',moments=0,mask='anomaly_kin_chans.mask',outfile='NGC5257_12CO10_anomaly_kin.mom0',chans='36~46')
makemask(mode='delete',inpmask='NGC5257_12CO10_anomaly_kin.mom0:mask0')
imagename=imageDir+'NGC5257_12CO10_combine_contsub_smooth.image/'
region=regionDir+'anomaly_kin.crtf'
flux10=imstat(imagename=imagename,region=region,chans='36~46')['flux'][0]
imagename=imageDir+'NGC5257_12CO10_combine_contsub_smooth.image/'
# outfile='NGC5257_12CO10_combine_contsub_smooth.mom0'
#immoments(imagename=imagename,moments=0,outfile=outfile,chans='10~60')
#call(['mv',outfile,imageDir])
rms=0.0016;Npts=387;chans=10
flux10_error=rms*10*np.sqrt(chans*Npts/beamarea_pix)

# subtract the anomaly_kin from the moment 0 map
image='NGC5257_12CO10_combine_contsub_smooth.mom0'
IM0=imageDir+'NGC5257_12CO10_combine_contsub_smooth.mom0'
IM1='NGC5257_12CO10_anomaly_kin.mom0'
#outfile=IM0.replace('.mom0','_subtr.mom0')
outfile=image.replace('.mom0','_subtr.mom0')
immath(imagename=[IM0,IM1],
       expr='IM0-IM1',
       outfile=outfile)

imagename=outfile
region=regionDir+'south_west_subtr.crtf'
beamarea_pix=1.1331*2.186*1.896/(0.09)
SW_10=imstat(imagename=imagename,region=region)['sum'][0]/beamarea_pix
Npts=353;chans=21
SW10_error=rms*10*np.sqrt(chans*Npts/beamarea_pix)


# 12CO 2-1 
inpmask=regionDir+'anomaly_kin.mask'
inpimage=imageDir+'NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid.image/'
output='anomaly_kin_chans.mask'
makemask(inpmask=inpmask,inpimage=inpimage,output=output,mode='expand')

imagename=inpimage
mask=output+'==1'
outfile='NGC5257_12CO21_anomaly_kin.mom0'
immoments(imagename=imagename,moments=0,mask=mask,outfile=outfile,chans='36~46')
makemask(mode='delete',inpmask='NGC5257_12CO21_anomaly_kin.mom0:mask0')
region=regionDir+'anomaly_kin.crtf'
flux21=imstat(imagename=imagename,region=region,chans='36~46')['flux'][0]
outfile='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid.mom0'
#immoments(imagename=imagename,moments=0,outfile=outfile,chans='10~60')
#call(['mv',outfile,imageDir])
rms=0.003;Npts=387;chans=10
flux21_error=rms*10*np.sqrt(chans*Npts/beamarea_pix)


image='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid.mom0'
IM0=imageDir+'NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid.mom0'
IM1='NGC5257_12CO21_anomaly_kin.mom0'
#outfile=IM0.replace('.mom0','_subtr.mom0')
outfile=image.replace('.mom0','_subtr.mom0')
immath(imagename=[IM0,IM1],
       expr='IM0-IM1',
       outfile=outfile)
imagename=outfile
region=regionDir+'south_west_subtr.crtf'
SW_21=imstat(imagename=imagename,region=region)['sum'][0]/beamarea_pix
Npts=353;chans=21
SW21_error=rms*10*np.sqrt(chans*Npts/beamarea_pix)

# 13CO 1-0 no signal
imagename=imageDir+'NGC5257_13CO10_12m_contsub_smooth.image'
region=regionDir+'anomaly_kin.crtf'
outfile='NGC5257_13CO10_12m_contsub_smooth.mom0'
upper=imstat(imagename=imagename,chans='18~23',region=region)['flux'][0]
#immoments(imagename=imagename,moments=0,outfile=outfile,chans='5~30')
#call(['mv',outfile,imageDir])

############################################################
# move the result files to certain directories. 

# # move the file to directory
# lists=glob.glob('*subtr.mom0')
# for image in lists:
#     call(['mv',image,imageDir])

ratio1=flux21/flux10*115.27**2/230.54**2
uncertainty1=math.sqrt((flux21_error/flux21)**2+(flux10_error/flux10)**2)*ratio1

ratio2=SW_21/SW_10*115.27**2/230.54**2
uncertainty2=math.sqrt((SW21_error/SW_21)**2+(SW10_error/SW_10)**2)*ratio2

