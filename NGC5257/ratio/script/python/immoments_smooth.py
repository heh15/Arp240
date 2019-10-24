from subprocess import call

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/ratio/'
imageDir=Dir+'image/'
scriptDir=Dir+'script/'
workDir=Dir+'anomaly_kin/'
regionDir=Dir+'region/'

os.chdir(imageDir)

imagename=imageDir+'NGC5257_12CO10_combine_contsub_smooth.image/'
outfile='NGC5257_12CO10_combine_contsub_smooth.mom0'
immoments(imagename=imagename,moments=0,outfile=outfile,chans='10~60')

outfile='NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid.mom0'
immoments(imagename=imagename,moments=0,outfile=outfile,chans='10~60')


immoments(imagename=imagename,moments=0,outfile=outfile,chans='5~30')



