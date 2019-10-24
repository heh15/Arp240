'''
Mar. 22nd, 2018
casa 5.1.1

Nov. 19th, 2018
casa 5.4.0

'''
import glob, os
import re
from distutils.dir_util import copy_tree

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/ratio/'
imageDir=Dir+'image/'
scriptDir=Dir+'script/'
regionDir=Dir+'region/'

############################################################
# main program
os.chdir(imageDir)

# data=glob.glob('../*CO*')
# imageDirs=[string+'/NGC5258/finalImage/' for string in data]
# for imageDir in imageDirs:
#     imagenames=glob.glob(imageDir+'*.image')
#     imagename=[string for string in imagenames\
#            if re.match('.*center.*',string)==None][0]
#     image_new=imagename.split('/')[-1]
#     copy_tree(imagename,'./NGC5258/startImage/'+image_new)

rmtables('*_smooth.image')
images=glob.glob('*CO*.image')

# imhead(imagename=imagename,mode='list',hdkey='beammajor')['perplanebeams']['median area beam']

for imagename in images:
    imsmooth(imagename=imagename,
             kernel='gauss',
             major='2.186arcsec',
             minor='1.896arcsec',
             pa='-87.314deg',
             targetres=True,
             outfile=imagename.replace('.image','_smooth.image'))

imagename=glob.glob('*12CO21*uvtaper_smooth.image')[0]
template=glob.glob('*12CO10*contsub_smooth.image')[0]
imregrid(imagename=imagename,template=template,output=imagename.replace('.image','_regrid.image'))


          
