'''
Mar. 22nd, 2018
casa 5.1.1

'''
import glob, os
import re
from distutils.dir_util import copy_tree

data=glob.glob('../*CO*')
imageDirs=[string+'/NGC5258/finalImage/' for string in data]
for imageDir in imageDirs:
    imagenames=glob.glob(imageDir+'*.image')
    imagename=[string for string in imagenames\
           if re.match('.*center.*',string)==None][0]
    image_new=imagename.split('/')[-1]
    copy_tree(imagename,'./NGC5258/startImage/'+image_new)


imhead(imagename=imagename,mode='list',hdkey='beammajor')['perplanebeams']['median area beam']


imsmooth(imagename=imagename,
         kernel='gauss',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.314deg',
         targetres=True,
         outfile='NGC5258_13CO10_12m_smooth.image')


          
