'''
Mar. 26th, 2018
casa 5.1.1

'''
import numpy as np
import os,glob
import scipy.ndimage as sni
import sys
import re

imageDir='../12CO10/continuum/'
imagenames=glob.glob(imageDir+'*spw*.image')
region_peak='/home/heh15/workingspace/Arp240/'\
        'continuum/south_source.crtf'
region_rms='/home/heh15/workingspace/Arp240/'\
            'continuum/south_source.crtf'
uncertainty_12CO={}


for imagename in imagenames:
    # dictionary name
    imagename_tmp=re.split('\_|\.|/',imagename)
    CO=[string for string in imagename_tmp\
        if re.match('.*CO.*',string)][0]
    spw=[string for string in imagename_tmp\
         if re.match('spw.*',string)][0]
    keyname=CO+'_'+spw
    keyname_tmp={}
    # measure flux and rms
    keyname_tmp['peak']=imstat(imagename=imagename,region=region_peak)\
                                  ['max']['0']
    keyname_tmp['rms']=imstat(imagename=imagename,region=region_rms)\
                                 ['rms']['0']
    keyname_tmp['uncertainty']=keyname_tmp['rms']/keyname_tmp['peak']
    uncertainty_12CO[keyname]=keyname_tmp
                                
    
