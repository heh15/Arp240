'''
Mar. 21st, 2018
'''
import numpy as np
import os,glob
import scipy.ndimage as sni
import sys
import re


imageDir='/home/heh15/practice/Arp240/'\
          'workingspace/CO13_10/continuum/'
region_init='/home/heh15/practice/Arp240/'\
        'workingspace/continuum/south_source.crtf'

origDir = os.getcwd()
os.chdir(imageDir)
imagenames=glob.glob("NGC5257_cont_12m_13CO_spw*.image")
flux_13CO={}

for imagename in imagenames:
    # dictionary name
    imagename_tmp=re.split('\_|\.',imagename)
    CO=[string for string in imagename_tmp\
        if re.match('.*CO.*',string)][0]
    spw=[string for string in imagename_tmp\
         if re.match('spw.*',string)][0]
    keyname=CO+'_'+spw
    # flux in beam
    beam=imfit(imagename=imagename,region=region_init)
    x_value=beam['results']['component0']['shape']\
             ['direction']['m0']['value']
    y_value=beam['results']['component0']['shape']\
             ['direction']['m1']['value']
    bmaj_value=beam['results']['component0']\
                ['shape']['majoraxis']['value']
    bmin_value=beam['results']['component0']['shape']\
                ['minoraxis']['value']
    pa_value=beam['results']['component0']['shape']\
              ['positionangle']['value']
    x=str(x_value)+'rad'
    y=str(y_value)+'rad'
    bmaj=str(bmaj_value)+'arcsec'
    bmin=str(bmin_value)+'arcsec'
    pa=str(pa_value)+'deg'
    region='ellipse[['+x+','+y+'],['+bmaj+','+bmin+'],'+pa+']'
    flux_13CO[keyname]=imstat(imagename=imagename,region=region)\
                        ['flux']['0']

os.chdir(origDir)


    


