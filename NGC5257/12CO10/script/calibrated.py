import numpy as np
import os
import scipy.ndimage as sni
import sys

data_12m = '/home/heh15/Data/Arp240/mnt/fhgfs/sliwa/ULIRGS/'\
           'Arp240/12CO10/12m/'\
           '2015.1.00804.S/'\
           'science_goal.uid___A001_X2d6_X25b/'\
           'group.uid___A001_X2d6_X25c/member.uid___A001_X2d6_X25d/'\
           'calibrated/arp240spw01_12m_co10.ms/'

data_7m = '/home/heh15/Data/Arp240/mnt/fhgfs/sliwa/ULIRGS/'\
          'Arp240/12CO10/7m/'\
          '2015.1.00804.S-12CO10/'\
          'science_goal.uid___A001_X2d6_X25b/group.uid___A001_X2d6_X25c/'\
          'member.uid___A001_X2d6_X25f/calibrated/arp240spw01-aca.ms/'

calibration_dir='/1/home/heh15/workingspace/Arp240/'\
                '12CO10/calibrated/'

combine_data='NGC5258_combine.ms'

origDir=os.getcwd()
os.chdir(calibration_dir)

split(vis=data_7m,outputvis='NGC5258_7m.ms',field='1',
      datacolumn='data',keepflags=False)
   
split(vis=data_12m,outputvis='NGC5258_12m.ms',field='1',
      datacolumn='data',keepflags=False)
 
concat(vis=['NGC5258_7m.ms','NGC5258_12m.ms'],
       concatvis=combine_data)

# split(vis=combine_data,spw='0',outputvis='NGC5258_combine_CO.ms')

fitspw='0:0~500,0:700~900,1,2:0~500;700~900,3'
spw='0,2'
uvcontsub(vis='NGC5258_combine.ms',fitorder=1,fitspw='0:0~500,0:700~900,1,2:0~500;700~900,3',spw='0,2')
os.rename('NGC5258_combine.ms.contsub','NGC5258_combine_CO.ms.contsub')

os.chdir(origDir)
