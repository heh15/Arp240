'''
Jan. 18th, 2019
'''

############################################################
# directory

Dir='/home/heh15/workingspace/Arp240/NGC5257/12CO21/'
imageDir=Dir+'casa5.4/'

default("simobserve")
project="NGC5257"

skymodel ='NGC5257_12CO21_combine_contsub_uvtaper.image/'

indirection='J2000 13h39m52.922s 0d50m24.1s'
incell="0.3arcsec"
inbright="0.00524"
incenter='225.46GHz'
inwidth="30MHz"

antennalist='sma.subcompact.cfg'
totaltime="12800s"

mapsize="1arcmin"
dryrun = False
simobserve()
