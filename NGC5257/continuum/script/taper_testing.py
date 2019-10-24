'''
This is very much still a work in progress so use at your own risk. brunettn
can help with questions but a more finished version should be coming soon
(haha what does soon mean?).

TODO
    -consider normalizing position angle differences when calculating closeness
     to target metric
    -figure out how to pass parameters into this code
      -along with this will be deciding how I want to run it (most convenient
       is probably from within CASA)
      -info script needs appears to be:
        *MS
        *storageDir
        *field
        *spwStr
        *spwList
        *antenna
        *observation
        *uvrange
        *nx
        *ny
        *cell
        *phasecenter
        *mode
        *robust
        *sourceBm
        *circularize
        *targetBm (if circularize=False)
        *startMaj
        *endMaj
        *nMaj
        *startMin
        *endMin
        *nMin
        *startPA
        *endPA
        *nPA
    -is the circularize targetBm always going to be useable (will the major
     axis have to increaes a bit in the final convolution, I thought so but
     first round of testing did not need to)
    -include option to just do it the old way (just an approximate grid) so
     there is more flexibility in how it is matched (e.g. very different PAs
     can mean having to make the beam larger than originally planned)
    -prevent ia.beamforconvolvedsize junk getting printed to the terminal


Hao: 
'''

import glob
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as so
import time

# quadratic function
def quad(x, a, b, c):
    return a*x**2 + b*x + c


#------------------------------------------------------------------------------#
# initialize some stuff
MS = '/home/heh15/Data/Arp240/arp240-110GHz.ms/'
scriptDir='/home/heh15/workingspace/Arp240/NGC5257/continuum/script/'
storageDir='/home/heh15/workingspace/Arp240/NGC5257/continuum/test_uvtaper/'

delmod(MS)

field = '0'
spwStr = '2,3'
spwList = np.arange(4)
antenna = ''
observation = ''
uvrange = ''
nx = 320
ny = 320
cell = '0.3arcsec'
phasecenter = 0
mode = 'mfs'
robust = 0.5

sourceBm = ['2.336arcsec', '1.765arcsec', '87.91deg']
targetBm = ['6.0arcsec', '6.0arcsec', '0deg']

circularize = False #this will be an argument to this code
if circularize:
    diam = sqrt(float(sourceBm[0].split('arcsec')[0])**2
                + 2.0*(float(cell.split('arcsec')[0])**2))
    targetBm = ['{:f}arcsec'.format(diam),
                '{:f}arcsec'.format(diam),
                sourceBm[2]]
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Get tapering first guess and define grid around it.
try:
    convKern = ia.beamforconvolvedsize(source=sourceBm, convolved=targetBm)
except RuntimeError:
    print '\nCannot reach target beam from the source beam.\n' \
          + 'Try experimenting with ia.beamforconvolvedsize first.'
    raise ValueError('See above.')
convKern = [convKern['major']['value'],
            convKern['minor']['value'],
            convKern['pa']['value']]

startMaj = convKern[0] - convKern[0]/4
endMaj = convKern[0] + convKern[0]/4
nMaj = 6
startMin = convKern[1] - convKern[1]/4
endMin = convKern[1] + convKern[1]/4
nMin = 6
startPA = convKern[2] - convKern[2]/4
endPA = convKern[2] + convKern[2]/4
nPA = 6
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Run the fast but approximate tapering grid.
# start the stopwatch
start = time.time()

# remember what files were in current directory before doing these tests
globBefore = set(glob.glob('./*'))

# initialize image within MS
delmod(vis=MS, field=field)
im.open(MS)
im.setvp(usedefaultvp=True)
im.selectvis(spw=spwStr, field=field, baseline=antenna,
             observation=observation, uvrange=uvrange)
im.defineimage(nx=nx, ny=ny, cellx=cell, celly=cell, phasecenter=phasecenter,
               mode=mode, spw=spwList)

# loop over all input taper parameters and keep output beam parameters
majAbscissa, minAbscissa, paAbscissa = \
                           np.meshgrid(np.linspace(startMaj, endMaj, num=nMaj),
                                       np.linspace(startMin, endMin, num=nMin),
                                       np.linspace(startPA, endPA, num=nPA))
resultMaj = list()
resultMin = list()
resultPA = list()
for i in np.arange(majAbscissa.shape[0]):
    for j in np.arange(majAbscissa.shape[1]):
        for k in np.arange(majAbscissa.shape[2]):
            im.weight('briggs', robust=robust)
            im.filter(type='gaussian',
                      bmaj=str(majAbscissa[i, j, k])+'arcsec',
                      bmin=str(minAbscissa[i, j, k])+'arcsec',
                      bpa=str(paAbscissa[i, j, k])+'deg')
            params = im.fitpsf(psf='')
            resultMaj.append(params[1]['value'])
            resultMin.append(params[2]['value'])
            resultPA.append(params[3]['value'])
im.done()
globAfter = set(glob.glob('./*'))
rmtables(list(globAfter - globBefore))
resultMaj = np.array(resultMaj)
resultMin = np.array(resultMin)
resultPA = np.array(resultPA)

# give us some execution stats
duration = np.round(time.time() - start, decimals=2)
print 'Taper test duration (approximate method):', duration, 'seconds'
print 'Individual test duration:', \
      np.round(duration/majAbscissa.size, decimals=2), 'seconds'
print 'Total individual tests run:', majAbscissa.size

# save inputs and outputs
writeMaj = majAbscissa.flatten()
writeMin = minAbscissa.flatten()
writePA = paAbscissa.flatten()
np.savetxt(storageDir+'approx_taper_output_'
           +time.strftime('%Y:%m:%d:%H:%M:%S', time.gmtime()),
           np.transpose([writeMaj, writeMin, writePA,
                         resultMaj, resultMin, resultPA]),
           delimiter='   ',
           header='in Maj   in Min   in PA   out Maj   out Min   out PA')
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Run the tapering grid extrema through tclean to get the actual output beams.
# start the stopwatch
start = time.time()

# loop over all input taper parameters and keep output beam parameters
majAbscissa, minAbscissa, paAbscissa = \
                              np.meshgrid(np.linspace(startMaj, endMaj, num=2),
                                          np.linspace(startMin, endMin, num=2),
                                          np.linspace(startPA, endPA, num=2))
resultMaj = list()
resultMin = list()
resultPA = list()
for i in np.arange(majAbscissa.shape[0]):
    for j in np.arange(majAbscissa.shape[1]):
        for k in np.arange(majAbscissa.shape[2]):
            delmod(vis=MS, field=field)
            tclean(vis=MS,
                   imagename='from_tclean',
                   field=field,
                   spw=spwStr,
                   antenna=antenna,
                   observation=observation,
                   uvrange=uvrange,
                   imsize=[nx, ny],
                   cell=cell,
                   phasecenter=phasecenter,
                   specmode='mfs',
                   gridder='mosaic',
                   weighting='briggs',
                   robust=robust,
                   uvtaper=['{}arcsec'.format(majAbscissa[i, j, k]),
                            '{}arcsec'.format(minAbscissa[i, j, k]),
                            '{}deg'.format(paAbscissa[i, j, k])],
                   niter=0,
                   calcpsf=True,
                   calcres=False,
                   restoration=False,
                   interactive=False)
            psf = imhead(imagename='from_tclean.psf')['restoringbeam']
            resultMaj.append(psf['major']['value'])
            resultMin.append(psf['minor']['value'])
            resultPA.append(psf['positionangle']['value'])
            rmtables('from_tclean.*')
delmod(vis=MS, field=field)
resultMaj = np.array(resultMaj)
resultMin = np.array(resultMin)
resultPA = np.array(resultPA)

# give us some execution stats
duration = np.round(time.time() - start, decimals=2)
print 'tclean taper test duration:', duration, 'seconds'
print 'Individual test duration:', \
      np.round(duration/majAbscissa.size, decimals=2), 'seconds'
print 'Total individual tests run:', majAbscissa.size

# save inputs and outputs
writeMaj = majAbscissa.flatten()
writeMin = minAbscissa.flatten()
writePA = paAbscissa.flatten()
np.savetxt(storageDir+'tclean_taper_output_'+ \
           time.strftime('%Y:%m:%d:%H:%M:%S', time.gmtime()), \
           np.transpose([writeMaj, writeMin, writePA, \
                         resultMaj, resultMin, resultPA]), \
           delimiter='   ', \
           header='in Maj   in Min   in PA   out Maj   out Min   out PA')
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Determine best tapering input given the tests.
# load all test results in the storage directory
outputFiles = glob.glob(storageDir+'approx_taper_output_*')
approxInMaj, approxInMin, approxInPA, approxMaj, approxMin, approxPA = \
                                        np.loadtxt(outputFiles[0], unpack=True)
for i in np.arange(1, len(outputFiles)):
    tmp = np.loadtxt(outputFiles[i])
    approxInMaj = np.append(approxInMaj, tmp[:, 0])
    approxInMin = np.append(approxInMin, tmp[:, 1])
    approxInPA = np.append(approxInPA, tmp[:, 2])
    approxMaj = np.append(approxMaj, tmp[:, 3])
    approxMin = np.append(approxMin, tmp[:, 4])
    approxPA = np.append(approxPA, tmp[:, 5])
outputFiles = glob.glob(storageDir+'tclean_taper_output_*')
tcleanInMaj, tcleanInMin, tcleanInPA, tcleanMaj, tcleanMin, tcleanPA = \
                                        np.loadtxt(outputFiles[0], unpack=True)
for i in np.arange(1, len(outputFiles)):
    tmp = np.loadtxt(outputFiles[i])
    tcleanInMaj = np.append(tcleanInMaj, tmp[:, 0])
    tcleanInMin = np.append(tcleanInMin, tmp[:, 1])
    tcleanInPA = np.append(tcleanInPA, tmp[:, 2])
    tcleanMaj = np.append(tcleanMaj, tmp[:, 3])
    tcleanMin = np.append(tcleanMin, tmp[:, 4])
    tcleanPA = np.append(tcleanPA, tmp[:, 5])

# match same inputs to approximate method and tclean
tmp1 = np.unique(tcleanInMaj)
tmp2 = np.unique(approxInMaj)
tmp3 = np.where(tmp1 == tmp2)

matchInds = list()
for i in np.arange(tcleanInMaj.shape[0]):
    matchInds.append(np.where((approxInMaj == tcleanInMaj[i])
                              & (approxInMin == tcleanInMin[i])
                              & (approxInPA == tcleanInPA[i])))
matchInds = np.array(matchInds).flatten()

# fit relation between approximate beam dimensions and those from tclean
MajPOpt, pCov = so.curve_fit(quad, approxMaj[matchInds], tcleanMaj)
MinPOpt, pCov = so.curve_fit(quad, approxMin[matchInds], tcleanMin)
PAPOpt, pCov = so.curve_fit(quad, approxPA[matchInds], tcleanPA)

# plot relations between approximate and tclean beam dimensions
if False:
    x = np.linspace(-180, 180, num=100000)
    fig, ((ax1, ax2, ax3),
          (ax4, ax5, ax6),
          (ax7, ax8, ax9)) = plt.subplots(3, 3)

    ax1.plot(approxMaj[matchInds], tcleanMaj, 'ko', zorder=2)
    ax1.plot([0, 2], [0, 2], 'k', zorder=1)
    ax1.plot(x, quad(x, *MajPOpt), 'r--', zorder=3)
    ax2.plot(approxMin[matchInds], tcleanMin, 'ko', zorder=2)
    ax2.plot([0, 2], [0, 2], 'k', zorder=1)
    ax2.plot(x, quad(x, *MinPOpt), 'r--', zorder=3)
    ax3.plot(approxPA[matchInds], tcleanPA, 'ko', zorder=2)
    ax3.plot([-100, 80], [-100, 80], 'k', zorder=1)
    ax3.plot(x, quad(x, *PAPOpt), 'r--', zorder=3)
    ax1.set_xlim([1.7, 2])
    ax1.set_ylim([1.7, 2])
    ax2.set_xlim([1.35, 1.75])
    ax2.set_ylim([1.35, 1.75])
    ax3.set_xlim([-100, 80])
    ax3.set_ylim([-100, 80])

    ax4.plot(approxMaj[matchInds], tcleanMaj-approxMaj[matchInds],
             'ko')
    ax4.plot(x, quad(x, *MajPOpt)-x, 'r--')
    ax5.plot(approxMin[matchInds], tcleanMin-approxMin[matchInds], 'ko')
    ax5.plot(x, quad(x, *MinPOpt)-x, 'r--')
    ax6.plot(approxPA[matchInds], tcleanPA-approxPA[matchInds], 'ko')
    ax6.plot(x, quad(x, *PAPOpt)-x, 'r--')
    ax4.set_xlim([1.7, 2])
    ax4.set_ylim([-0.008, 0.002])
    ax5.set_xlim([1.35, 1.75])
    ax5.set_ylim([-0.008, 0])
    ax6.set_xlim([-100, 80])
    ax6.set_ylim([-2.5, 2])

    ax7.plot(approxMaj[matchInds],
             tcleanMaj-quad(approxMaj[matchInds], *MajPOpt), 'ko')
    ax7.axhline(0, color='k')
    ax8.plot(approxMin[matchInds],
             tcleanMin-quad(approxMin[matchInds], *MinPOpt), 'ko')
    ax8.axhline(0, color='k')
    ax9.plot(approxPA[matchInds],
             tcleanPA-quad(approxPA[matchInds], *PAPOpt), 'ko')
    ax9.axhline(0, color='k')
    ax7.set_xlim([1.7, 2])
    ax7.set_ylim([-0.006, 0.003])
    ax8.set_xlim([1.35, 1.75])
    ax8.set_ylim([-0.004, 0.005])
    ax9.set_xlim([-100, 80])
    ax9.set_ylim([-1, 1.5])

    plt.show()

# get "corrected" approximate beam dimensions using the fits
corrMaj = quad(approxMaj, *MajPOpt)
corrMin = quad(approxMin, *MinPOpt)
corrPA = quad(approxPA, *PAPOpt)

# calculate differences from target beam shape
metric = np.sqrt((corrMaj - float(targetBm[0].split('arc')[0]))**2 + \
                 (corrMin - float(targetBm[1].split('arc')[0]))**2 + \
                 (corrPA - float(targetBm[2].split('deg')[0]))**2)

# find closest match where axes are less than target and imsmooth will be
# successful with the given pixel size
metricSortInds = np.argsort(metric)
minIndMetric = 0
for i in np.arange(metricSortInds.shape[0]):
    fromTclean = False
    # assume we're using corrected output parameters
    source = ['{:f}arcsec'.format(corrMaj[metricSortInds[i]]),
              '{:f}arcsec'.format(corrMin[metricSortInds[i]]),
              '{:f}deg'.format(corrPA[metricSortInds[i]])]
    # check if tclean was run on these inputs
    for j in np.arange(tcleanInMaj.shape[0]):
        if (approxInMaj[metricSortInds[i]] == tcleanInMaj[j]) and \
           (approxInMin[metricSortInds[i]] == tcleanInMin[j]) and \
           (approxInPA[metricSortInds[i]] == tcleanInPA[j]):
            source = ['{:f}arcsec'.format(tcleanMaj[j]),
                      '{:f}arcsec'.format(tcleanMin[j]),
                      '{:f}deg'.format(tcleanPA[j])]
            fromTclean = True
            tcleanOutInd = j
    # check if post-taper convolution looks possible
    try:
        kernel = ia.beamforconvolvedsize(source=source, convolved=targetBm)
    except RuntimeError:
        continue
    # check if cell size will permit this convolution
    kernMin = kernel['minor']['value']
    pixDiag = 2*(float(cell.split('arc')[0])**2)
    if kernMin > pixDiag:
        minIndMetric = metricSortInds[i]
        break

# report findings
print 'Target:', targetBm
print "Optimal tapering input: " \
      + "['{}arcsec', '{}arcsec', '{}deg']".format(approxInMaj[minIndMetric],
                                                   approxInMin[minIndMetric],
                                                   approxInPA[minIndMetric])
if fromTclean:
    print "Exact tapering output (from tclean): " \
          + "['{}arcsec', '{}arcsec', '{}deg']".format(tcleanMaj[tcleanOutInd],
                                                       tcleanMin[tcleanOutInd],
                                                       tcleanPA[tcleanOutInd])
else:
    print "Estimated tapering output: " \
          + "['{}arcsec', '{}arcsec', '{}deg']".format(corrMaj[minIndMetric],
                                                       corrMin[minIndMetric],
                                                       corrPA[minIndMetric])
#------------------------------------------------------------------------------#
