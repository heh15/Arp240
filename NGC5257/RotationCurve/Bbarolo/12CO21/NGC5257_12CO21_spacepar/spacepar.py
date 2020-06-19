# A script to produce output plots for SPACEPAR task.
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
matplotlib.rc('xtick',direction='in')
matplotlib.rc('ytick',direction='in')
matplotlib.rc('font',family='sans-serif',serif='Helvetica',size=10)

f = fits.open('../NGC5257_12CO21_spacepar/spacepar_NGC_5257.fits')[0]
d, h = f.data, f.header
p1, p2 = h['CTYPE1'], h['CTYPE2']
p1u, p2u, ru = h['CUNIT1'], h['CUNIT2'], h['CUNIT3']
p1ran = (np.arange(0,d.shape[2])+1-h['CRPIX1'])*h['CDELT1']+h['CRVAL1']
p2ran = (np.arange(0,d.shape[1])+1-h['CRPIX2'])*h['CDELT2']+h['CRVAL2']
rings = (np.arange(0,d.shape[0])+1-h['CRPIX3'])*h['CDELT3']+h['CRVAL3']

nrad = d.shape[0]
ncols = 4
nrows = int(np.ceil(nrad/float(ncols)))
ext = [p1ran[0],p1ran[-1],p2ran[0],p2ran[-1]]
cmap = plt.get_cmap('nipy_spectral') #plt.get_cmap('gnuplot')

fig = plt.figure(figsize=(8,8))
x_axis_len, y_axis_len = 0.27, 0.27 
x_sep, y_sep = 0.07, 0.08 
count = 0
axis, bottom_corner = [], [0.1,0.7]
for i in range (nrows):
	bottom_corner[0] = 0.1
	for j in range (ncols):
		if (count>=nrad): break
		axis.append(fig.add_axes([bottom_corner[0],bottom_corner[1],x_axis_len,y_axis_len]))
		bottom_corner[0]+=x_axis_len+x_sep
		count += 1
	bottom_corner[1]-=(y_axis_len+y_sep)

for i in range (nrad):
	nr = int(i/ncols) + 1
	nc = i - (nr-1)*ncols + 1
	toplot = d[i]
	a = np.unravel_index(np.argmin(toplot),toplot.shape)
	p1min, p2min = p1ran[a[1]], p2ran[a[0]]
	ax = axis[i]
	ax.set_xlim(ext[0],ext[1])
	ax.set_ylim(ext[2],ext[3])
	ax.imshow(toplot,origin='lower',extent=ext,aspect='auto',cmap=cmap)
	ax.plot(p1min,p2min,'x',mew=2,ms=8,c='w')
	radstr = 'R = %.2f %s'%(rings[i],ru)
	minstr = 'min = (%.1f %s, %.1f %s)'%(p1min,p1u,p2min,p2u)
	ax.text(0.01,1.1,radstr,transform=ax.transAxes)
	ax.text(0.01,1.03,minstr,transform=ax.transAxes)
	if nc==1: ax.set_ylabel(p2+' ('+p2u+')')
	if (nr==nrows) or (nr==nrows-1 and nrad%ncols!=0 and nc>nrad%ncols):
		ax.set_xlabel(p1+' ('+p1u+')')

fig.savefig('../NGC5257_12CO21_spacepar/spacepar_NGC_5257.pdf',bbox_inches='tight')
