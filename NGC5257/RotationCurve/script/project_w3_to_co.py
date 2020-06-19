# Run this locally
# from ~/Dropbox/mac/wise_w3_vs_co
import glob
import os
import astropy
import numpy as np
# import matplotlib
# matplotlib.use("macOSX")
import matplotlib.pyplot as pl
import matplotlib.image as mpimg
from astropy import wcs
from astropy.wcs import WCS
from astropy.io import fits
from astropy.table import Table, Column
from astropy.cosmology import FlatLambdaCDM
from astropy import coordinates
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS

from PyPDF2 import PdfFileMerger

from reproject import reproject_interp
from reproject import reproject_exact

from skimage.transform import rescale, resize
from skimage.measure import block_reduce

import pandas as pd
import pickle
import reproject_califa as repc

from scipy.optimize import curve_fit

import linmix
from ltsfit.lts_linefit import lts_linefit

fun = lambda x, m, b: m * x + b


def galname_short2long(galname):
    gal = galname
    if galname[:3] == 'UGC':
        if len(
                glob.glob(
                    '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
                    % (galname, ))) == 0:
            gal = 'UGC0' + galname[3:]
        if len(
                glob.glob(
                    '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
                    % (gal, ))) == 0:
            gal = 'UGC00' + galname[3:]
    if galname[:2] == 'IC':
        if len(
                glob.glob(
                    '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
                    % (galname, ))) == 0:
            gal = 'IC0' + galname[2:]
    if galname[:3] == 'NGC':
        if len(
                glob.glob(
                    '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
                    % (galname, ))) == 0:
            gal = 'NGC0' + galname[3:]

    if len(
            glob.glob(
                '/Users/ryan/venus/shared_data/califa/DR3-stack/%s/%s_result.pk'
                % (gal, gal))) == 0:
        print("Nothing exists for this galaxy")
        return None
    return gal


def run_linmix(x, y, xerr, yerr, parallelize=True):
    # if parallelize==True:
    #     print("Parallel")
    #     print("x: ",x)
    #     print("xerr: ",xerr)
    #     print("y: ",y)
    #     print("yerr: ",yerr)
    lm_result = linmix.LinMix(x,
                              y,
                              xerr,
                              yerr,
                              K=3,
                              parallelize=parallelize,
                              nchains=8)
    lm_result.run_mcmc(silent=True)
    chains = np.vstack([lm_result.chain['alpha'], lm_result.chain['beta']])
    # print(np.average(chains[0]), np.average(chains[1]))
    result = dict()
    result['chains'] = chains
    result['intercept'] = np.average(chains[0])
    result['intercept_err'] = np.std(chains[0])
    result['slope'] = np.average(chains[1])
    result['slope_err'] = np.std(chains[1])
    return result


def plot_fit(t,
             a,
             b,
             a_err=0,
             b_err=0,
             s=None,
             pivot=0,
             ax=None,
             log=False,
             color='b',
             lw=2,
             alpha=0.5,
             **kwargs):
    """
    alpha is used to shade the uncertainties from a_err and b_err
    **kwargs is passed to pl.plot() for the central line only
    the error band has zorder=-10
    """
    if log:
        if pivot == 0:
            pivot = 1
        y = lambda A, B: 10**A * (t / pivot)**B
    else:
        y = lambda A, B: A + B * (t - pivot)
    if ax is None:
        ax = pl
    # the length may vary depending on whether it's a default color
    # (e.g., 'r' or 'orange') or an rgb(a) color, etc, but as far as
    # I can think of none of these would have length 2.
    if len(color) != 2:
        color = (color, color)
    print('in lnr.plot: color =', color)
    ax.plot(t, y(a, b), ls='-', color=color[0], lw=lw, **kwargs)
    if a_err != 0 or b_err != 0:
        # to make it compatible with either one or two values
        a_err = np.array([a_err]).flatten()
        b_err = np.array([b_err]).flatten()
        if a_err.size == 1:
            a_err = [a_err, a_err]
        if b_err.size == 1:
            b_err = [b_err, b_err]
        err = [
            y(a - a_err[0], b - b_err[0]),
            y(a - a_err[0], b + b_err[1]),
            y(a + a_err[1], b - b_err[0]),
            y(a + a_err[1], b + b_err[1])
        ]
        ylo = np.min(err, axis=0)
        yhi = np.max(err, axis=0)
        ax.fill_between(t,
                        ylo,
                        yhi,
                        color=color[1],
                        alpha=alpha,
                        lw=0,
                        edgecolor='none',
                        zorder=-10)
    if s:
        if log:
            ax.plot(t, (1 + s) * y(a, b), ls='--', color=color[0], lw=lw)
            ax.plot(t, y(a, b) / (1 + s), ls='--', color=color[0], lw=lw)
        else:
            ax.plot(t, y(a, b) + s, ls='--', color=color[0], lw=lw)
            ax.plot(t, y(a, b) - s, ls='--', color=color[0], lw=lw)
    return


def scatterplot(x,
                y,
                xlabel,
                ylabel,
                xerr=None,
                yerr=None,
                ax=None,
                label=None):
    if (xerr is not None) or (yerr is not None):
        pl.errorbar(x,
                    y,
                    xerr=xerr,
                    yerr=yerr,
                    marker='o',
                    markersize=4.,
                    markerfacecolor='none',
                    markeredgecolor='k',
                    markeredgewidth=0.5,
                    capsize=3,
                    linestyle='none',
                    label=label)
    else:
        pl.scatter(x,
                   y,
                   marker='o',
                   s=4.,
                   facecolor='none',
                   edgewidth=0.5,
                   linestyle='none',
                   label=label)
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)


# H0=70.2, Omega_M = 0.275, Lambda0 = 0.725

# New file I found on Tony Wong's github -- note there are different versions with different things
edge_co_params = pd.read_csv(
    '/Users/ryan/Dropbox/mac/tsinghua/edge_pydb-master/dat_glob/derived/build/EDGE_COparameters.csv',
    skiprows=18)
#use it for offsets of the images and inclinations
names_ecopars, xoffs, yoffs, inc_ecopars = np.array(
    edge_co_params["Name"]), np.array(edge_co_params[" xoff"]), np.array(
        edge_co_params[" yoff"]), np.array(edge_co_params[' Inc'])

cosmo = FlatLambdaCDM(H0=70.2, Om0=0.275, Tcmb0=2.725)
redshifts = np.loadtxt(
    '/Users/ryan/Downloads/stz349_supplemental_files/TableA1.csv',
    usecols=5,
    delimiter=',',
    skiprows=1)
names = np.loadtxt(
    '/Users/ryan/Downloads/stz349_supplemental_files/TableA1.csv',
    usecols=0,
    delimiter=',',
    skiprows=1,
    dtype=str)
bam = np.loadtxt('/Users/ryan/Downloads/stz349_supplemental_files/TableA1.csv',
                 usecols=1,
                 delimiter=',',
                 skiprows=1,
                 dtype=str)

# fnames_mom0 = glob.glob('/Users/ryan/venus/shared_data/edge/moment0_W3_26April19/*mom0.fits')
fnames_mom0 = glob.glob(
    '/Users/ryan/venus/shared_data/edge/moment0_W3_30April19/*mom0.fits')

# edge_list = glob.glob('/Users/ryan/venus/shared_data/edge/signal_cubes/*.co.cmmsk.fits')
# galnames = [ nm.split('/')[-1].split('.')[0] for nm in edge_list ]
edge_list = glob.glob('/Users/ryan/venus/shared_data/edge/W3/image_half/*')
galnames = [nm.split('/')[-1] for nm in edge_list]

califa_basic = fits.open('/Users/ryan/Dropbox/mac/CALIFA_1_MS_basic.fits')
califa_basic_names = califa_basic[1].data.field('REALNAME')
califa_basic_z = califa_basic[1].data.field('redshift')

califa_es = fits.open('/Users/ryan/Dropbox/mac/CALIFA_1_ES_basic.fits')
califa_es_names = califa_basic[1].data.field('REALNAME')
califa_es_z = califa_basic[1].data.field('redshift')

# galname = 'ARP220'
# galname = 'NGC5000'
# galname = 'NGC4211NED02'
# galname = 'NGC6361'
# galname = 'UGC10710'

leda = Table.read(
    '/usr/local/lib/python3.7/site-packages/edge_pydb/dat_glob/external/edge_leda.csv',
    format='ascii.ecsv')

# PA and INC from Becca's table
rfpars = Table.read(
    '/Users/ryan/Dropbox/mac/tsinghua/edge_pydb-master/dat_glob/derived/edge_rfpars.csv',
    format='ascii.ecsv')

edge_califa = Table.read(
    '/usr/local/lib/python3.7/site-packages/edge_pydb/dat_glob/external/edge_califa.csv',
    format='ascii.ecsv')
edge_califa_names = np.array(edge_califa['Name'])

gal_dict = dict()
Lsun = 3.839e33  # Lsun in erg/s

fname_mom0_jypb_sun_rebinned = lambda galname: '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_sun_rebin6_jypb_mom0.fits' % (
    galname, )
fname_noise_map = lambda galname: '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_noise_jypbkms.fits' % (
    galname, )
fname_noise_map_total = lambda galname: '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_total_noise_rebin6_jypbkms.fits' % (
    galname, )

# def popen(galname,label):
#     with open('/Users/ryan/venus/home/edgecalifa_processed/%s_%s.pkl'%(galname,label)) as f:
#         lst = pickle.load(f)
#     return lst


def get_3color(galname, ra, dec):
    w = wcs.WCS(naxis=2)
    gal = galname
    if galname == 'NGC4211NED02':
        gal = 'NGC4211B'
    # c = coordinates.get_icrs_coordinates(gal) # center of galaxy coordinate
    linwcs = lambda x, y, n: ((x - y) / n, (x + y) / 2)
    cdeltaX, crvalX = linwcs(ra - 1.2 / 60. / 2., ra + 1.2 / 60. / 2., 512)
    cdeltaY, crvalY = linwcs(dec - 1.2 / 60. / 2., dec + 1.2 / 60. / 2., 512)
    # what is the center pixel of the XY grid.
    w.wcs.crpix = [512. / 2, 512. / 2]
    # what is the galactic coordinate of that pixel.
    w.wcs.crval = [crvalX, crvalY]
    # what is the pixel scale in lon, lat.
    w.wcs.cdelt = np.array([cdeltaX, cdeltaY])
    w.wcs.ctype = ['RA---SIN', 'DEC--SIN']
    # ax = pl.subplot(2,5,1,projection=w)
    image = mpimg.imread(
        glob.glob('/Users/ryan/venus/home/edge_sdss_thumbnails/%s*.jpg' %
                  (galname, ))[0])
    return image, w


def do_all_fits(x, y, x_err, y_err, parallelize=True):
    temp_dir = '/Users/ryan/Dropbox/mac/wise_w3_vs_co/tmp/'
    pfit_flag = 0
    lts_fit_reverse = -1
    lts_fit_forward = -1
    fit_mask = -1
    fit_result = -1
    fit_result_reverse = -1
    fit_result_forward = -1
    fit_result_reverse_masked = -1
    fit_result_forward_masked = -1
    cfit = -1
    pfit = -1
    if x.size > 6:
        try:
            lts_fit_reverse = lts_linefit(y,
                                          x,
                                          y_err,
                                          x_err,
                                          plot=True,
                                          fontsize=8)
            fit_mask_reverse = lts_fit_reverse.mask
        except:
            print("LTS reverse failed")
            fit_mask_reverse = None
            lts_fit_reverse = -1

        # Pickle the result to save memory
        lts_reverse_tmp = dict()
        lts_reverse_tmp['lts_fit_reverse'] = lts_fit_reverse
        lts_reverse_tmp['fit_mask_reverse'] = fit_mask_reverse
        with open(temp_dir + 'lts_fit_reverse.pk', 'wb') as p:
            pickle.dump(lts_reverse_tmp, p)

        del lts_reverse_tmp, lts_fit_reverse

        try:
            lts_fit_forward = lts_linefit(x, y, x_err, y_err, plot=False)
            fit_mask_forward = lts_fit_forward.mask
        except:
            print("LTS forward failed")
            fit_mask_forward = None
            lts_fit_forward = -1

        # Pickle the result to save memory
        lts_forward_tmp = dict()
        lts_forward_tmp['lts_fit_forward'] = lts_fit_forward
        lts_forward_tmp['fit_mask_forward'] = fit_mask_forward
        with open(temp_dir + 'lts_fit_forward.pk', 'wb') as p:
            pickle.dump(lts_forward_tmp, p)
        del lts_forward_tmp, lts_fit_forward

        # Linmix fit on unmasked
        fit_result_reverse = run_linmix(y,
                                        x,
                                        y_err,
                                        x_err,
                                        parallelize=parallelize)

        # Pickle to save memory
        fit_result_reverse = {'fit_result_reverse': fit_result_reverse}
        with open(temp_dir + 'fit_result_reverse.pk', 'wb') as p:
            pickle.dump(fit_result_reverse, p)
        del fit_result_reverse

        fit_result_forward = run_linmix(x,
                                        y,
                                        x_err,
                                        y_err,
                                        parallelize=parallelize)
        # Pickle to save memory
        fit_result_forward = {'fit_result_forward': fit_result_forward}
        with open(temp_dir + 'fit_result_forward.pk', 'wb') as p:
            pickle.dump(fit_result_forward, p)
        del fit_result_forward

        # Construct mask
        if (fit_mask_forward is None) and (fit_mask_reverse is None):
            print("Both LTS fits failed")
            fit_mask = -1
        if (fit_mask_forward is not None) and (fit_mask_reverse is None):
            fit_mask = fit_mask_forward
        if (fit_mask_forward is not None) and (fit_mask_reverse is not None):
            fit_mask = fit_mask_forward & fit_mask_reverse
        if (fit_mask_forward is None) and (fit_mask_reverse is not None):
            fit_mask = fit_mask_reverse

        # Pass non-outlier points to linmix
        if type(fit_mask) == np.ndarray:
            fit_result_forward_masked = run_linmix(x[fit_mask],
                                                   y[fit_mask],
                                                   x_err[fit_mask],
                                                   y_err[fit_mask],
                                                   parallelize=parallelize)
            # Pickle to save memory
            fit_result_forward_masked = {
                'fit_result_forward_masked': fit_result_forward_masked
            }
            with open(temp_dir + 'fit_result_forward_masked.pk', 'wb') as p:
                pickle.dump(fit_result_forward_masked, p)
            del fit_result_forward_masked

            fit_result_reverse_masked = run_linmix(y[fit_mask],
                                                   x[fit_mask],
                                                   y_err[fit_mask],
                                                   x_err[fit_mask],
                                                   parallelize=parallelize)
            # Pickle to save memory
            fit_result_reverse_masked = {
                'fit_result_reverse_masked': fit_result_reverse_masked
            }
            with open(temp_dir + 'fit_result_reverse_masked.pk', 'wb') as p:
                pickle.dump(fit_result_reverse_masked, p)
            del fit_result_reverse_masked

        else:
            fit_result_reverse_masked = -1
            fit_result_forward_masked = -1

            # Pickle to save memory
            fit_result_forward_masked = {
                'fit_result_forward_masked': fit_result_forward_masked
            }
            with open(temp_dir + 'fit_result_forward_masked.pk', 'wb') as p:
                pickle.dump(fit_result_forward_masked, p)
            del fit_result_forward_masked

            # Pickle to save memory
            fit_result_reverse_masked = {
                'fit_result_reverse_masked': fit_result_reverse_masked
            }
            with open(temp_dir + 'fit_result_reverse_masked.pk', 'wb') as p:
                pickle.dump(fit_result_reverse_masked, p)
            del fit_result_reverse_masked

        pfit_flag = 1  # did do fit
        cfit = -1
        pfit = -1

        # Load everything back again
        lts_fit_reverse = pickle.load(
            open(temp_dir + 'lts_fit_reverse.pk', 'rb'))['lts_fit_reverse']
        lts_fit_forward = pickle.load(
            open(temp_dir + 'lts_fit_forward.pk', 'rb'))['lts_fit_forward']
        fit_result_reverse = pickle.load(
            open(temp_dir + 'fit_result_reverse.pk',
                 'rb'))['fit_result_reverse']
        fit_result_forward = pickle.load(
            open(temp_dir + 'fit_result_forward.pk',
                 'rb'))['fit_result_forward']
        fit_result_reverse_masked = pickle.load(
            open(temp_dir + 'fit_result_reverse_masked.pk',
                 'rb'))['fit_result_reverse_masked']
        fit_result_forward_masked = pickle.load(
            open(temp_dir + 'fit_result_forward_masked.pk',
                 'rb'))['fit_result_forward_masked']

        os.system('rm -rf %s*' % (temp_dir, ))

    # Save to final dictionary
    result = dict()
    result['pfit_flag'] = pfit_flag
    result['lts_fit_reverse'] = lts_fit_reverse
    result['lts_fit_forward'] = lts_fit_forward
    result['fit_mask'] = fit_mask
    result['fit_result'] = fit_result
    result['fit_result_reverse'] = fit_result_reverse
    result['fit_result_forward'] = fit_result_forward
    result['fit_result_reverse_masked'] = fit_result_reverse_masked
    result['fit_result_forward_masked'] = fit_result_forward_masked
    result['cfit'] = cfit
    result['pfit'] = pfit
    return result


def do_gal(galname, gal_dict, skip_fitting=False):
    gal = galname
    vstring = '_co_smooth_wise_v2_rebin6_mom0_Sun_v2.K.fits'
    gal = galname_short2long(galname)

    # if galname[:3] == 'UGC':
    #     if len(
    #             glob.glob(
    #                 '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
    #                 % (galname, ))) == 0:
    #         gal = 'UGC0' + galname[3:]
    #     if len(
    #             glob.glob(
    #                 '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
    #                 % (gal, ))) == 0:
    #         gal = 'UGC00' + galname[3:]
    # if galname[:2] == 'IC':
    #     if len(
    #             glob.glob(
    #                 '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
    #                 % (galname, ))) == 0:
    #         gal = 'IC0' + galname[2:]
    # if galname[:3] == 'NGC':
    #     if len(
    #             glob.glob(
    #                 '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
    #                 % (galname, ))) == 0:
    #         gal = 'NGC0' + galname[3:]
    if (gal is None) or (gal in ['UGC05359', 'UGC03253']):
        return gal_dict
    else:
        # mom0 = fits.open('/Users/ryan/venus/shared_data/edge/co_cubes_w3conv_aniano_6Jun19/%s%s'%(gal,vstring))
        # mom0 = fits.open('/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'%(gal,))

        mom0 = fits.open(fname_mom0_jypb_sun_rebinned(gal))

        # mom0_map = mom0[0].data[0,0,:,:]

        # mdat, mhdr = fits.getdata('/Users/ryan/venus/shared_data/edge/co_cubes_w3conv_aniano_6Jun19/%s%s'%(gal,vstring_v2), header=True)
        mdat, mhdr = fits.getdata(
            '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
            % (gal, ),
            header=True)

        alpha_co = 3.2

        beam_area = 2 * np.pi * 6.6**2 / (4 * 2 * np.log(2))
        pix_area = 6**2
        pix_per_beam = beam_area / pix_area

        mom0_map = alpha_co * mom0[0].data[0] / pix_per_beam
        mom0_error = alpha_co * fits.open(
            fname_noise_map_total(gal))[0].data[0] / pix_per_beam

        # Make 2 other mom0 maps:
        # 1. Using metallicity-dependent alpha_co (only over star-forming spaxels)
        # 2. Using constant alpha_co but only over star-forming spaxels
        alpha_co_met_sf_only = pickle.load(
            open(repc.fname_alpha_co_stacked(gal), 'rb'))['alpha_co']
        alpha_co_sf_only = np.ones(alpha_co_met_sf_only.shape) * alpha_co
        alpha_co_sf_only[np.isnan(alpha_co_met_sf_only)] = np.nan

        mom0_map_met_sf_only = alpha_co_met_sf_only * (mom0[0].data[0] /
                                                       pix_per_beam)
        mom0_map_met_sf_only_err = alpha_co_met_sf_only * fits.open(
            fname_noise_map_total(gal))[0].data[0] / pix_per_beam
        mom0_map_sf_only = alpha_co_sf_only * (mom0[0].data[0] / pix_per_beam)
        mom0_map_sf_only_err = alpha_co_sf_only * fits.open(
            fname_noise_map_total(gal))[0].data[0] / pix_per_beam

        mom0_shape = mom0_map.shape

        w = wcs.WCS(mhdr)
        # w.wcs.cdelt = np.array([mom0[0].header['CDELT1']*6 , mom0[0].header['CDELT2']*6 ])
        # w.wcs.crpix = w.wcs.crpix/6.
        # w.wcs.crpix = [x, y]

        idx = np.where(leda['Name'] == gal)[0][0]
        ractr = leda['ledaRA'].quantity[idx].to(u.deg).value
        dcctr = leda['ledaDE'].quantity[idx].to(u.deg).value
        r25 = leda['ledaD25'].quantity[idx].to(u.arcsec).value / 2.
        # dist_pc = leda['ledaDistMpc'].quantity[idx].to(u.Mpc).value * 1e6
        dist_pc = edge_califa['caDistMpc'][edge_califa_names == gal].to(
            u.pc).value
        dist_m = (dist_pc * u.pc).to(u.m).value  # * 3.086e16

        x, y = w.wcs_world2pix(ractr, dcctr, 0)

        idx = np.where(rfpars['Name'] == gal)[0][0]
        pa = rfpars['rfPA'][idx]
        inc = rfpars['rfInc'][idx]
        print("Galaxy: " + gal)
        print("Inclination: " + str(inc))
        print("Deprojection factor: " + str(np.cos(np.radians(inc))))

        wcs_out = wcs.WCS(naxis=2)
        wcs_out.wcs.crpix = [x, y]
        # what is the galactic coordinate of that pixel.
        wcs_out.wcs.crval = [ractr, dcctr]
        # what is the pixel scale in lon, lat.
        wcs_out.wcs.cdelt = np.array([mhdr['CDELT1'], mhdr['CDELT2']])
        wcs_out.wcs.ctype = ['RA---SIN', 'DEC--SIN']

        w = wcs.WCS(naxis=2)

        image_w3 = fits.open(
            '/Users/ryan/venus/shared_data/edge/W3/image_half/%s/final_36_flux.fits'
            % (galname, ))[0].data
        uncertainty_w3 = fits.open(
            '/Users/ryan/venus/shared_data/edge/W3/image_half/%s/uncer_36_flux.fits'
            % (galname, ))[0].data
        image_w3_v2 = fits.open(
            '/Users/ryan/venus/shared_data/edge/W3/image_half1/%s/final_36_flux.fits'
            % (galname, ))[1].data
        ra_w3, dec_w3 = image_w3_v2['RA'].reshape(
            image_w3.shape), image_w3_v2['DEC'].reshape(image_w3.shape)

        pix_area_pc2 = (6. * dist_pc / 206265.)**2

        z = edge_califa['caZgas'][edge_califa_names == gal]
        mom0_map *= 2453. * (dist_pc / 1e6)**2 * np.cos(
            np.radians(inc)) / (1. + z) / pix_area_pc2
        mom0_error *= 2453. * (dist_pc / 1e6)**2 * np.cos(
            np.radians(inc)) / (1. + z) / pix_area_pc2

        mom0_map_met_sf_only *= 2453. * (dist_pc / 1e6)**2 * np.cos(
            np.radians(inc)) / (1. + z) / pix_area_pc2
        mom0_map_sf_only *= 2453. * (dist_pc / 1e6)**2 * np.cos(
            np.radians(inc)) / (1. + z) / pix_area_pc2
        mom0_map_met_sf_only_err *= 2453. * (dist_pc / 1e6)**2 * np.cos(
            np.radians(inc)) / (1. + z) / pix_area_pc2
        mom0_map_sf_only_err *= 2453. * (dist_pc / 1e6)**2 * np.cos(
            np.radians(inc)) / (1. + z) / pix_area_pc2

        nx, ny = image_w3.shape

        # what is the galactic coordinate of that pixel.
        w.wcs.crval = [
            ra_w3[int((nx) / 2.), int((ny) / 2.)], dec_w3[int(
                (nx) / 2.), int((ny) / 2.)]
        ]
        w.wcs.crpix = [int((nx) / 2.), int((ny) / 2.)]
        # what is the pixel scale in lon, lat.
        w.wcs.cdelt = np.array([-1.375, 1.375]) / 3600.
        # w.wcs.cdelt = np.array([-6., 6.])/3600.
        w.wcs.ctype = ['RA---SIN', 'DEC--SIN']

        wcs_out_tmp = wcs.WCS(naxis=2)
        wcs_out_tmp.wcs.crpix = [int((nx) / 2.), int((ny) / 2.)]
        wcs_out_tmp.wcs.crval = [ractr, dcctr]
        # what is the pixel scale in lon, lat.
        wcs_out_tmp.wcs.cdelt = np.array([-1.375, 1.375]) / 3600.
        wcs_out_tmp.wcs.ctype = ['RA---SIN', 'DEC--SIN']

        wcs_wise = wcs.WCS(naxis=2)
        wcs_wise.wcs.crpix = [
            int((nx) / 2.) + 1, int((ny) / 2.) + 1
        ]  #[mom0[0].header['CRPIX1'],mom0[0].header['CRPIX2']]
        wcs_wise.wcs.crval = [
            ra_w3[int((nx) / 2.), int((ny) / 2.)], dec_w3[int(
                (nx) / 2.), int((ny) / 2.)]
        ]
        wcs_wise.wcs.cdelt = np.array([-1.375, 1.375]) / 3600.
        wcs_wise.wcs.ctype = ['RA---SIN', 'DEC--SIN']

        image_w3 *= (dist_pc / 1e6)**2 * np.cos(
            np.radians(inc)) * 10**0.8477 / pix_area_pc2
        image_w3_proj, footprint = reproject_exact((image_w3, wcs_wise),
                                                   wcs.WCS(mhdr),
                                                   mom0_map.shape,
                                                   parallel=False)
        pixel_size_arcsec_input = wcs_wise.wcs.cdelt[1]
        pixel_size_arcsec_output = wcs.WCS(mhdr).wcs.cdelt[1]
        image_w3_proj *= (pixel_size_arcsec_output /
                          pixel_size_arcsec_input)**2

        # Do the same for uncertainty
        uncertainty_w3 *= (dist_pc / 1e6)**2 * np.cos(
            np.radians(inc)) * 10**0.8477 / pix_area_pc2

        uncertainty_w3_proj, footprint = reproject_exact(
            (uncertainty_w3**2, wcs_wise),
            wcs.WCS(mhdr),
            mom0_map.shape,
            parallel=False)
        uncertainty_w3_proj = np.sqrt(uncertainty_w3_proj)
        uncertainty_w3_proj *= (pixel_size_arcsec_output /
                                pixel_size_arcsec_input)**2

        uncertainty_w3_proj = np.sqrt(5 * uncertainty_w3_proj**2 +
                                      (0.0414 * image_w3_proj)**2)

        pl.figure(figsize=(22, 4))

        sdss_3color = get_3color(gal, ractr, dcctr)

        ax = pl.subplot(141, projection=sdss_3color[1])
        ax.imshow(sdss_3color[0])
        ax.set_autoscale_on(False)
        pl.title(galname)

        ax = pl.subplot(142, projection=wcs_out)
        ax.imshow(mom0_map[:, ::-1])
        pl.title("EDGE CO(1-0) map, W3 resolution")

        ax = pl.subplot(143, projection=wcs_out)
        image_w3_proj[np.isnan(image_w3_proj)] = 0.
        ax.imshow(image_w3_proj[:, ::-1], cmap='hot')
        pl.title(r"$WISE$ W3 map")
        pl.subplot(144)

        # Initialize x's and y's (and their errors) for fitting
        ys_deg = image_w3_proj.flatten()
        ys_err = 0.434 * uncertainty_w3_proj.flatten() / ys_deg

        # Sigma_H2 using constant alpha_co, all pixels
        xs_deg = mom0_map.flatten()
        xs_err = 0.434 * mom0_error.flatten() / xs_deg

        # Sigma_H2 using metallicity-dependent alpha_co, SF pixels only
        xs_deg_met_sf_only = mom0_map_met_sf_only.flatten()
        xs_deg_met_sf_only_err = 0.434 * mom0_map_met_sf_only_err.flatten(
        ) / xs_deg_met_sf_only

        # Sigma_H2 using constant alpha_co, SF pixels only
        xs_deg_sf_only = mom0_map_sf_only.flatten()
        xs_deg_sf_only_err = 0.434 * mom0_map_sf_only_err.flatten(
        ) / xs_deg_sf_only

        # Points using metallicity-dependent alpha_co, SF spaxels only
        x_fit_met = np.log10(xs_deg_met_sf_only[xs_deg_met_sf_only > 0])
        x_fit_met_err = xs_deg_met_sf_only_err[xs_deg_met_sf_only > 0]
        y_fit_met = np.log10(ys_deg[xs_deg_met_sf_only > 0])
        y_fit_met_err = ys_err[xs_deg_met_sf_only > 0]

        x_sum_met = np.log10(np.sum(10**x_fit_met))
        x_sum_met_err = 0.434 * np.sqrt(
            np.sum(mom0_map_met_sf_only_err.flatten()[xs_deg_met_sf_only > 0]**
                   2)) / 10**x_sum_met
        y_sum_met = np.log10(np.sum(10**y_fit_met))
        y_sum_met_err = 0.434 * np.sqrt(
            np.sum(uncertainty_w3_proj.flatten()[xs_deg_met_sf_only > 0]**
                   2)) / 10**y_sum_met

        good_met = (~np.isinf(x_fit_met)) & (~np.isinf(y_fit_met)) & (
            ~np.isinf(x_fit_met_err)) & (~np.isinf(y_fit_met_err)) & (
                ~np.isnan(x_fit_met)) & (~np.isnan(y_fit_met)) & (
                    ~np.isnan(x_fit_met_err)) & (~np.isnan(y_fit_met_err))

        x_fit_met = x_fit_met[good_met]
        x_fit_met_err = x_fit_met_err[good_met]
        y_fit_met = y_fit_met[good_met]
        y_fit_met_err = y_fit_met_err[good_met]

        # Points using constant alpha_co but only on SFing spaxels
        x_fit_sf = np.log10(xs_deg_sf_only[xs_deg_sf_only > 0])
        x_fit_sf_err = xs_deg_sf_only_err[xs_deg_sf_only > 0]
        y_fit_sf = np.log10(ys_deg[xs_deg_sf_only > 0])
        y_fit_sf_err = ys_err[xs_deg_sf_only > 0]

        x_sum_sf = np.log10(np.sum(10**x_fit_sf))
        x_sum_sf_err = 0.434 * np.sqrt(
            np.sum(mom0_map_sf_only_err.flatten()[xs_deg_sf_only > 0]**
                   2)) / 10**x_sum_sf
        y_sum_sf = np.log10(np.sum(10**y_fit_sf))
        y_sum_sf_err = 0.434 * np.sqrt(
            np.sum(uncertainty_w3_proj.flatten()[xs_deg_sf_only > 0]**
                   2)) / 10**y_sum_sf

        good_sf = (~np.isinf(x_fit_sf)) & (~np.isinf(y_fit_sf)) & (
            ~np.isinf(x_fit_sf_err)) & (~np.isinf(y_fit_sf_err)) & (
                ~np.isnan(x_fit_sf)) & (~np.isnan(y_fit_sf)) & (
                    ~np.isnan(x_fit_sf_err)) & (~np.isnan(y_fit_sf_err))

        x_fit_sf = x_fit_sf[good_sf]
        x_fit_sf_err = x_fit_sf_err[good_sf]
        y_fit_sf = y_fit_sf[good_sf]
        y_fit_sf_err = y_fit_sf_err[good_sf]

        # Points using constant alpha_co on all spaxels

        X = np.log10(xs_deg[xs_deg > 0])
        X_err = xs_err[xs_deg > 0]
        Y = np.log10(ys_deg[xs_deg > 0])
        Y_err = ys_err[xs_deg > 0]
        Xsum = np.log10(np.sum(10**(X)))
        Ysum = np.log10(np.sum(10**(Y)))

        Xsum_err = 0.434 * np.sqrt(np.sum(mom0_error.flatten()[xs_deg > 0]**
                                          2)) / 10**Xsum
        Ysum_err = 0.434 * np.sqrt(
            np.sum(uncertainty_w3_proj.flatten()[xs_deg > 0]**2)) / 10**Ysum

        pfit = np.poly1d(np.ones(2))
        pfit_flag = 0  # Didn't do fit

        good = (~np.isinf(X)) & (~np.isinf(Y)) & (~np.isinf(X_err)) & (
            ~np.isinf(Y_err)) & (~np.isnan(X)) & (~np.isnan(Y)) & (
                ~np.isnan(X_err)) & (~np.isnan(Y_err))

        Xf = X[good]
        Yf = Y[good]
        Xf_err = X_err[good]
        Yf_err = Y_err[good]

        if skip_fitting == False:
            result_fixed_alpha = do_all_fits(Xf, Yf, Xf_err, Yf_err)

            result_met_alpha_sf_only = do_all_fits(x_fit_met, y_fit_met,
                                                   x_fit_met_err,
                                                   y_fit_met_err)

            result_fixed_alpha_sf_only = do_all_fits(x_fit_sf, y_fit_sf,
                                                     x_fit_sf_err,
                                                     y_fit_sf_err)

        # pfit_flag = 0
        # lts_fit_reverse = -1
        # lts_fit_forward = -1
        # fit_mask = -1
        # fit_result = -1
        # fit_result_reverse = -1
        # fit_result_forward = -1
        # fit_result_reverse_masked = -1
        # fit_result_forward_masked = -1
        # cfit = -1
        # pfit = -1
        # if Xf.size > 6:
        #     try:
        #         lts_fit_reverse = lts_linefit(Yf,
        #                                       Xf,
        #                                       Yf_err,
        #                                       Xf_err,
        #                                       plot=True,
        #                                       fontsize=8)
        #         fit_mask_reverse = lts_fit_reverse.mask
        #     except:
        #         print("LTS reverse failed")
        #         fit_mask_reverse = None
        #         lts_fit_reverse = -1
        #
        #     try:
        #         lts_fit_forward = lts_linefit(Xf,
        #                                       Yf,
        #                                       Xf_err,
        #                                       Yf_err,
        #                                       plot=False)
        #         fit_mask_forward = lts_fit_forward.mask
        #     except:
        #         print("LTS forward failed")
        #         fit_mask_forward = None
        #         lts_fit_forward = -1
        #
        #     # Construct mask
        #     if (fit_mask_forward is None) and (fit_mask_reverse is None):
        #         print("Both LTS fits failed")
        #         fit_mask = -1
        #     if (fit_mask_forward is not None) and (fit_mask_reverse is None):
        #         fit_mask = fit_mask_forward
        #     if (fit_mask_forward is not None) and (fit_mask_reverse is
        #                                            not None):
        #         fit_mask = fit_mask_forward & fit_mask_reverse
        #     if (fit_mask_forward is None) and (fit_mask_reverse is not None):
        #         fit_mask = fit_mask_reverse
        #
        #     # Linmix fit on unmasked
        #     fit_result_reverse = run_linmix(Yf,
        #                                     Xf,
        #                                     Yf_err,
        #                                     Xf_err,
        #                                     parallelize=False)
        #     fit_result_forward = run_linmix(Xf,
        #                                     Yf,
        #                                     Xf_err,
        #                                     Yf_err,
        #                                     parallelize=False)
        #
        #     # Pass non-outlier points to linmix
        #     if type(fit_mask) == np.ndarray:
        #         fit_result_forward_masked = run_linmix(Xf[fit_mask],
        #                                                Yf[fit_mask],
        #                                                Xf_err[fit_mask],
        #                                                Yf_err[fit_mask],
        #                                                parallelize=False)
        #         fit_result_reverse_masked = run_linmix(Yf[fit_mask],
        #                                                Xf[fit_mask],
        #                                                Yf_err[fit_mask],
        #                                                Xf_err[fit_mask],
        #                                                parallelize=False)
        #     else:
        #         fit_result_reverse_masked = -1
        #         fit_result_forward_masked = -1
        #
        #     pfit_flag = 1  # did do fit
        #     cfit = -1
        #     pfit = -1

        # gal_dict[gal] = {'LCO_map' : mom0_map,\ # L_CO map in units of K km/s pc^2
        # 'Lwise_map' : array,\ # L_12um map in units of Lsun
        # 'log10_Lwise' : Y,\ # Array of log L_12um (Lsun)
        # 'log10_LCO' : X,\ # Array of log L_CO (K km/s pc^2)
        # 'incl' : inc,\ # inclination angle in degrees
        # 'dist_Mpc' : dist_pc/1.e6,\ # Distance in Mpc
        # 'log10_Lwise_tot' : Ysum,\ # log integrated value of L_12um (Lsun)
        # 'log10_LCO_tot' : Xsum} # log integrated value of L_CO (K km/s pc^2)

        # 'LCO_map' # L_CO map in units of K km/s pc^2
        # 'Lwise_map' # L_12um map in units of Lsun
        # 'log10_Lwise' # Array of log L_12um (Lsun)
        # 'log10_LCO' # Array of log L_CO (K km/s pc^2)
        # 'incl'  # inclination angle in degrees
        # 'dist_Mpc'  # Distance in Mpc
        # 'log10_Lwise_tot'  # log integrated value of L_12um (Lsun)
        # 'log10_LCO_tot'  # log integrated value of L_CO (K km/s pc^2)
        # result_fixed_alpha = do_all_fits(Xf, Yf, Xf_err, Yf_err)
        #
        # result_met_alpha_sf_only = do_all_fits(x_fit_met, y_fit_met, x_fit_met_err, y_fit_met_err)
        #
        # result_fixed_alpha_sf_only = do_all_fits(x_fit_sf, y_fit_sf, x_fit_sf_err, y_fit_sf_err)

        if skip_fitting == False:
            gal_dict[gal] = {'LCO_map' : mom0_map ,\
            'LCO_map_err' : mom0_error ,\
            'Lwise_map_err' : uncertainty_w3_proj ,\
            'Lwise_map' : image_w3_proj ,\
            'log10_Lwise' : Y ,\
            'log10_Lwise_err' : Y_err ,\
            'log10_LCO' : X ,\
            'log10_LCO_err' : X_err ,\
            'incl' : inc ,\
            'dist_Mpc' : dist_pc / 1.e6 ,\
            'log10_Lwise_tot' : Ysum ,\
            'log10_Lwise_tot_err' : Ysum_err ,\
            'log10_LCO_tot' : Xsum ,\
            'log10_LCO_tot_err' : Xsum_err ,\
            'polyfit' : result_fixed_alpha['pfit'] ,\
            'polyfit_flag' : result_fixed_alpha['pfit_flag'],\
            'linmix_fit_reverse' : result_fixed_alpha['fit_result_reverse'],\
            'linmix_fit_forward' : result_fixed_alpha['fit_result_forward'],\
            'linmix_fit_reverse_masked' : result_fixed_alpha['fit_result_reverse_masked'],\
            'linmix_fit_forward_masked' : result_fixed_alpha['fit_result_forward_masked'],\
            'lts_fit_reverse' : result_fixed_alpha['lts_fit_reverse'],\
            'lts_fit_forward' : result_fixed_alpha['lts_fit_forward'],\
            'lts_fit_mask' : result_fixed_alpha['fit_mask'],\
            'x_fit' : Xf,\
            'y_fit' : Yf,\
            'x_err_fit' : Xf_err,\
            'y_err_fit' : Yf_err,\
            'fits_fixed_alpha_sf_only' : result_fixed_alpha_sf_only ,\
            'x_fit_sf' : x_fit_sf ,\
            'y_fit_sf' : y_fit_sf ,\
            'x_fit_sf_err': x_fit_sf_err ,\
            'y_fit_sf_err' : y_fit_sf_err ,\
            'result_met_alpha_sf_only' : result_met_alpha_sf_only ,\
            'x_fit_met' : x_fit_met ,\
            'y_fit_met' : y_fit_met ,\
            'x_fit_met_err': x_fit_met_err ,\
            'y_fit_met_err' : y_fit_met_err
            }

            if (Xf.size <= 6) or (type(result_fixed_alpha['fit_mask']) !=
                                  np.array):
                pl.errorbar(Y,
                            X,
                            xerr=Y_err,
                            yerr=X_err,
                            markersize=4.,
                            markeredgecolor='k',
                            capsize=3,
                            linestyle='none')
            pl.scatter(Ysum, Xsum, marker='o', facecolor='none', c='b')

            if result_fixed_alpha['pfit_flag'] == 1:
                plot_fit(np.linspace(-2, 4, 10),
                         result_fixed_alpha['fit_result_reverse']['intercept'],
                         result_fixed_alpha['fit_result_reverse']['slope'],
                         a_err=result_fixed_alpha['fit_result_reverse']
                         ['intercept_err'],
                         b_err=result_fixed_alpha['fit_result_reverse']
                         ['slope_err'])
                # for i_chain in np.linspace(0, fit_result['chains'][0].size - 1,
                #                            100).astype(int):
                #     pl.plot(
                #         np.linspace(-2, 4, 10),
                #         fit_result['chains'][0, i_chain] +
                #         np.linspace(-2, 4, 10) * fit_result['chains'][1, i_chain],
                #         color='r',
                #         alpha=0.1)
                #
            # pl.ylim(-2, 4)
            # pl.ylim(-2, 3)
            # pl.xlabel(r"Integrated CO intensity [Jy/beam km/s]")
            # pl.xlabel(r"log(CO luminosity) [K km/s pc$^2$/pixel]")
            # pl.ylabel(r"log($L_{WISE\> W3}$) [erg/s/pixel]")
            pl.ylabel(r"$\log \> \Sigma(\mathrm{H_2})$ [$M_\odot$ pc$^{-2}$]")
            pl.xlabel(
                r"$\log \> \Sigma(12\mu\mathrm{m})$ [L$_\odot$ pc$^{-2}$]")
            pl.title(galname + '; incl.=' + str(inc))
            pl.tight_layout()
            # pl.savefig('/Users/ryan/Dropbox/mac/wise_w3_vs_co/figs_v2/%s.pdf'%(galname,))
            pl.savefig('/Users/ryan/Dropbox/mac/wise_w3_vs_co/figs_v3/%s.pdf' %
                       (galname, ))
            # print("Done galaxy")
            pl.close()

        else:
            gal_dict[gal] = {'LCO_map' : mom0_map ,\
            'LCO_map_err' : mom0_error ,\
            'Lwise_map_err' : uncertainty_w3_proj ,\
            'Lwise_map' : image_w3_proj ,\
            'log10_Lwise' : Y ,\
            'log10_Lwise_err' : Y_err ,\
            'log10_LCO' : X ,\
            'log10_LCO_err' : X_err ,\
            'incl' : inc ,\
            'dist_Mpc' : dist_pc / 1.e6 ,\
            'log10_Lwise_tot' : Ysum ,\
            'log10_Lwise_tot_err' : Ysum_err ,\
            'log10_LCO_tot' : Xsum ,\
            'log10_LCO_tot_err' : Xsum_err ,\
            'x_fit' : Xf,\
            'y_fit' : Yf,\
            'x_err_fit' : Xf_err,\
            'y_err_fit' : Yf_err,\
            'x_fit_sf' : x_fit_sf ,\
            'y_fit_sf' : y_fit_sf ,\
            'x_fit_sf_err': x_fit_sf_err ,\
            'y_fit_sf_err' : y_fit_sf_err ,\
            'x_fit_met' : x_fit_met ,\
            'y_fit_met' : y_fit_met ,\
            'x_fit_met_err': x_fit_met_err ,\
            'y_fit_met_err' : y_fit_met_err
            }

    return gal_dict


fname_l12_lco_dict = '/Users/ryan/Dropbox/mac/wise_w3_vs_co/l12_lco_dict_v5.pk'  # Without Aniano kernels, until that is fixed, and with more conservative s/n cuts in mom0 map making
fname_l12_lco_dict_nofits = '/Users/ryan/Dropbox/mac/wise_w3_vs_co/l12_lco_dict_nofits_v5.pk'  # Without Aniano kernels, until that is fixed, and with more conservative s/n cuts in mom0 map making

# gal_dict = pickle.load(open(fname_l12_lco_dict, 'rb'))

if __name__ == '__main__':
    import sys
    import os

    if sys.argv[1] == 'gal_dict':
        pick_up_from = None
        if len(sys.argv) > 2:
            caught_up = False
            pick_up_from = sys.argv[2]
            if pick_up_from == 'skip_fit':
                print("Skipping fitting")
            else:
                print("Resuming from %s" % (pick_up_from, ))

        i = 1
        n = len(galnames)
        for galname in galnames:  #['UGC10043']: #galnames:# ['IC2487']:
            gal = galname_short2long(galname)
            if gal == None:
                i += 1
                continue

            # UGC4029

            if (pick_up_from is not None) and (pick_up_from != 'skip_fit'):
                if galname == pick_up_from:
                    caught_up = True
                if caught_up == False:
                    i += 1
                    continue

            print("=========")
            print("%s: %i / %i" % (galname, i, n))
            print("=========")
            # gal_dict = do_gal(galname, gal_dict)
            if pick_up_from == 'skip_fit':
                os.system(
                    'python3 project_w3_to_co.py gal_dict_i %s skip_fit' %
                    (galname, ))
            else:
                os.system('python3 project_w3_to_co.py gal_dict_i %s do_fit' %
                          (galname, ))
            # gal_dict = do_gal(galname, dict())
            # with open('/Users/ryan/Dropbox/mac/wise_w3_vs_co/gal_dicts/%s.pk'%(galname,), 'wb') as p:
            #     pickle.dump(gal_dict, p)
            # p.close()
            # del gal_dict
            i += 1

        # with open(fname_l12_lco_dict, 'wb') as p:
        #     pickle.dump(gal_dict, p)

    if sys.argv[1] == 'gal_dict_i':
        galname = sys.argv[2]
        skip_fit = sys.argv[3]
        if skip_fit == 'skip_fit':
            skip_fit = True
        else:
            skip_fit = False

        gal_dict = do_gal(galname, dict(), skip_fitting=skip_fit)
        with open(
                '/Users/ryan/Dropbox/mac/wise_w3_vs_co/gal_dicts/%s.pk' %
            (galname, ), 'wb') as p:
            pickle.dump(gal_dict, p)
        p.close()
        del gal_dict

    if sys.argv[1] == 'combine_gal_dicts':
        fname_out = fname_l12_lco_dict
        if len(sys.argv) > 2:
            skip_fit = sys.argv[2]
            if skip_fit == 'skip_fit':
                fname_out = fname_l12_lco_dict_nofits

        gal_dict = dict()
        for galname in galnames:
            gal = galname_short2long(galname)
            if gal == None:
                continue
            else:
                gal_dict_i = pickle.load(
                    open(
                        '/Users/ryan/Dropbox/mac/wise_w3_vs_co/gal_dicts/%s.pk'
                        % (galname, ), 'rb'))
                gal_dict[gal] = gal_dict_i[gal]
        # Now save the combined dictionary
        with open(fname_out, 'wb') as p:
            pickle.dump(gal_dict, p)

    if sys.argv[1] == 'book':
        version = str(10)
        if len(sys.argv) > 2:
            version = sys.argv[2]
        merger = PdfFileMerger()
        for pdf in glob.glob(
                '/Users/ryan/Dropbox/mac/wise_w3_vs_co/figs_v3/*.pdf'):
            merger.append(open(pdf, 'rb'))
        with open(
                '/Users/ryan/Dropbox/mac/wise_w3_vs_co/booklet_edge_co_wise_v%s.pdf'
                % (version, ), 'wb') as fout:
            merger.write(fout)

# merger = PdfFileMerger()
# for pdf in glob.glob('/Users/ryan/Dropbox/mac/wise_w3_vs_co/figs/*.pdf'):
#     merger.append(open(pdf, 'rb'))
# with open('~/Dropbox/mac/wise_w3_vs_co/booklet_edge_co_wise_v1.pdf',
#           'wb') as fout:
#     merger.write(fout)
#
#
# import itertools
#
# gal_dict_i = gal_dict['IC0944']
# # gal_dict_i = gal_dict['UGC09537']
#
# pixel_size_arcsec = 6.
# pixel_area_pc2 = (pixel_size_arcsec * gal_dict_i['dist_Mpc'] * 1e6 /
#                   206265.)**2
#
# # Convert x and y from log10 to linear units
# x = (10**gal_dict_i['x_fit'])
# x_err = gal_dict_i['x_err_fit'] * x / 0.434
# y = (10**gal_dict_i['y_fit'])
# y_err = gal_dict_i['y_err_fit'] * y / 0.434
#
# # Convert units from luminosity to surface density
# x /= pixel_area_pc2
# x_err /= pixel_area_pc2
# y /= pixel_area_pc2
# y_err /= pixel_area_pc2
#
# # Linear fit on x and y in units of log10(surface density)
# fit_result = run_linmix(np.log10(x), np.log10(y), 0.434 * x_err / x,
#                         0.434 * y_err / y)
#
# n_pix_fit = x.size
# indices = np.arange(n_pix_fit, dtype=int)
#
# x_coadd = np.zeros(n_pix_fit)
# x_coadd_err = np.zeros(n_pix_fit)
# y_coadd = np.zeros(n_pix_fit)
# y_coadd_err = np.zeros(n_pix_fit)
# k = 1
# alphas = [1, 0.1]
# sizes = [5, 1]
# slope_best = fit_result['slope']
# intercept_best = fit_result['intercept']
# xmin, xmax = -1, 2  #7,7.5
# ymin, ymax = -.5, .5
# for j in [1, 4]:  #range(4, 6): #n_pix_fit+1):
#     print(j)
#     index_combos = np.array(list(itertools.combinations(indices, j)))
#     print("Done " + str(j))
#     # if index_combos.shape[0] > 1e3:
#     # index_combos = index_combos[np.linspace(0, index_combos.shape[0]-1, 1000, dtype=int), :]
#     x_sum_each_combo = np.array([np.sum(x[combo])
#                                  for combo in index_combos]) / np.float64(j)
#     x_sum_err_each_combo = np.array(
#         [np.sqrt(np.sum((x_err[combo])**2))
#          for combo in index_combos]) / np.float64(j)
#
#     y_sum_each_combo = np.array([np.sum(y[combo])
#                                  for combo in index_combos]) / np.float64(j)
#     y_sum_err_each_combo = np.array(
#         [np.sqrt(np.sum((y_err[combo])**2))
#          for combo in index_combos]) / np.float64(j)
#
#     x_wt = 1. / x_sum_err_each_combo**2
#     y_wt = 1. / y_sum_err_each_combo**2
#     x_sum_wt = np.sum(x_wt)
#     y_sum_wt = np.sum(y_wt)
#
#     # x_coadd[j-1] = np.sum(x_sum_each_combo*x_wt)/x_sum_wt
#     # x_coadd_err[j-1] = 1./np.sqrt(x_sum_wt)
#     # y_coadd[j-1] = np.sum(y_sum_each_combo*y_wt)/y_sum_wt
#     # y_coadd_err[j-1] = 1./np.sqrt(y_sum_wt)
#     x_coadd[j - 1] = np.average(x_sum_each_combo)
#     x_coadd_err[j - 1] = np.sqrt(np.sum(x_sum_err_each_combo**
#                                         2)) / x_sum_each_combo.size
#     y_coadd[j - 1] = np.average(y_sum_each_combo)
#     y_coadd_err[j - 1] = np.sqrt(np.sum(y_sum_err_each_combo**
#                                         2)) / y_sum_each_combo.size
#
#     pl.subplot(1, 2, k)
#     pl.scatter(np.log10(x_sum_each_combo),
#                np.log10(y_sum_each_combo),
#                s=sizes[k - 1],
#                alpha=alphas[k - 1])
#     if k == 1:
#         pl.errorbar(np.log10(x_sum_each_combo),
#                     np.log10(y_sum_each_combo),
#                     xerr=0.434 * x_sum_err_each_combo / x_sum_each_combo,
#                     yerr=0.434 * y_sum_err_each_combo / y_sum_each_combo,
#                     marker='o',
#                     markersize=4.,
#                     markerfacecolor='none',
#                     markeredgecolor='k',
#                     markeredgewidth=0.5,
#                     capsize=3,
#                     linestyle='none')
#     else:
#         pl.scatter(np.log10(x_sum_each_combo),
#                    np.log10(y_sum_each_combo),
#                    s=sizes[k - 1],
#                    alpha=alphas[k - 1],
#                    c=np.arange(x_sum_each_combo.size))
#
#     pl.errorbar(np.log10(x_coadd[j - 1]),
#                 np.log10(y_coadd[j - 1]),
#                 xerr=0.434 * x_coadd_err[j - 1] / x_coadd[j - 1],
#                 yerr=0.434 * y_coadd_err[j - 1] / y_coadd[j - 1],
#                 marker='o',
#                 markersize=4.,
#                 markerfacecolor='none',
#                 markeredgecolor='k',
#                 markeredgewidth=0.5,
#                 capsize=3,
#                 linestyle='none')
#     pl.ylabel(r"$\log\> L_{12}/A$ ($L_\odot$ pc$^{-2}$)")
#     pl.xlabel(r"$\log\> L_\mathrm{CO}/A$ (K km/s pc$^2$/pc$^2$)")
#     pl.plot(np.linspace(xmin, xmax, 10),
#             np.linspace(xmin, xmax, 10) * slope_best + intercept_best,
#             color='k',
#             linestyle='--')
#     pl.xlim(xmin, xmax)
#     pl.ylim(ymin, ymax)
#     k += 1
#
# A, B = np.log10(x_sum_each_combo), np.log10(y_sum_each_combo)
# combo_stripe_1 = np.arange(x_sum_each_combo.size)[(0.55 < A) & (0.6 > A)]
# combo_stripe_2 = np.arange(x_sum_each_combo.size)[(0.18 < B) & (0.19 > B)]
# combo_stripe = np.intersect1d(combo_stripe_1, combo_stripe_2)
# combo = index_combos[combo_stripe[0]]
#
# map_combo = np.zeros(gal_dict_i['LCO_map'].shape)
# map_combo[np.isnan(gal_dict_i['LCO_map']) == False] = 1
# idcs_map_combo = np.arange(map_combo.size)
# idcs_map_combo = idcs_map_combo[map_combo.flatten() == 1]
#
# map_combo_flat = map_combo.flatten()
# map_combo_flat[idcs_map_combo[combo]] = 2
# map_combo = map_combo_flat.reshape(map_combo.shape)
#
# pl.figure()
# pl.subplot(121)
# pl.imshow(map_combo[::-1, ::-1], origin='lower')
# pl.subplot(122)
# scatterplot(np.log10(x[combo]),
#             np.log10(y[combo]),
#             xlabel=r"$\log\> L_\mathrm{CO}/A$ (K km/s pc$^2$/pc$^2$)",
#             ylabel=r"$\log\> L_{12}/A$ ($L_\odot$ pc$^{-2}$)",
#             xerr=0.434 * x_err[combo] / x[combo],
#             yerr=0.434 * y_err[combo] / y[combo])
# pl.plot(np.linspace(xmin, xmax, 10),
#         np.linspace(xmin, xmax, 10) * slope_best + intercept_best,
#         color='k',
#         linestyle='--')
# k = 2
# pl.scatter(np.log10(x_sum_each_combo),
#            np.log10(y_sum_each_combo),
#            s=sizes[k - 1],
#            alpha=alphas[k - 1],
#            c=np.arange(x_sum_each_combo.size))
# pl.scatter(np.log10(np.average(x[combo])),
#            np.log10(np.average(y[combo])),
#            linestyle='None')
#
# Xfit = X[(np.isinf(X) == False) & (np.isinf(Y) == False)]
# Yfit = Y[(np.isinf(X) == False) & (np.isinf(Y) == False)]
# Xfit_err = X_err[(np.isinf(X) == False) & (np.isinf(Y) == False)]
# Yfit_err = Y_err[(np.isinf(X) == False) & (np.isinf(Y) == False)]
#
# Xf, Yf = Xfit[(np.isnan(Xfit) == False)
#               & (np.isnan(Yfit) == False)], Yfit[(np.isnan(Xfit) == False)
#                                                  & (np.isnan(Yfit) == False)]
# Xf_err, Yf_err = Xfit_err[(np.isnan(Xfit) == False) & (
#     np.isnan(Yfit) == False)], Yfit_err[(np.isnan(Xfit) == False)
#                                         & (np.isnan(Yfit) == False)]
#
# # pl.scatter( np.log10(x_coadd)+np.log10(pixel_area_pc2), np.log10(y_coadd)+np.log10(pixel_area_pc2), c=indices )
# pl.errorbar(np.log10(x_coadd),
#             np.log10(y_coadd),
#             xerr=0.434 * x_coadd_err / x_coadd,
#             yerr=0.434 * y_coadd_err / y_coadd,
#             marker='o',
#             markersize=4.,
#             markerfacecolor='none',
#             markeredgecolor='k',
#             markeredgewidth=0.5,
#             capsize=3,
#             linestyle='none')
# slope_best = gal_dict_i['linmix_fit']['slope']
# intercept_best = gal_dict_i['linmix_fit']['intercept']
# xmin, xmax = 0, 1  #7,7.5
# pl.plot(np.linspace(xmin, xmax, 10),
#         np.linspace(xmin, xmax, 10) * slope_best + intercept_best,
#         color='k',
#         linestyle='--')
# pl.ylabel(r"$\log\> L_{12}/A$ ($L_\odot$ pc$^{-2}$)")
# pl.xlabel(r"$\log\> L_\mathrm{CO}/A$ (K km/s pc$^2$/pc$^2$)")
#
# # fname_fit_surf_dens_dict = '/Users/ryan/Dropbox/mac/wise_w3_vs_co/fit_surf_dens.pk' # Without Aniano kernels, until that is fixed, and with more conservative s/n cuts in mom0 map making
# # with open(fname_fit_surf_dens_dict, 'wb') as p:
# #     pickle.dump(fit_results_surfdens, p)
# #
# # fname_fit_surf_dens_dict = '/Users/ryan/Dropbox/mac/wise_w3_vs_co/fit_surf_dens.pk' # Without Aniano kernels, until that is fixed, and with more conservative s/n cuts in mom0 map making
# # fit_results_surfdens = pickle.load(open(fname_fit_surf_dens_dict, 'rb'))
# #
#
# # When you obtain the mag from raw flux(1 > 18 mag for W3), the code of computing luminosity is:
# # cconst = 2.9979246  ; * 10^10 cm s^-1
# # wband = [3.3526, 4.6028, 11.5608, 22.0883]  ;micro meter
# # f0nu = [309.54, 171.787, 31.674, 8.363] ;[306.681, 170.663, 29.0448, 8.2839]  ;Jy
# # bw = [1.7506, 1.4653, 1.1327, 2.4961]  ; bandwidth
# # c1 = -19. ; 1Jy = 10^(-19) erg s^-1 m^-2 Hz^1
# # c2 = 13.  ; bw* 10^13. bandwidth, in units of Hz
# # c3 = alog10(4. * !pi) + 2.* alog10(3.08567) + 44.  ; Mpc=3.08*10^22 m
# # lglumsun = alog10(3.846)+33  ; solar lum. ergs s^-1
# # lgfconst = alog10(f0nu) + alog10(bw) + c1 + c2
# # lumd = lumdist(z, H0=70.2, Omega_M = 0.275, Lambda0 = 0.725)
# # lgf0 = lgfconst[band-1]
# # lglum = lgf0 - 0.4*mag + 2.*alog10(lumd) + c3 - lglumsun
# #
# # thisf0nu = f0nu[band-1]
# # flux = thisf0nu * 10.^(- 0.4 * mag)   ; these two lines for flux in Jansky.
#
# pl.figure()
# bb = 0
# aa = 0
# Lsun = 3.839e33  # Lsun in erg/s
# for g in gal_dict.keys():
#     if g in names:
#         if bam[names == g] == 'B':
#             if bb == 0:
#                 pl.scatter(np.log10(
#                     np.sum(10**(gal_dict[g]['log10_Lwise'] - np.log10(Lsun))))
#                            + np.log10(2.31e-11),
#                            np.log10(np.sum(10**(gal_dict[g]['log10_LCO']))),
#                            s=4.,
#                            c='b',
#                            marker='s',
#                            facecolor='none')
#                 pl.scatter(gal_dict[g]['log10_Lwise'] - np.log10(Lsun) +
#                            np.log10(2.31e-11),
#                            gal_dict[g]['log10_LCO'],
#                            s=10.,
#                            c='b',
#                            marker='+',
#                            label='Barred')
#                 bb = 1
#             else:
#                 pl.scatter(np.log10(
#                     np.sum(10**(gal_dict[g]['log10_Lwise'] - np.log10(Lsun))))
#                            + np.log10(2.31e-11),
#                            np.log10(np.sum(10**(gal_dict[g]['log10_LCO']))),
#                            s=4.,
#                            c='b',
#                            marker='s',
#                            facecolor='none')
#                 pl.scatter(gal_dict[g]['log10_Lwise'] - np.log10(Lsun) +
#                            np.log10(2.31e-11),
#                            gal_dict[g]['log10_LCO'],
#                            s=10.,
#                            c='b',
#                            marker='+')
#
#         if bam[names == g] == 'A':
#             if aa == 0:
#                 pl.scatter(np.log10(
#                     np.sum(10**(gal_dict[g]['log10_Lwise'] - np.log10(Lsun))))
#                            + np.log10(2.31e-11),
#                            np.log10(np.sum(10**(gal_dict[g]['log10_LCO']))),
#                            s=4.,
#                            c='g',
#                            marker='s',
#                            facecolor='none')
#                 pl.scatter(gal_dict[g]['log10_Lwise'] - np.log10(Lsun) +
#                            np.log10(2.31e-11),
#                            gal_dict[g]['log10_LCO'],
#                            s=10.,
#                            c='g',
#                            marker='+',
#                            label='Unbarred')
#                 aa = 1
#             else:
#                 pl.scatter(np.log10(
#                     np.sum(10**(gal_dict[g]['log10_Lwise'] - np.log10(Lsun))))
#                            + np.log10(2.31e-11),
#                            np.log10(np.sum(10**(gal_dict[g]['log10_LCO']))),
#                            s=4.,
#                            c='g',
#                            marker='s',
#                            facecolor='none')
#                 pl.scatter(gal_dict[g]['log10_Lwise'] - np.log10(Lsun) +
#                            np.log10(2.31e-11),
#                            gal_dict[g]['log10_LCO'],
#                            s=10.,
#                            c='g',
#                            marker='+')
#     # pl.scatter( np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise']+np.log10(1.375/6.)-np.log10(Lsun)))), np.log10(np.sum( 10**(gal_dict[g]['log10_LCO']))), s=4. , c='b', marker='s')
#     # pl.scatter(gal_dict[g]['log10_Lwise']+np.log10(1.375/6.)-np.log10(Lsun), gal_dict[g]['log10_LCO'], s=4. , c='k')
#
# pl.figure()
#
# WISE_COR = 1.  #2.2636905e-12 #2.31e-11 # 2.31e-11#2.31e-11
# logL12_all = np.zeros(len(gal_dict.keys()))
# logLCO_all = np.zeros(len(gal_dict.keys()))
# logL12_all_err = np.zeros(len(gal_dict.keys()))
# logLCO_all_err = np.zeros(len(gal_dict.keys()))
# logL12_all_pix = []
# logLCO_all_pix = []
# logL12_all_pix_err = []
# logLCO_all_pix_err = []
# fit_result = np.zeros([len(gal_dict.keys()), 2])
# i = 0
# pp = 0
# lpix = 'Pixels'
# lglob = 'Totals'
# lline = 'Pixel-pixel fit'
#
# galaxies_to_plot = list(gal_dict.keys())
# # galaxies_to_plot.remove('UGC05359')
# # galaxies_to_plot.remove('UGC03253')
#
# for g in galaxies_to_plot:
#     ltemp_pix = ''
#     ltemp_glob = ''
#     ltemp_line = ''
#     if pp == 0:
#         ltemp_pix = lpix
#         ltemp_glob = lglob
#         ltemp_line = lline
#         pp = 1
#     X = gal_dict[g]['log10_LCO']
#     Xerr = gal_dict[g]['log10_LCO_err']
#     # Y = gal_dict[g]['log10_Lwise'] -np.log10(Lsun)+np.log10(WISE_COR)
#     Y = gal_dict[g][
#         'log10_Lwise']  # -np.log10(Lsun) + np.log10(2.2636905e-12)#- np.log10(2.2636905e-12)
#     Yerr = gal_dict[g]['log10_Lwise_err']
#
#     Ysum = gal_dict[g][
#         'log10_Lwise_tot']  # np.log10(np.sum( 10**(gal_dict[g]['log10_LCO'])))
#     Ysum_err = gal_dict[g]['log10_Lwise_tot_err']
#     # Ysum = np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise']-np.log10(Lsun)+np.log10(WISE_COR))))
#     Xsum = gal_dict[g][
#         'log10_LCO_tot']  # np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise'])))
#     Xsum_err = gal_dict[g]['log10_LCO_tot_err']
#     if g != 'UGC09476':
#         # pl.scatter( Xsum, Ysum, s=10. , edgecolor='b', marker='s', facecolor='none', label=ltemp_glob)
#         # pl.scatter( X, Y, s=5. , alpha=0.5, edgecolor='k', facecolor='none', label=ltemp_pix)
#         if Xsum > 4 and Xsum < 12:
#             pl.errorbar(Xsum,
#                         Ysum,
#                         xerr=Xsum_err,
#                         yerr=Ysum_err,
#                         markersize=10.,
#                         linestyle='none',
#                         alpha=0.5,
#                         capsize=5,
#                         marker='s',
#                         ecolor='b',
#                         markeredgecolor='b',
#                         markerfacecolor='none',
#                         label=ltemp_glob)
#             pl.errorbar(X,
#                         Y,
#                         xerr=Xerr,
#                         yerr=Yerr,
#                         markersize=5.,
#                         linestyle='none',
#                         alpha=0.2,
#                         capsize=5,
#                         markeredgecolor='k',
#                         ecolor='k',
#                         markerfacecolor='none',
#                         label=ltemp_pix)
#
#     # if X.size > 1:
#     #     # popt, pcov = curve_fit(fun, X[(np.isinf(X)==False)&(np.isinf(Y)==False)], Y[(np.isinf(X)==False)&(np.isinf(Y)==False)], p0=[1,0])
#     #     pfit = np.poly1d(np.polyfit(X, Y, 1))
#     #     # print(popt)
#     #     # fit_result[i] = popt
#     #     # pl.plot( np.linspace(5,10,10), np.linspace(5,10,10)*popt[0] + popt[1], linewidth=0.5, alpha=0.5, c='r', label=ltemp_line)
#     #     pl.plot( np.linspace(5,10,10), pfit(np.linspace(5,10,10)), linewidth=0.5, alpha=0.5, c='r', label=ltemp_line)
#     logL12_all[
#         i] = Ysum  #np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise']-np.log10(Lsun)+np.log10(WISE_COR))))
#     logLCO_all[i] = Xsum  # np.log10(np.sum( 10**(gal_dict[g]['log10_LCO'])))
#
#     logL12_all_err[
#         i] = Ysum_err  #np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise']-np.log10(Lsun)+np.log10(WISE_COR))))
#     logLCO_all_err[
#         i] = Xsum_err  # np.log10(np.sum( 10**(gal_dict[g]['log10_LCO'])))
#
#     logL12_all_pix += Y[(np.isinf(X) == False) & (np.isinf(Y) == False)].tolist(
#     )  #np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise']-np.log10(Lsun)+np.log10(WISE_COR))))
#     logLCO_all_pix += X[(np.isinf(X) == False)
#                         & (np.isinf(Y) == False)].tolist(
#                         )  # np.log10(np.sum( 10**(gal_dict[g]['log10_LCO'])))
#
#     logL12_all_pix_err += Yerr[(np.isinf(X) == False)
#                                & (np.isinf(Y) == False)].tolist()
#     logLCO_all_pix_err += Xerr[(np.isinf(X) == False)
#                                & (np.isinf(Y) == False)].tolist()
#
#     i += 1
#
# pl.xlabel(r"$\log \> L(\mathrm{CO(1-0)})$ [K km/s pc$^2$]", fontsize=15)
# pl.ylabel(r"$\log \> L(12\mu\mathrm{m})$ [L$_\odot$]", fontsize=15)
# # pl.plot(np.linspace(5,10,10), (np.linspace(5,10,10)+0.06)/0.97, label='Gao+ (in prep.)')
# pl.plot(np.linspace(5, 10, 10),
#         np.linspace(5, 10, 10) * 1.13 - .967,
#         label='Jiang+15')
# pl.plot(np.linspace(5, 10, 10),
#         np.linspace(5, 10, 10) * 1.13 + (np.log10(4.35) - 1.49) / .88,
#         label='Jiang+15')
#
# ks_yang = [.98, 1.03]
# bs_yang = [-.14, -.64]
# for i in range(0, 2):
#     pl.plot(np.linspace(5, 10, 10),
#             (np.linspace(5, 10, 10) - bs_yang[i]) / ks_yang[i])
#
# l12_yang = np.array([8.711581, 8.7366333, 7.9821650, 8.5337734, 9.4445164])
# names_yang = ['NGC2487', 'NGC4644', 'NGC5485', 'NGC5520', 'UGC05111']
#
# i = 0
# for name in names_yang:
#     print("Name: %s" % (name, ))
#     l12_yg_i = l12_yang[i]
#     l12_rc_i = np.log10(np.sum(gal_dict[name]['Lwise_map']))
#     l12_rc_co_i = gal_dict[name]['log10_Lwise_tot']
#     lco_rc_co_i = gal_dict[name]['log10_LCO_tot']
#     pl.scatter(lco_rc_co_i, l12_yg_i, label=name)
#
#     i += 1
#
# l12_yang = np.array([8.711581, 8.7366333, 7.9821650, 8.5337734, 9.4445164])
# names_yang = ['NGC2487', 'NGC4644', 'NGC5485', 'NGC5520', 'UGC05111']
#
# i = 0
# for name in names_yang:
#     print("Name: %s" % (name, ))
#     l12_yg_i = l12_yang[i]
#     l12_rc_i = np.log10(np.sum(gal_dict[name]['Lwise_map']))
#     l12_rc_co_i = gal_dict[name]['log10_Lwise_tot']
#
#     idx = np.where(rfpars['Name'] == name)[0][0]
#     pa = rfpars['rfPA'][idx]
#     inc = rfpars['rfInc'][idx]
#
#     l12_rc_i -= np.log10(np.cos(np.radians(inc)))
#     l12_rc_co_i -= np.log10(np.cos(np.radians(inc)))
#
#     print("    Yang: %3.3f;   Me: %3.3f, %3.3f " %
#           (l12_yang[i], gal_dict[name]['log10_Lwise_tot'],
#            np.log10(np.sum(gal_dict[name]['Lwise_map']))))
#     print("    Ratio:  %3.3f, %3.3f " %
#           ((l12_yang[i] - gal_dict[name]['log10_Lwise_tot']),
#            l12_yang[i] - np.log10(np.sum(gal_dict[name]['Lwise_map']))))
#
#     mdat, mhdr = fits.getdata(
#         '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
#         % (name, ),
#         header=True)
#     w = wcs.WCS(mhdr)
#     idx = np.where(leda['Name'] == name)[0][0]
#     ractr = leda['ledaRA'].quantity[idx].to(u.deg).value
#     dcctr = leda['ledaDE'].quantity[idx].to(u.deg).value
#     x, y = w.wcs_world2pix(ractr, dcctr, 0)
#     wcs_out = wcs.WCS(naxis=2)
#     wcs_out.wcs.crpix = [x, y]
#     # what is the galactic coordinate of that pixel.
#     wcs_out.wcs.crval = [ractr, dcctr]
#     # what is the pixel scale in lon, lat.
#     wcs_out.wcs.cdelt = np.array([mhdr['CDELT1'], mhdr['CDELT2']])
#     wcs_out.wcs.ctype = ['RA---SIN', 'DEC--SIN']
#
#     pl.subplot(5, 2, 2 * i + 1, projection=wcs_out)
#     pl.imshow(gal_dict[name]['Lwise_map'])
#
#     pl.subplot(5, 2, 2 * i + 2)
#     ax = pl.gca()
#     ax.axis('off')
#     ymax = 1.
#     pl.text(0, ymax, name, fontsize=15)
#     pl.text(0,
#             ymax - 0.15,
#             r"$\log_{10}L_{12,\mathrm{YG}}/L_\odot=%3.3f$" % (l12_yg_i, ),
#             fontsize=12)
#     pl.text(0,
#             ymax - 2 * 0.15,
#             r"$\log_{10}L_{12,\mathrm{RC}}/L_\odot=%3.3f$ (all pix)" %
#             (l12_rc_i, ),
#             fontsize=12)
#     pl.text(0,
#             ymax - 3 * 0.15,
#             r"$\log_{10}L_{12,\mathrm{RC}}/L_\odot=%3.3f$ (CO pix)" %
#             (l12_rc_co_i, ),
#             fontsize=12)
#     pl.text(0,
#             ymax - 4 * 0.15,
#             r"YG-RC$=%3.3f$ (all pix)" % (l12_yg_i - l12_rc_i, ),
#             fontsize=12)
#     pl.text(0,
#             ymax - 5 * 0.15,
#             r"YG-RC$=%3.3f$ (CO pix)" % (l12_yg_i - l12_rc_co_i, ),
#             fontsize=12)
#
#     i += 1
#
# import astropy.constants.si as _si
#
# c_si = _si.c.value  # m/s
# c2overkb = c_si**2 / _si.k_B.value  # c^2 / Boltzmann in SI units
# c3overkb = c_si * c2overkb  # c^3 / Boltzmann in SI units
# lsun_kkmspc2_const = _si.L_sun.value * c3overkb / (4.0 * pi * u.pc.to(u.m)**2)
# frest = 115.3 * 1e9
# kkmspc2_to_lsun = 2e3 * frest**3 / lsun_kkmspc2_const
# lsun_to_kkmspc2 = 5e-4 * lsun_kkmspc2_const * frest**(-3)
#
# # Pixel-pixel fit
# good_pixel_indices = np.arange(len(logL12_all_pix))
# good_pixel_indices = good_pixel_indices[
#     np.array(logL12_all_pix)[good_pixel_indices] > 5.]
# good_pixel_indices = good_pixel_indices[
#     np.array(logLCO_all_pix)[good_pixel_indices] > 5.]
#
# fit_pixels = run_linmix(
#     np.array(logLCO_all_pix)[good_pixel_indices],
#     np.array(logL12_all_pix)[good_pixel_indices],
#     np.array(logLCO_all_pix_err)[good_pixel_indices],
#     np.array(logL12_all_pix_err)[good_pixel_indices])
#
# # Global fit
# good_pixel_indices = np.arange(len(logL12_all))
# good_pixel_indices = good_pixel_indices[
#     np.array(logL12_all)[good_pixel_indices] > 5.]
# good_pixel_indices = good_pixel_indices[
#     np.array(logLCO_all)[good_pixel_indices] > 5.]
# fit_global = run_linmix(
#     np.array(logLCO_all)[good_pixel_indices],
#     np.array(logL12_all)[good_pixel_indices],
#     np.array(logLCO_all_err)[good_pixel_indices],
#     np.array(logL12_all_err)[good_pixel_indices])
#
# print("Pixel-pixel fit: (%2.2f +/- %2.2f) x + (%2.2f +/- %2.2f)" %
#       (fit_pixels['slope'], fit_pixels['slope_err'], fit_pixels['intercept'],
#        fit_pixels['intercept_err']))
# print("Global fit: (%2.2f +/- %2.2f) x + (%2.2f +/- %2.2f)" %
#       (fit_global['slope'], fit_global['slope_err'], fit_global['intercept'],
#        fit_global['intercept_err']))
#
# pl.plot(np.linspace(5, 10, 10),
#         np.linspace(5, 10, 10) * fit_global['slope'] + fit_global['intercept'],
#         label='Global fit')
# pl.plot(np.linspace(5, 10, 10),
#         np.linspace(5, 10, 10) * fit_pixels['slope'] + fit_pixels['intercept'],
#         label='Pixel fit')
#
# np.sqrt(
#     np.average(
#         (np.array(logL12_all)[good_pixel_indices] -
#          fit_global['slope'] * np.array(logLCO_all)[good_pixel_indices] -
#          fit_global['intercept'])**2))
#
# np.array(galaxies_to_plot)[good_pixel_indices]
# # 82 galaxies
# # 30 barred
# # Mostly late-type spirals
# # 14 interacting
# # 72-128 Mpc distance (2 kpc pixel size on avg. at 6 arcsec)
# k = 0
# i_leda = []
# for g in np.array(galaxies_to_plot)[good_pixel_indices]:
#     if g in np.array(leda['Name']):
#         print(g)
#         i_leda.append(np.where(np.array(leda['Name']) == g)[0][0])
#         k += 1
# i_leda = np.array(i_leda)
#
# from collections import Counter
#
# # a = ['a', 'a', 'a', 'a', 'b', 'b', 'c', 'c', 'c', 'd', 'e', 'e', 'e', 'e', 'e']
# # morph_counts = Counter(list(leda['ledaMorph'][i_leda]))
# # df = pd.DataFrame.from_dict(morph_counts, orient='index') #, columns=np.array(['S0', 'S0-a', 'Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd', 'Sd']))#orient='index')
# # df.plot(kind='bar')
#
# morph_counts = Counter(np.array(leda['ledaMultiple'][i_leda]))
# df = pd.DataFrame.from_dict(
#     morph_counts, orient='index'
# )  #, columns=np.array(['S0', 'S0-a', 'Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd', 'Sd']))#orient='index')
# df.plot(kind='bar')
#
# morph_counts = Counter(np.array(leda['ledaDistMpc'][i_leda]))
# df = pd.DataFrame.from_dict(
#     morph_counts, orient='index'
# )  #, columns=np.array(['S0', 'S0-a', 'Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd', 'Sd']))#orient='index')
# df.plot(kind='bar')
#
# np.average(6. * np.array(leda['ledaDistMpc'][i_leda]) * 1e3 / 206265.)
#
# redshifts = np.loadtxt(
#     '/Users/ryan/Downloads/stz349_supplemental_files/TableA1.csv',
#     usecols=5,
#     delimiter=',',
#     skiprows=1)
# names = np.loadtxt(
#     '/Users/ryan/Downloads/stz349_supplemental_files/TableA1.csv',
#     usecols=0,
#     delimiter=',',
#     skiprows=1,
#     dtype=str)
# bam = np.loadtxt('/Users/ryan/Downloads/stz349_supplemental_files/TableA1.csv',
#                  usecols=1,
#                  delimiter=',',
#                  skiprows=1,
#                  dtype=str)
#
# table_a2 = pd.read_csv(
#     '/Users/ryan/Downloads/stz349_supplemental_files/tableA2.csv')
# galnames_a2 = np.array(table_a2['name'])
# delta_d4000 = np.array(table_a2['Delta Dn4000'])
# delta_ew_hdelta = np.array(table_a2['Delta EW Hdelta_A'])
# delta_ew_halpha = np.array(table_a2['Delta log EW Halpha'])
#
# k = 0
# i_table_a2 = []
# i_plot = []
# for g in np.array(galaxies_to_plot)[good_pixel_indices]:
#     if g in galnames_a2:
#         print(g)
#         i_table_a2.append(np.where(galnames_a2 == g)[0][0])
#         i_plot.append(k)
#         k += 1
# i_table_a2 = np.array(i_table_a2)
# i_plot = np.array(i_plot)
#
# x_plot, y_plot, xerr_plot, yerr_plot = np.array(
#     logLCO_all)[good_pixel_indices][i_plot], np.array(
#         logL12_all)[good_pixel_indices][i_plot], np.array(
#             logLCO_all_err)[good_pixel_indices][i_plot], np.array(
#                 logL12_all_err)[good_pixel_indices][i_plot]
#
# pl.figure()
# pl.subplot(131)
# pl.scatter(x_plot, y_plot, c=delta_ew_halpha[i_table_a2])
# pl.xlabel(r"$\log \> L(\mathrm{CO(1-0)})$ [K km/s pc$^2$]", fontsize=15)
# pl.ylabel(r"$\log \> L(12\mu\mathrm{m})$ [L$_\odot$]", fontsize=15)
# pl.plot(np.linspace(5, 10, 10), (np.linspace(5, 10, 10) + 0.06) / 0.97,
#         label='Gao+ (in prep.)')
# pl.colorbar(label=r'$\Delta \log$ EW(H$\alpha$)')
# pl.subplot(132)
# pl.scatter(x_plot, y_plot, c=delta_ew_hdelta[i_table_a2])
# pl.xlabel(r"$\log \> L(\mathrm{CO(1-0)})$ [K km/s pc$^2$]", fontsize=15)
# pl.ylabel(r"$\log \> L(12\mu\mathrm{m})$ [L$_\odot$]", fontsize=15)
# pl.plot(np.linspace(5, 10, 10), (np.linspace(5, 10, 10) + 0.06) / 0.97,
#         label='Gao+ (in prep.)')
# pl.colorbar(label=r'$\Delta$ EW(H$\delta_A$) (\\A)')
# pl.subplot(133)
# pl.scatter(x_plot, y_plot, c=delta_d4000[i_table_a2])
# pl.xlabel(r"$\log \> L(\mathrm{CO(1-0)})$ [K km/s pc$^2$]", fontsize=15)
# pl.ylabel(r"$\log \> L(12\mu\mathrm{m})$ [L$_\odot$]", fontsize=15)
# pl.plot(np.linspace(5, 10, 10), (np.linspace(5, 10, 10) + 0.06) / 0.97,
#         label='Gao+ (in prep.)')
# pl.colorbar(label=r'$\Delta$D$_n$(4000)')
#
# import sklearn
# from sklearn.decomposition import PCA
# X = np.vstack([
#     x_plot, y_plot, delta_ew_halpha[i_table_a2], delta_ew_hdelta[i_table_a2],
#     delta_d4000[i_table_a2]
# ]).T
# X -= np.average(X, axis=0)
# # X /= np.sum(X, axis=0)
# pca = PCA(n_components=5)
# pca.fit(X)
# print(pca.explained_variance_ratio_)
# # pca.components_
# pc5 = pca.components_[4]
#
# pl.scatter(
#     -1. * (X[:, 1] * pc5[1] + X[:, 2] * pc5[2] + X[:, 3] * pc5[3] +
#            X[:, 4] * pc5[4]) / pc5[0], X[:, 0] + np.average(X[:, 0]))
# pl.xlabel(
#     r"$\mathrm{PC5}(\log L_\mathrm{12}, \Delta D(4000), \Delta \mathrm{EW\> H\alpha}, \Delta \mathrm{EW\>H\delta_A} )$"
# )
# pl.ylabel(r"$\log L_\mathrm{CO}$ [K km/s pc$^2$]")
#
# # edge_califa = Table.read('/usr/local/lib/python3.7/site-packages/edge_pydb/dat_glob/external/edge_califa.csv',format='ascii.ecsv')
#
# x_plot, y_plot, xerr_plot, yerr_plot = np.array(
#     logLCO_all)[good_pixel_indices], np.array(
#         logL12_all)[good_pixel_indices], np.array(
#             logLCO_all_err)[good_pixel_indices], np.array(
#                 logL12_all_err)[good_pixel_indices]
#
# print("Global fit: (%2.2f +/- %2.2f) x + (%2.2f +/- %2.2f)" %
#       (fit_global['slope'], fit_global['slope_err'], fit_global['intercept'],
#        fit_global['intercept_err']))
#
# dist_xy = np.abs(fit_global['slope'] * x_plot - y_plot +
#                  fit_global['intercept']) / np.sqrt(1 + fit_global['slope']**2)
#
# k = 0
# i_edge_califa = []
# for g in np.array(galaxies_to_plot)[good_pixel_indices]:
#     if g in list(edge_califa['Name']):
#         print(g)
#         i_edge_califa.append(
#             np.where(np.array(edge_califa['Name']) == g)[0][0])
#         k += 1
# i_edge_califa = np.array(i_edge_califa)
#
# mstar_edge = np.array(edge_califa['caMass'][i_edge_califa])
# sfr_edge = np.array(edge_califa['caSFR'][i_edge_califa])
# age_edge = np.array(edge_califa['caAge'][i_edge_califa])
#
# # X = np.vstack([x_plot, y_plot, mstar_edge, sfr_edge, age_edge]).T
# X = np.vstack([dist_xy, mstar_edge, sfr_edge, age_edge]).T
# X = np.vstack([X[:37], X[38:]])
# X -= np.average(X, axis=0)
# # X /= np.sum(X, axis=0)
# pca = PCA(n_components=4)
# pca.fit(X)
# print(pca.explained_variance_ratio_)
# # pca.components_
# pc5 = pca.components_[3]
#
# # pl.scatter( -1.*(X[:,1]*pc5[1] + X[:,2]*pc5[2] + X[:,3]*pc5[3] +X[:,4]*pc5[4])/pc5[0] , X[:,0]+np.average(X[:,0]) )
# pl.scatter(-1. * (X[:, 1] * pc5[1] + X[:, 2] * pc5[2] + X[:, 3] * pc5[3]) /
#            pc5[0], X[:, 0])  #+np.average(X[:,0]) )
# pl.xlabel(
#     r"$\mathrm{PC5}(\log L_\mathrm{12}, \log M_*, \log\mathrm{SFR}, \log \mathrm{Age} )$"
# )
# pl.ylabel(r"$\log L_\mathrm{CO}$ [K km/s pc$^2$]")
#
# # pl.scatter( X
#
# Xsum = logLCO_all
# Ysum = logL12_all
# # pfit_sum = np.poly1d(np.polyfit(Xsum[(np.isinf(Xsum)==False)&(np.isinf(Ysum)==False)], Ysum[(np.isinf(Xsum)==False)&(np.isinf(Ysum)==False)], 1))
# cfit_sum = curve_fit(
#     fun, Xsum[(np.isinf(Xsum) == False) & (np.isinf(Ysum) == False)],
#     Ysum[(np.isinf(Xsum) == False) & (np.isinf(Ysum) == False)])
# pfit_sum = np.poly1d(cfit_sum[0])
#
# Xpix = np.array(logLCO_all_pix)
# Ypix = np.array(logL12_all_pix)
# Xpix = Xpix[Ypix > 5]
# Ypix = Ypix[Ypix > 5]
#
# # pfit_pix = np.poly1d(np.polyfit(Xpix[(np.isinf(Xpix)==False)&(np.isinf(Ypix)==False)], Ypix[(np.isinf(Xpix)==False)&(np.isinf(Ypix)==False)], 1))
# # pfit_pix = np.poly1d(np.polyfit(Xpix[(np.isinf(Xpix)==False)&(np.isinf(Ypix)==False)], Ypix[(np.isinf(Xpix)==False)&(np.isinf(Ypix)==False)], 1))
# cfit_pix = curve_fit(
#     fun, Xpix[(np.isinf(Xpix) == False) & (np.isinf(Ypix) == False)],
#     Ypix[(np.isinf(Xpix) == False) & (np.isinf(Ypix) == False)])
# pfit_pix = np.poly1d(cfit_pix[0])
#
# pl.plot(np.linspace(5, 10, 10),
#         pfit_sum(np.linspace(5, 10, 10)),
#         label='Global fit',
#         c='g')
# pl.plot(np.linspace(5, 10, 10),
#         pfit_pix(np.linspace(5, 10, 10)),
#         label='Pixel-pixel fit',
#         c='r')
#
# pl.legend(loc='best')
#
# pl.xlim(5, 10)
# pl.ylim(5, 10)
#
# # pl.savefig('/Users/ryan/Dropbox/mac/wise_w3_vs_co/l12_lco.pdf')
# pl.savefig('/Users/ryan/Dropbox/mac/wise_w3_vs_co/l12_lco_23Aug19.pdf')
#
#
# galaxy_dict = dict( galname = list(gal_dict.keys()),\
#     logL12_total = logL12_all,\
#     logLCO_total = logLCO_all)
#
# df = pd.DataFrame(galaxy_dict)
# cols = ['galname', 'logL12_total', 'logLCO_total']
#
# df.to_csv("/Users/ryan/Dropbox/mac/tsinghua/edge_all_logL12_logLCO.csv",
#           sep=",",
#           float_format='%.9g',
#           header=True,
#           index=False,
#           columns=cols)
#
#
# def coordinate_meshgrid(wcs, naxis1, naxis2, center, pixel_size_kpc):
#     x = np.arange(naxis1)
#     y = np.arange(naxis2)
#     X, Y = np.meshgrid(x, y)
#     ra, dec = wcs.wcs_pix2world(X.flatten(), Y.flatten(), 0)
#     coords = np.dstack((ra, dec))
#     pixel_size_arcsec = 6.
#     i = 0
#     res = []
#     for c1 in coords:
#         c1 = SkyCoord(ra[i] * u.deg, dec[i] * u.deg, frame=ICRS)
#         sep = c1.separation(
#             center).arcsecond / pixel_size_arcsec * pixel_size_kpc
#         print(sep, ra[i], dec[i])
#         if sep <= 10.:
#             pixcoords = w.wcs_world2pix(ra[i], dec[i], 0)
#             pixcoords = [int(j) for j in pixcoords]
#             print(pixcoords)
#             res.append(pixcoords)
#         i += 1
#     return res
#
#
# pl.figure()
#
# WISE_COR = 1.  #2.2636905e-12 #2.31e-11 # 2.31e-11#2.31e-11
# logL12_all = np.zeros(len(gal_dict.keys()))
# logLCO_all = np.zeros(len(gal_dict.keys()))
# logL12_all_err = np.zeros(len(gal_dict.keys()))
# logLCO_all_err = np.zeros(len(gal_dict.keys()))
#
# logL12_all_pix = []
# logLCO_all_pix = []
# logL12_all_pix_err = []
# logLCO_all_pix_err = []
#
# fit_result = np.zeros([len(gal_dict.keys()), 2])
# i = 0
# pp = 0
# lpix = 'Pixels'
# lglob = 'Totals'
# lline = 'Pixel-pixel fit'
#
# galaxies_to_plot = list(gal_dict.keys())
# # galaxies_to_plot.remove('UGC05359')
# # galaxies_to_plot.remove('UGC03253')
#
# bb = 0
# aa = 0
# barred_unbarred_pix = []
# barred_unbarred_global = []
#
# for g in galaxies_to_plot:
#     ltemp_pix = ''
#     ltemp_glob = ''
#     ltemp_line = ''
#     if pp == 0:
#         ltemp_pix = lpix
#         ltemp_glob = lglob
#         ltemp_line = lline
#         pp = 1
#
#     mdat, mhdr = fits.getdata(
#         '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
#         % (g, ),
#         header=True)
#     w = wcs.WCS(mhdr)
#     idx = np.where(leda['Name'] == g)[0][0]
#     ractr = leda['ledaRA'].quantity[idx].to(u.deg).value
#     dcctr = leda['ledaDE'].quantity[idx].to(u.deg).value
#     r25 = leda['ledaD25'].quantity[idx].to(u.arcsec).value / 2.
#     dist_pc = leda['ledaDistMpc'].quantity[idx].to(u.Mpc).value * 1e6
#     dist_m = dist_pc * 3.086e16
#
#     x, y = w.wcs_world2pix(ractr, dcctr, 0)
#
#     lco_map = gal_dict[g]['LCO_map']
#     pixsize_kpc = 6. * dist_pc * 1.e-3 / 206265.
#
#     # i_mask = []
#     # if pixsize_kpc >= 1.:
#     #     # look at central pixel
#     #     lco_central = lco_map[x,y]
#     #     if np.isnan(lco_central):
#     #         print("Central pixel isn't good")
#     #     else:
#
#     # List of indices
#     center = SkyCoord(ractr * u.deg, dcctr * u.deg, frame=ICRS)
#     # naxis1 = mdat.shape[0]
#     # naxis2 = mdat.shape[1]
#     # i_mask = coordinate_meshgrid(w, naxis1, naxis2, center, pixsize_kpc)
#     aperture_temp = EllipticalAperture((x, y),
#                                        a=1. / pixsize_kpc,
#                                        b=1. / pixsize_kpc,
#                                        theta=0.)
#     m_gal_tmp = aperture_temp.to_mask(method='center')[0]
#     mask_central = m_gal_tmp.to_image(shape=gal_dict[g]['LCO_map'].shape)
#     pre_array = gal_dict[g]['LCO_map'][gal_dict[g]['LCO_map'] > 0]
#     pre_mask = mask_central[gal_dict[g]['LCO_map'] > 0]
#     i_mask = np.where(pre_mask == 1)[0]
#
#     X = gal_dict[g]['log10_LCO'][i_mask]
#     Xerr = gal_dict[g]['log10_LCO_err'][i_mask]
#     # Y = gal_dict[g]['log10_Lwise'] -np.log10(Lsun)+np.log10(WISE_COR)
#     Y = gal_dict[g]['log10_Lwise'][
#         i_mask]  # -np.log10(Lsun) + np.log10(2.2636905e-12)#- np.log10(2.2636905e-12)
#     Yerr = gal_dict[g]['log10_Lwise_err'][i_mask]
#
#     Ysum = gal_dict[g][
#         'log10_Lwise_tot']  # np.log10(np.sum( 10**(gal_dict[g]['log10_LCO'])))
#     Ysum_err = gal_dict[g]['log10_Lwise_tot_err']
#     # Ysum = np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise']-np.log10(Lsun)+np.log10(WISE_COR))))
#     Xsum = gal_dict[g][
#         'log10_LCO_tot']  # np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise'])))
#     Xsum_err = gal_dict[g]['log10_LCO_tot_err']
#     if g in names:
#         if bam[names == g] == 'B':
#             if bb == 0:
#                 if g != 'UGC09476':
#                     if Xsum > 4 and Xsum < 12:
#                         pl.errorbar(Xsum,
#                                     Ysum,
#                                     xerr=Xsum_err,
#                                     yerr=Ysum_err,
#                                     markersize=10.,
#                                     linestyle='none',
#                                     alpha=0.5,
#                                     capsize=3,
#                                     marker='s',
#                                     ecolor='b',
#                                     markeredgecolor='b',
#                                     markerfacecolor='none',
#                                     label='Barred (global)')
#                         pl.errorbar(X,
#                                     Y,
#                                     xerr=Xerr,
#                                     yerr=Yerr,
#                                     markersize=5.,
#                                     linestyle='none',
#                                     alpha=0.2,
#                                     capsize=3,
#                                     markeredgecolor='b',
#                                     ecolor='b',
#                                     markerfacecolor='none',
#                                     label='Barred (pixels)')
#                         bb = 1
#             else:
#                 if Xsum > 4 and Xsum < 12:
#                     pl.errorbar(Xsum,
#                                 Ysum,
#                                 xerr=Xsum_err,
#                                 yerr=Ysum_err,
#                                 markersize=10.,
#                                 linestyle='none',
#                                 alpha=0.5,
#                                 capsize=3,
#                                 marker='s',
#                                 ecolor='b',
#                                 markeredgecolor='b',
#                                 markerfacecolor='none')
#                     pl.errorbar(X,
#                                 Y,
#                                 xerr=Xerr,
#                                 yerr=Yerr,
#                                 markersize=5.,
#                                 linestyle='none',
#                                 alpha=0.2,
#                                 capsize=3,
#                                 markeredgecolor='b',
#                                 ecolor='b',
#                                 markerfacecolor='none')
#         if bam[names == g] == 'A':
#             if aa == 0:
#                 if g != 'UGC09476':
#                     if Xsum > 4 and Xsum < 12:
#                         pl.errorbar(Xsum,
#                                     Ysum,
#                                     xerr=Xsum_err,
#                                     yerr=Ysum_err,
#                                     markersize=10.,
#                                     linestyle='none',
#                                     alpha=0.5,
#                                     capsize=3,
#                                     marker='s',
#                                     ecolor='g',
#                                     markeredgecolor='g',
#                                     markerfacecolor='none',
#                                     label='Unbarred (global)')
#                         pl.errorbar(X,
#                                     Y,
#                                     xerr=Xerr,
#                                     yerr=Yerr,
#                                     markersize=5.,
#                                     linestyle='none',
#                                     alpha=0.2,
#                                     capsize=3,
#                                     markeredgecolor='g',
#                                     ecolor='g',
#                                     markerfacecolor='none',
#                                     label='Unbarred (pixels)')
#                         aa = 1
#             else:
#                 if Xsum > 4 and Xsum < 12:
#                     pl.errorbar(Xsum,
#                                 Ysum,
#                                 xerr=Xsum_err,
#                                 yerr=Ysum_err,
#                                 markersize=10.,
#                                 linestyle='none',
#                                 alpha=0.5,
#                                 capsize=3,
#                                 marker='s',
#                                 ecolor='g',
#                                 markeredgecolor='g',
#                                 markerfacecolor='none')
#                     pl.errorbar(X,
#                                 Y,
#                                 xerr=Xerr,
#                                 yerr=Yerr,
#                                 markersize=5.,
#                                 linestyle='none',
#                                 alpha=0.2,
#                                 capsize=3,
#                                 markeredgecolor='g',
#                                 ecolor='g',
#                                 markerfacecolor='none')
#
#     logL12_all[
#         i] = Ysum  #np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise']-np.log10(Lsun)+np.log10(WISE_COR))))
#     logLCO_all[i] = Xsum  # np.log10(np.sum( 10**(gal_dict[g]['log10_LCO'])))
#
#     logL12_all_err[
#         i] = Ysum_err  #np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise']-np.log10(Lsun)+np.log10(WISE_COR))))
#     logLCO_all_err[
#         i] = Xsum_err  # np.log10(np.sum( 10**(gal_dict[g]['log10_LCO'])))
#
#     logL12_all_pix += Y[(np.isinf(X) == False) & (np.isinf(Y) == False)].tolist(
#     )  #np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise']-np.log10(Lsun)+np.log10(WISE_COR))))
#     logLCO_all_pix += X[(np.isinf(X) == False)
#                         & (np.isinf(Y) == False)].tolist(
#                         )  # np.log10(np.sum( 10**(gal_dict[g]['log10_LCO'])))
#
#     logL12_all_pix_err += Yerr[(np.isinf(X) == False)
#                                & (np.isinf(Y) == False)].tolist()
#     logLCO_all_pix_err += Xerr[(np.isinf(X) == False)
#                                & (np.isinf(Y) == False)].tolist()
#
#     if bam[names == g] == 'A':
#         barred_unbarred_pix += ['A'] * X[(np.isinf(X) == False)
#                                          & (np.isinf(Y) == False)].size
#         barred_unbarred_global.append('A')
#     else:
#         barred_unbarred_pix += ['B'] * X[(np.isinf(X) == False)
#                                          & (np.isinf(Y) == False)].size
#         barred_unbarred_global.append('B')
#
#     i += 1
#
# barred_unbarred_global = np.array(barred_unbarred_global)
# barred_unbarred_pix = np.array(barred_unbarred_pix)
#
# pl.xlabel(r"$\log \> L(\mathrm{CO(1-0)})$ [K km/s pc$^2$]", fontsize=15)
# pl.ylabel(r"$\log \> L(12\mu\mathrm{m})$ [L$_\odot$]", fontsize=15)
# pl.plot(np.linspace(5, 10, 10), (np.linspace(5, 10, 10) + 0.06) / 0.97,
#         label='Gao+ (in prep.)')
#
# # Fit to global values (barred)
# Xsum = logLCO_all[barred_unbarred_global == 'B']
# Ysum = logL12_all[barred_unbarred_global == 'B']
# cfit_sum_barred = curve_fit(
#     fun, Xsum[(np.isinf(Xsum) == False) & (np.isinf(Ysum) == False)],
#     Ysum[(np.isinf(Xsum) == False) & (np.isinf(Ysum) == False)])
# pfit_sum_barred = np.poly1d(cfit_sum_barred[0])
#
# # Fit to global values (unbarred)
# Xsum = logLCO_all[barred_unbarred_global == 'A']
# Ysum = logL12_all[barred_unbarred_global == 'A']
# cfit_sum_unbarred = curve_fit(
#     fun, Xsum[(np.isinf(Xsum) == False) & (np.isinf(Ysum) == False)],
#     Ysum[(np.isinf(Xsum) == False) & (np.isinf(Ysum) == False)])
# pfit_sum_unbarred = np.poly1d(cfit_sum_unbarred[0])
#
# # Fit to pixels (barred)
# Xpix = np.array(logLCO_all_pix)[barred_unbarred_pix == 'B']
# Ypix = np.array(logL12_all_pix)[barred_unbarred_pix == 'B']
# Xpix = Xpix[Ypix > 5]
# Ypix = Ypix[Ypix > 5]
# # pfit_pix = np.poly1d(np.polyfit(Xpix[(np.isinf(Xpix)==False)&(np.isinf(Ypix)==False)], Ypix[(np.isinf(Xpix)==False)&(np.isinf(Ypix)==False)], 1))
# # pfit_pix = np.poly1d(np.polyfit(Xpix[(np.isinf(Xpix)==False)&(np.isinf(Ypix)==False)], Ypix[(np.isinf(Xpix)==False)&(np.isinf(Ypix)==False)], 1))
# cfit_pix_barred = curve_fit(
#     fun, Xpix[(np.isinf(Xpix) == False) & (np.isinf(Ypix) == False)],
#     Ypix[(np.isinf(Xpix) == False) & (np.isinf(Ypix) == False)])
# pfit_pix_barred = np.poly1d(cfit_pix_barred[0])
#
# # Fit to pixels (unbarred)
# Xpix = np.array(logLCO_all_pix)[barred_unbarred_pix == 'A']
# Ypix = np.array(logL12_all_pix)[barred_unbarred_pix == 'A']
# Xpix = Xpix[Ypix > 5]
# Ypix = Ypix[Ypix > 5]
# # pfit_pix = np.poly1d(np.polyfit(Xpix[(np.isinf(Xpix)==False)&(np.isinf(Ypix)==False)], Ypix[(np.isinf(Xpix)==False)&(np.isinf(Ypix)==False)], 1))
# # pfit_pix = np.poly1d(np.polyfit(Xpix[(np.isinf(Xpix)==False)&(np.isinf(Ypix)==False)], Ypix[(np.isinf(Xpix)==False)&(np.isinf(Ypix)==False)], 1))
# cfit_pix_unbarred = curve_fit(
#     fun, Xpix[(np.isinf(Xpix) == False) & (np.isinf(Ypix) == False)],
#     Ypix[(np.isinf(Xpix) == False) & (np.isinf(Ypix) == False)])
# pfit_pix_unbarred = np.poly1d(cfit_pix_unbarred[0])
#
# pl.plot(np.linspace(5, 10, 10),
#         pfit_sum_barred(np.linspace(5, 10, 10)),
#         label='Global fit (barred)',
#         c='b')
# pl.plot(np.linspace(5, 10, 10),
#         pfit_sum_unbarred(np.linspace(5, 10, 10)),
#         label='Global fit (unbarred)',
#         c='g')
# pl.plot(np.linspace(5, 10, 10),
#         pfit_pix_barred(np.linspace(5, 10, 10)),
#         label='Pixel-pixel fit (barred)',
#         c='b',
#         linestyle='--')
# pl.plot(np.linspace(5, 10, 10),
#         pfit_pix_unbarred(np.linspace(5, 10, 10)),
#         label='Pixel-pixel fit (unbarred)',
#         c='g',
#         linestyle='--')
#
# pl.legend(loc='best', fontsize='x-small')
#
# pl.xlim(5, 10)
# pl.ylim(5, 10)
#
# # pl.savefig('/Users/ryan/Dropbox/mac/wise_w3_vs_co/l12_lco.pdf')
# pl.savefig(
#     '/Users/ryan/Dropbox/mac/wise_w3_vs_co/l12_lco_barred_unbarred_29Aug19.pdf'
# )
#
# pl.figure()
# leda_barred = list(leda['Name'][leda['ledaBar'] == 'B'])
#
# WISE_COR = 2 * 2.31e-11  #2.31e-11
# logL12_all = np.zeros(len(gal_dict.keys()))
# logLCO_all = np.zeros(len(gal_dict.keys()))
# barred_unbarred = np.zeros(len(gal_dict.keys()))
# logL12_all_pix_b = []
# logLCO_all_pix_b = []
# logL12_all_pix_nb = []
# logLCO_all_pix_nb = []
#
# fit_result = np.zeros([len(gal_dict.keys()), 2])
# i = 0
# pp = 0
# lpix = 'Pixels'
# lglob = 'Total'
# lline = 'Pixel-pixel fits'
# for g in gal_dict.keys():
#     ltemp_pix = ''
#     ltemp_glob = ''
#     ltemp_line = ''
#     if pp == 0:
#         ltemp_pix = lpix
#         ltemp_glob = lglob
#         ltemp_line = lline
#         pp = 1
#
#     cc = 'g'
#     barred_unbarred[i] = 0
#     if g in leda_barred:
#         cc = 'b'
#         barred_unbarred[i] = 1
#
#     X = gal_dict[g]['log10_LCO']
#     Y = gal_dict[g]['log10_Lwise'] - np.log10(Lsun) + np.log10(WISE_COR)
#     Xsum = np.log10(np.sum(10**(gal_dict[g]['log10_LCO'])))
#     Ysum = np.log10(
#         np.sum(10**(gal_dict[g]['log10_Lwise'] - np.log10(Lsun) +
#                     np.log10(WISE_COR))))
#     if g != 'UGC09476':
#         pl.scatter(Xsum,
#                    Ysum,
#                    s=10.,
#                    edgecolor=cc,
#                    marker='s',
#                    facecolor='none',
#                    label=ltemp_glob)
#         pl.scatter(X,
#                    Y,
#                    s=5.,
#                    alpha=0.5,
#                    edgecolor=cc,
#                    facecolor='none',
#                    label=ltemp_pix)
#
#     logL12_all[
#         i] = Ysum  #np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise']-np.log10(Lsun)+np.log10(WISE_COR))))
#     logLCO_all[i] = Xsum  # np.log10(np.sum( 10**(gal_dict[g]['log10_LCO'])))
#
#     if barred_unbarred[i] == 1:
#         logL12_all_pix_b += Y[(np.isinf(X) == False) & (
#             np.isinf(Y) == False
#         )].tolist(
#         )  #np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise']-np.log10(Lsun)+np.log10(WISE_COR))))
#         logLCO_all_pix_b += X[
#             (np.isinf(X) == False) & (np.isinf(Y) == False)].tolist(
#             )  # np.log10(np.sum( 10**(gal_dict[g]['log10_LCO'])))
#     else:
#         logL12_all_pix_nb += Y[(np.isinf(X) == False) & (
#             np.isinf(Y) == False
#         )].tolist(
#         )  #np.log10(np.sum( 10**(gal_dict[g]['log10_Lwise']-np.log10(Lsun)+np.log10(WISE_COR))))
#         logLCO_all_pix_nb += X[
#             (np.isinf(X) == False) & (np.isinf(Y) == False)].tolist(
#             )  # np.log10(np.sum( 10**(gal_dict[g]['log10_LCO'])))
#
#     i += 1
#
# pl.xlabel(r"$\log \> L(\mathrm{CO(1-0)})$ [K km/s pc$^2$]", fontsize=15)
# pl.ylabel(r"$\log \> L(12\mu\mathrm{m})$ [L$_\odot$]", fontsize=15)
# pl.plot(np.linspace(5, 10, 10), (np.linspace(5, 10, 10) + 0.06) / 0.97,
#         label='Gao+ (in prep.)')
#
# Xsum = logLCO_all[barred_unbarred == 1]
# Ysum = logL12_all[barred_unbarred == 1]
# # pfit_sum = np.poly1d(np.polyfit(Xsum[(np.isinf(Xsum)==False)&(np.isinf(Ysum)==False)], Ysum[(np.isinf(Xsum)==False)&(np.isinf(Ysum)==False)], 1))
# cfit_sum_b = curve_fit(
#     fun, Xsum[(np.isinf(Xsum) == False) & (np.isinf(Ysum) == False)],
#     Ysum[(np.isinf(Xsum) == False) & (np.isinf(Ysum) == False)])
# pfit_sum_b = np.poly1d(cfit_sum_b[0])
#
# Xsum = logLCO_all[barred_unbarred == 0]
# Ysum = logL12_all[barred_unbarred == 0]
# # pfit_sum = np.poly1d(np.polyfit(Xsum[(np.isinf(Xsum)==False)&(np.isinf(Ysum)==False)], Ysum[(np.isinf(Xsum)==False)&(np.isinf(Ysum)==False)], 1))
# cfit_sum_nb = curve_fit(
#     fun, Xsum[(np.isinf(Xsum) == False) & (np.isinf(Ysum) == False)],
#     Ysum[(np.isinf(Xsum) == False) & (np.isinf(Ysum) == False)])
# pfit_sum_nb = np.poly1d(cfit_sum_nb[0])
#
# # Xsum = logLCO_all
# # Ysum = logL12_all
# # # pfit_sum = np.poly1d(np.polyfit(Xsum[(np.isinf(Xsum)==False)&(np.isinf(Ysum)==False)], Ysum[(np.isinf(Xsum)==False)&(np.isinf(Ysum)==False)], 1))
# # cfit_sum = curve_fit(fun, Xsum[(np.isinf(Xsum)==False)&(np.isinf(Ysum)==False)], Ysum[(np.isinf(Xsum)==False)&(np.isinf(Ysum)==False)])
# # pfit_sum = np.poly1d(cfit_sum[0])
#
# Xpix = np.array(logLCO_all_pix_nb)
# Ypix = np.array(logL12_all_pix_nb)
# Xpix = Xpix[Ypix > 5]
# Ypix = Ypix[Ypix > 5]
# cfit_pix_nb = curve_fit(
#     fun, Xpix[(np.isinf(Xpix) == False) & (np.isinf(Ypix) == False)],
#     Ypix[(np.isinf(Xpix) == False) & (np.isinf(Ypix) == False)])
# pfit_pix_nb = np.poly1d(cfit_pix_nb[0])
#
# Xpix = np.array(logLCO_all_pix_b)
# Ypix = np.array(logL12_all_pix_b)
# Xpix = Xpix[Ypix > 5]
# Ypix = Ypix[Ypix > 5]
# cfit_pix_b = curve_fit(
#     fun, Xpix[(np.isinf(Xpix) == False) & (np.isinf(Ypix) == False)],
#     Ypix[(np.isinf(Xpix) == False) & (np.isinf(Ypix) == False)])
# pfit_pix_b = np.poly1d(cfit_pix_b[0])
#
# pl.plot(np.linspace(5, 10, 10),
#         pfit_sum_b(np.linspace(5, 10, 10)),
#         label='Global fit (barred)',
#         c='k')
# pl.plot(np.linspace(5, 10, 10),
#         pfit_pix_b(np.linspace(5, 10, 10)),
#         label='Pixel-pixel fit (barred)',
#         c='k',
#         linestyle='--')
#
# pl.plot(np.linspace(5, 10, 10),
#         pfit_sum_nb(np.linspace(5, 10, 10)),
#         label='Global fit (unbarred)',
#         c='r')
# pl.plot(np.linspace(5, 10, 10),
#         pfit_pix_nb(np.linspace(5, 10, 10)),
#         label='Pixel-pixel fit (unbarred)',
#         c='r',
#         linestyle='--')
#
# pl.legend(loc='best')
#
# pl.savefig('/Users/ryan/Dropbox/mac/wise_w3_vs_co/l12_lco_barred_unbarred.pdf')
#
# array = np.ones(image_w3.shape)
# array = rescale(array, 1.375 / 6., anti_aliasing=True)
# array = resize(array, (mom0[0].data.shape[0], mom0[0].data.shape[1]),
#                anti_aliasing=True)
#
# # Write CSV file for matching with NSA
#
# edge_id_1 = []
# edge_id_2 = []
# RA_DEG = []
# DEC_DEG = []
#
# for galname in galnames:
#     gal = galname
#     if galname[:3] == 'UGC':
#         if len(
#                 glob.glob(
#                     '/Users/ryan/venus/shared_data/edge/co_cubes_w3conv_aniano_6Jun19/%s_co_smooth_wise_rebin6_mom0_Sun.pbcor.K.fits'
#                     % (galname, ))) == 0:
#             gal = 'UGC0' + galname[3:]
#         if len(
#                 glob.glob(
#                     '/Users/ryan/venus/shared_data/edge/co_cubes_w3conv_aniano_6Jun19/%s_co_smooth_wise_rebin6_mom0_Sun.pbcor.K.fits'
#                     % (gal, ))) == 0:
#             gal = 'UGC00' + galname[3:]
#     if galname[:2] == 'IC':
#         if len(
#                 glob.glob(
#                     '/Users/ryan/venus/shared_data/edge/co_cubes_w3conv_aniano_6Jun19/%s_co_smooth_wise_rebin6_mom0_Sun.pbcor.K.fits'
#                     % (galname, ))) == 0:
#             gal = 'IC0' + galname[2:]
#     if galname[:3] == 'NGC':
#         if len(
#                 glob.glob(
#                     '/Users/ryan/venus/shared_data/edge/co_cubes_w3conv_aniano_6Jun19/%s_co_smooth_wise_rebin6_mom0_Sun.pbcor.K.fits'
#                     % (galname, ))) == 0:
#             gal = 'NGC0' + galname[3:]
#
#     gtmp = galname
#     if galname == 'NGC4211NED02':
#         gtmp = 'NGC4211B'
#     c = coordinates.get_icrs_coordinates(gtmp)
#     edge_id_1.append(galname)
#     edge_id_2.append(gal)
#     RA_DEG.append(c.ra.deg)
#     DEC_DEG.append(c.dec.deg)
#
#
# galaxy_dict = dict( edge_id_1 = edge_id_1,\
#     edge_id_2 = edge_id_2,\
#     ra_deg = RA_DEG,\
#     dec_deg = DEC_DEG)
#
# df = pd.DataFrame(galaxy_dict)
# cols = ['edge_id_1', 'edge_id_2', 'ra_deg',\
# 'dec_deg', 'ba', 'z']
# df.to_csv("/Users/ryan/Dropbox/mac/tsinghua/edge_all_ids_radecs.csv",
#           sep=",",
#           float_format='%.9g',
#           header=True,
#           index=False,
#           columns=cols)
