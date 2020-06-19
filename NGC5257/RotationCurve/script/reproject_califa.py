import sys
import pickle
import glob
import numpy as np
import matplotlib.pyplot as pl
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import match_coordinates_sky

from reproject import reproject_exact

# galname = 'NGC5000'
# stellar mass : maps_corr[10].header (solar masses)
# Halpha flux maps_corr[11].data[15] (1E-16 erg/s/cm^2)

chi2_upper_limit = 20.

fname_l12_lco_dict = '/Users/ryan/Dropbox/mac/wise_w3_vs_co/l12_lco_dict_v5.pk'  # Without Aniano kernels, until that is fixed, and with more conservative s/n cuts in mom0 map making
gal_dict = pickle.load(open(fname_l12_lco_dict, 'rb'))


def extcorr(ha_flux, hb_flux):
    # ebv = 1.97 * np.log10((ha_flux / hb_flux).astype(np.float64) / 2.86)
    # a_ha = 3.33 * ebv
    a_ha = 5.86 * np.log10((ha_flux / hb_flux).astype(np.float64) / 2.86)
    a_ha[a_ha < 0] = 0.
    return 10**(0.4 * a_ha)


# Filename functions
def fname_mstar(galname):
    """
    Filename of the pickled stellar mass dict for a galaxy
    """
    out_dir = '/Users/ryan/venus/shared_data/califa/DR3-Niu/%s' % (galname, )
    if len(glob.glob(out_dir)) == 0:
        out_dir = '/Users/ryan/venus/shared_data/califa/DR3-V500-Niu/%s' % (
            galname, )

    out_name = out_dir + '/mstar_reproj.pk'
    return out_name


def fname_mstar_stacked(galname):
    """
    Filename of the pickled stellar mass dict for a galaxy
    """
    out_dir = '/Users/ryan/venus/shared_data/califa/DR3-Niu/%s' % (galname, )
    if len(glob.glob(out_dir)) == 0:
        out_dir = '/Users/ryan/venus/shared_data/califa/DR3-V500-Niu/%s' % (
            galname, )

    out_name = out_dir + '/mstar_stacked_reproj.pk'
    return out_name


def fname_alpha_co_stacked(galname):
    """
    Filename of the pickled stellar mass dict for a galaxy
    """
    out_dir = '/Users/ryan/venus/shared_data/califa/DR3-Niu/%s' % (galname, )
    if len(glob.glob(out_dir)) == 0:
        out_dir = '/Users/ryan/venus/shared_data/califa/DR3-V500-Niu/%s' % (
            galname, )

    out_name = out_dir + '/alphaco_stacked_reproj.pk'
    return out_name


def fname_halpha(galname):
    """
    Filename of the pickled F_Ha dict for a galaxy (contains F_Ha and error)
    """
    out_dir = '/Users/ryan/venus/shared_data/califa/DR3-Niu/%s' % (galname, )
    if len(glob.glob(out_dir)) == 0:
        out_dir = '/Users/ryan/venus/shared_data/califa/DR3-V500-Niu/%s' % (
            galname, )

    out_name = out_dir + '/halpha_reproj.pk'
    return out_name


def fname_halpha_stacked(galname):
    """
    Filename of the pickled F_Ha dict for a galaxy (contains F_Ha and error)
    """
    out_dir = '/Users/ryan/venus/shared_data/califa/DR3-Niu/%s' % (galname, )
    if len(glob.glob(out_dir)) == 0:
        out_dir = '/Users/ryan/venus/shared_data/califa/DR3-V500-Niu/%s' % (
            galname, )

    out_name = out_dir + '/halpha_stacked_reproj.pk'
    return out_name


def fname_metallicity_stacked(galname):
    """
    Filename of the pickled metallicity dict for a galaxy
    """
    out_dir = '/Users/ryan/venus/shared_data/califa/DR3-Niu/%s' % (galname, )
    if len(glob.glob(out_dir)) == 0:
        out_dir = '/Users/ryan/venus/shared_data/califa/DR3-V500-Niu/%s' % (
            galname, )

    out_name = out_dir + '/metallicity_stacked_reproj.pk'
    return out_name


def bpt(opened_fits):
    f_ha = opened_fits[11].data[15]
    f_ha_err = opened_fits[15].data[15]
    f_hb = opened_fits[11].data[8]
    f_hb_err = opened_fits[15].data[8]
    f_o3 = opened_fits[11].data[10]  # OIII 5007
    f_o3_err = opened_fits[15].data[10]
    f_n2 = opened_fits[11].data[16]  # NII 6583
    f_n2_err = opened_fits[15].data[16]

    sn_cut = 2.  # require s/n > this number for all 4 lines
    sn_ha = np.zeros(f_ha.shape)
    sn_ha[np.nonzero(f_ha_err)] = f_ha[np.nonzero(f_ha_err)] / f_ha_err[
        np.nonzero(f_ha_err)]

    sn_hb = np.zeros(f_hb.shape)
    sn_hb[np.nonzero(f_hb_err)] = f_hb[np.nonzero(f_hb_err)] / f_hb_err[
        np.nonzero(f_hb_err)]

    sn_o3 = np.zeros(f_o3.shape)
    sn_o3[np.nonzero(f_o3_err)] = f_o3[np.nonzero(f_o3_err)] / f_o3_err[
        np.nonzero(f_o3_err)]

    sn_n2 = np.zeros(f_n2.shape)
    sn_n2[np.nonzero(f_n2_err)] = f_n2[np.nonzero(f_n2_err)] / f_n2_err[
        np.nonzero(f_n2_err)]

    sngood = np.ones(f_ha.shape)
    sngood[sn_ha < sn_cut] = 0.
    sngood[sn_hb < sn_cut] = 0.
    sngood[sn_o3 < sn_cut] = 0.
    sngood[sn_n2 < sn_cut] = 0.

    bpt_x = np.zeros(f_ha.shape)
    bpt_y = np.zeros(f_ha.shape)

    bpt_x[sngood != 0] = np.log10((f_n2[sngood != 0]) / (f_ha[sngood != 0]))
    bpt_y[sngood != 0] = np.log10((f_o3[sngood != 0]) / (f_hb[sngood != 0]))

    k03 = lambda bpt_x: 0.61 / (bpt_x - 0.05) + 1.3
    kewley = lambda bpt_x: 0.61 / (bpt_x - 0.47) + 1.19
    kewley_bpt_y = kewley(bpt_x)
    yes_agn = np.zeros(f_ha.shape)
    yes_agn[bpt_x >= 0.47] = 1
    yes_agn[(kewley_bpt_y <= bpt_y) & (sngood != 0)] = 1
    # yes_agn[sngood == 0] = 1.
    return yes_agn


# Array getters
def get_map_halpha(opened_fits, do_balmer=False):
    """
    Args:
        opened_fits : Fits file containing Halpha
        do_balmer : Apply Balmer decrement to the map before returning
    Returns:
        map_line_flux (2d array) : Line flux in units of 10^-16 erg/s/cm^2
    """
    map_line_flux = opened_fits[11].data[15]

    if do_balmer == True:
        ew_ha = opened_fits[14].data[15].astype(np.float64)

        f_ha = opened_fits[11].data[15].astype(np.float64)
        f_ha_err = opened_fits[15].data[15].astype(np.float64)

        f_hb = opened_fits[11].data[8].astype(np.float64)
        f_hb_err = opened_fits[15].data[8].astype(np.float64)

        # snr_ha = f_ha / f_ha_err
        # snr_ha[np.isnan(snr_ha)] = 0
        # snr_ha[np.isinf(snr_ha)] = 0
        #
        # snr_ha[ew_ha < 6] = 0
        #
        # snr_hb = f_hb / f_hb_err
        # snr_hb[np.isnan(snr_hb)] = 0
        # snr_hb[np.isinf(snr_hb)] = 0

        ext_corr = extcorr(f_ha, f_hb)
        map_line_flux *= ext_corr

        yes_agn = bpt(opened_fits)
        map_line_flux[yes_agn == 1] = np.nan
        map_line_flux[ew_ha < 6] = np.nan

    else:
        chi2 = opened_fits[2].data
        map_line_flux[chi2 > chi2_upper_limit] = 0

    map_line_flux[np.isinf(map_line_flux)] = np.nan
    return map_line_flux


def get_map_halpha_err(opened_fits, do_balmer=False):
    """
    Args:
        opened_fits : Fits file containing Halpha
        do_balmer : Apply Balmer decrement to the map before returning
    Returns:
        map_line_flux (2d array) : Line flux error in units of 10^-16 erg/s/cm^2
    """
    map_line_flux = opened_fits[15].data[15]

    if do_balmer == True:
        ew_ha = opened_fits[14].data[15].astype(np.float64)

        f_ha = opened_fits[11].data[15].astype(np.float64)
        f_ha_err = opened_fits[15].data[15].astype(np.float64)

        f_hb = opened_fits[11].data[8].astype(np.float64)
        f_hb_err = opened_fits[15].data[8].astype(np.float64)

        # snr_ha = f_ha / f_ha_err
        # snr_ha[np.isnan(snr_ha)] = 0
        # snr_ha[np.isinf(snr_ha)] = 0
        #
        # snr_ha[ew_ha < 6] = 0
        #
        # snr_hb = f_hb / f_hb_err
        # snr_hb[np.isnan(snr_hb)] = 0
        # snr_hb[np.isinf(snr_hb)] = 0

        ext_corr = extcorr(f_ha, f_hb)
        map_line_flux *= ext_corr

        yes_agn = bpt(opened_fits)
        map_line_flux[yes_agn == 1] = np.nan
        map_line_flux[ew_ha < 6] = np.nan
    else:
        chi2 = opened_fits[2].data
        map_line_flux[chi2 > chi2_upper_limit] = 0
    map_line_flux[np.isinf(map_line_flux)] = np.nan
    return map_line_flux


def get_map_stellar_mass(opened_fits):
    """
    Args:
        opened_fits : Fits file containing extinction-corrected stellar mass
    Returns:
        map_stellar_mass (2d array) : Stellar mass in units of Msun
    """
    map_stellar_mass = opened_fits[10].data
    chi2 = opened_fits[2].data
    map_stellar_mass[chi2 > chi2_upper_limit] = np.nan
    return map_stellar_mass


def reproj_califa(galname, get_map, do_balmer=False, is_uncertainty=False):
    """
    Args:
        get_map : Function that will take an opened fits file and give the
            desired map array.
        do_balmer : Whether or not to apply the Balmer decrement to the map.
            Note that if do_balmer is True, get_map must also have a do_balmer
            argument (True/False).
        is_uncertainty : Whether this is a map of uncertainty or not (True/False)
    Returns:
        map_param_reproj : map of the desired parameter reprojected
            into the same pixel size and array size as the corresponding CO map.
    """
    fname_niu_map = '/Users/ryan/venus/shared_data/califa/DR3-Niu/%s/califa-%s-ppxf-Maps-corr.fits.gz' % (
        galname, galname)
    if do_balmer == True:
        fname_niu_map = '/Users/ryan/venus/shared_data/califa/DR3-Niu/%s/califa-%s-ppxf-Maps.fits' % (
            galname, galname)
    if len(glob.glob(fname_niu_map)) == 0:
        fname_niu_map = '/Users/ryan/venus/shared_data/califa/DR3-V500-Niu/%s/califa-%s-ppxf-Maps-corr.fits.gz' % (
            galname, galname)
        if do_balmer == True:
            fname_niu_map = '/Users/ryan/venus/shared_data/califa/DR3-V500-Niu/%s/califa-%s-ppxf-Maps.fits.gz' % (
                galname, galname)

    if len(glob.glob(fname_niu_map)) == 0:
        print("This galaxy doesn't have a V500 data cube")
        return 0

    maps_corr = fits.open(fname_niu_map)
    if do_balmer == True:
        map_param = get_map(maps_corr, do_balmer)
    else:
        map_param = get_map(maps_corr)

    if is_uncertainty == True:
        map_param = map_param**2

    w = wcs.WCS(maps_corr[0].header)
    w = w.dropaxis(2)

    co_map, co_header = fits.getdata(
        '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
        % (galname, ),
        header=True)
    co_wcs = wcs.WCS(co_header)
    co_shape = co_map.shape

    map_param[map_param == 0] = np.nan

    mask = np.zeros(map_param.shape)
    mask[np.isnan(map_param)] = 1
    mask[map_param == 0] = 1
    mask = mask.astype(bool)

    map_param_reproj, footprint = reproject_exact((map_param, w),
                                                  co_wcs,
                                                  co_shape,
                                                  parallel=False)
    if is_uncertainty == True:
        map_param_reproj = np.sqrt(map_param_reproj)
    area_fac = (co_wcs.wcs.cdelt[1] / w.wcs.cd[1, 1])**2
    map_param_reproj *= area_fac

    # Check for not fully sampled pixels
    msk = np.ones(map_param.shape)
    msk[np.isnan(map_param)] = 0
    msk[map_param == 0] = 0.
    msk_reproj, footprint = reproject_exact((msk, w),
                                            co_wcs,
                                            co_shape,
                                            parallel=False)
    # msk_reproj *= area_fac
    msk_reproj[msk_reproj != 1.] = 0.

    # Mask these pixels in the reprojected map
    map_param_reproj[msk_reproj == 0] = np.nan
    return map_param_reproj  #, msk, msk_reproj


def reproj_califa_mapping(galname):
    """
    Returns:
        res (dict) : 'pixel_map_co_proj' is the array of indices for each pixel
            in the desired output projection (for me this is 6 arcsec pixels).
            Each element of 'pixel_map' (an array, the same shape as the
            input map) gives the index of the corresponding output pixel.

    """
    fname_niu_map = '/Users/ryan/venus/shared_data/califa/DR3-Niu/%s/califa-%s-ppxf-Maps-corr.fits.gz' % (
        galname, galname)
    if len(glob.glob(fname_niu_map)) == 0:
        fname_niu_map = '/Users/ryan/venus/shared_data/califa/DR3-V500-Niu/%s/califa-%s-ppxf-Maps-corr.fits.gz' % (
            galname, galname)

    if len(glob.glob(fname_niu_map)) == 0:
        print("This galaxy doesn't have a V500 data cube")
        # return 0

    maps_corr = fits.open(fname_niu_map)

    w = wcs.WCS(maps_corr[0].header)
    w = w.dropaxis(2)

    co_map, co_header = fits.getdata(
        '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
        % (galname, ),
        header=True)
    co_wcs = wcs.WCS(co_header)
    co_shape = co_map.shape

    # Begin new method (currently pretty slow)
    ny_co, nx_co = co_shape
    xs_co, ys_co = np.meshgrid(np.arange(nx_co), np.arange(ny_co))
    coords_co = pixel_to_skycoord(xs_co, ys_co, co_wcs)
    coords_co_flat = coords_co.flatten()
    pixel_labels_co = np.arange(xs_co.size)  # .reshape(co_shape)

    ny_ifu, nx_ifu = maps_corr[1].data.shape
    xs_ifu, ys_ifu = np.meshgrid(np.arange(nx_ifu), np.arange(ny_ifu))
    coords_ifu = pixel_to_skycoord(xs_ifu, ys_ifu, w)
    coords_ifu = coords_ifu.transform_to('fk5')
    pixel_map_arr = np.full((nx_ifu, ny_ifu), np.nan).flatten()

    i_ifu = 0
    npix_ifu = coords_ifu.flatten().size
    dra, ddec = np.zeros(npix_ifu), np.zeros(npix_ifu)
    i_co_all = np.zeros(npix_ifu, dtype=int)
    i_co, d2d, d3d = match_coordinates_sky(coords_ifu.flatten(),
                                           coords_co_flat)
    dra, ddec = (coords_ifu.flatten()).spherical_offsets_to(
        coords_co_flat[i_co])
    print(dra, ddec)
    dra = dra.arcsec
    ddec = ddec.arcsec

    # for c in coords_ifu.flatten():
    #     dra_i, ddec_i = c.spherical_offsets_to(coords_co_flat[i_co[i_ifu]])
    #     dra[i_ifu] = dra_i.arcsec
    #     ddec[i_ifu] = ddec_i.arcsec
    #     i_co_all[i_ifu] = i_co[i_ifu]
    #     i_ifu += 1

    good = (-3 <= dra) & (dra < 3) & (-3 <= ddec) & (ddec < 3)

    # if (-3 <= dra) and (dra < 3):
    #     if (-3 <= ddec) and (ddec < 3):
    #         pixel_map_arr[i_ifu] = pixel_labels_co[i_co]
    # i_ifu += 1

    pixel_map_arr[good] = pixel_labels_co[i_co[good]]

    pixel_map_arr = pixel_map_arr.reshape(maps_corr[1].data.shape).astype(
        np.int)
    pixel_labels_co = pixel_labels_co.reshape(co_shape)
    # End new method

    # mask = np.ones(maps_corr[1].data.shape)
    # mask_reproj, footprint_v1 = reproject_exact((mask, w),
    #                                          co_wcs,
    #                                          co_shape,
    #                                          parallel=False)
    # mask_reproj_v2 = mask_reproj.copy()
    # mask_reproj[mask_reproj != 1] = -1 # np.nan
    # mask_reproj[mask_reproj == 1] = np.arange(
    #     mask_reproj[mask_reproj == 1].size).reshape(
    #         mask_reproj[mask_reproj == 1].shape)
    #
    # mask_reproj_back, footprint = reproject_exact((mask_reproj, co_wcs),
    #                                               w,
    #                                               mask.shape,
    #                                               parallel=False)
    # mask_reproj_back_before_rounding = mask_reproj_back.copy()
    # mask_reproj_back = mask_reproj_back.astype(int)
    res = dict()
    # res['mask'] = mask
    res['pixel_map_co_proj'] = pixel_labels_co  #mask_reproj
    res['pixel_map'] = pixel_map_arr  #mask_reproj_back

    # mask_reproj_back_again, footprint = reproject_exact((mask_reproj_back, w),
    #                                               co_wcs,
    #                                               co_shape,
    #                                               parallel=False)
    #
    # ignore = np.unique( mask_reproj[np.abs(mask_reproj-mask_reproj_back_again)>0.5] )
    # #
    # # mask_reproj_v2[(mask_reproj_v2.astype(int) != 1) | (footprint_v1.astype(int) != 1)] = np.nan
    # # mask_reproj_v2[~np.isnan(mask_reproj_v2)] = 1
    # # # mask_reproj_v2[mask_reproj_v2 == 1] = np.arange(
    # # #     mask_reproj_v2[mask_reproj_v2 == 1].size).reshape(
    # # #         mask_reproj_v2[mask_reproj_v2 == 1].shape)
    # #
    # # mask_reproj_back_v2, footprint = reproject_exact((mask_reproj_v2, co_wcs),
    # #                                               w,
    # #                                               mask.shape,
    # #                                               parallel=False)
    # # mask_reproj_back_v2 = mask_reproj_back_v2.astype(int)
    # #
    # # mask_reproj_back_before_rounding > mask_reproj_back
    #
    # res['ignore'] = ignore #np.unique(mask_reproj_back[np.isnan(mask_reproj_back_v2) | (mask_reproj_back_v2.astype(int) != 1)])
    return res


def reproj_any_to_co(galname, wcs_in, map_in_shape):
    """
    Args:
        wcs_in : WCS of the input map
        map_in_shape: shape of the input map
    Returns:
        dictionary of the pixel mapping between input and output
    """
    # fname_niu_map = '/Users/ryan/venus/shared_data/califa/DR3-Niu/%s/califa-%s-ppxf-Maps-corr.fits.gz' % (
    #     galname, galname)
    # if len(glob.glob(fname_niu_map)) == 0:
    #     fname_niu_map = '/Users/ryan/venus/shared_data/califa/DR3-V500-Niu/%s/califa-%s-ppxf-Maps-corr.fits.gz' % (
    #         galname, galname)
    #
    # if len(glob.glob(fname_niu_map)) == 0:
    #     print("This galaxy doesn't have a V500 data cube")
    #     # return 0
    #
    # maps_corr = fits.open(fname_niu_map)
    #
    # w = wcs.WCS(maps_corr[0].header)
    # w = w.dropaxis(2)

    # w = wcs_in

    co_map, co_header = fits.getdata(
        '/Users/ryan/Dropbox/mac/wise_w3_vs_co/%s_co_smooth_wise_v2_rebin6_mom0.fits'
        % (galname, ),
        header=True)
    co_wcs = wcs.WCS(co_header)
    co_shape = co_map.shape

    # Begin new method (currently pretty slow)
    ny_co, nx_co = co_shape
    xs_co, ys_co = np.meshgrid(np.arange(nx_co), np.arange(ny_co))
    coords_co = pixel_to_skycoord(xs_co, ys_co, co_wcs)
    coords_co_flat = coords_co.flatten()
    pixel_labels_co = np.arange(xs_co.size)  # .reshape(co_shape)

    ny_in, nx_in = map_in_shape  #maps_corr[1].data.shape
    xs_in, ys_in = np.meshgrid(np.arange(nx_in), np.arange(ny_in))
    coords_in = pixel_to_skycoord(xs_in, ys_in, wcs_in)
    coords_in = coords_in.transform_to('fk5')
    pixel_map_arr = np.full((nx_in, ny_in), np.nan).flatten()

    i_in = 0
    npix_in = coords_in.flatten().size
    dra, ddec = np.zeros(npix_in), np.zeros(npix_in)
    i_co_all = np.zeros(npix_in, dtype=int)
    i_co, d2d, d3d = match_coordinates_sky(coords_in.flatten(), coords_co_flat)
    dra, ddec = (coords_in.flatten()).spherical_offsets_to(
        coords_co_flat[i_co])
    print(dra, ddec)
    dra = dra.arcsec
    ddec = ddec.arcsec

    good = (-3 <= dra) & (dra < 3) & (-3 <= ddec) & (ddec < 3)

    pixel_map_arr[good] = pixel_labels_co[i_co[good]]

    pixel_map_arr = pixel_map_arr.reshape(map_in_shape).astype(np.int)
    pixel_labels_co = pixel_labels_co.reshape(co_shape)

    res = dict()
    res['pixel_map_co_proj'] = pixel_labels_co
    res['pixel_map'] = pixel_map_arr

    return res


def reproj_halpha(galname, do_balmer=True):
    """
    Wrapper on reproj_califa to reproject H alpha flux and uncertainty,
    applying Balmer decrement by default, convert to SFR surface density
    in Msun/yr/kpc^2, and saves the result into a dictionary.
    """
    halpha = reproj_califa(galname, get_map_halpha, do_balmer=do_balmer)
    if type(halpha) == int:
        print("Not saving anything for this galaxy (%s)" % (galname, ))
        return 0
    else:
        halpha_err = reproj_califa(galname,
                                   get_map_halpha_err,
                                   do_balmer=do_balmer,
                                   is_uncertainty=True)

        gal_dict_i = gal_dict[galname]
        d_mpc = gal_dict_i['dist_Mpc'] * u.Mpc
        d_cm = np.float64((d_mpc).to(u.cm).value)
        l_ha_to_sfr = 5.3e-42
        f_ha_to_erg_s_cm2 = 1e-16
        pc2_to_kpc2 = 1e6
        pix_size_arcsec = 6.
        d_kpc = np.float64(d_mpc.to(u.kpc).value)
        pix_area_kpc2 = (pix_size_arcsec * d_kpc / 206265.)**2
        cos_i = np.cos(np.radians(gal_dict_i['incl']))

        units_factor = cos_i * 4. * np.pi * d_cm**2 * l_ha_to_sfr * f_ha_to_erg_s_cm2 / pix_area_kpc2

        res = dict()
        # res['halpha_flux'] = halpha * units_factor
        # res['halpha_flux_err'] = halpha_err * units_factor
        res['sigma_sfr'] = halpha * units_factor
        res['sigma_sfr_err'] = halpha_err * units_factor

        out_name = fname_halpha(galname)

        with open(out_name, 'wb') as p:
            pickle.dump(res, p)


def reproj_mstar(galname):
    """
    Wrapper on reproj_califa to reproject mstar,
    and saves the result into a dictionary.
    The result is stellar mass surface density in Msun/pc^2
    """
    mstar = reproj_califa(galname, get_map_stellar_mass)
    if type(mstar) == int:
        print("Not saving anything for this galaxy (%s)" % (galname, ))
        return 0
    else:
        gal_dict_i = gal_dict[galname]
        d_mpc = gal_dict_i['dist_Mpc'] * u.Mpc
        d_cm = np.float64((d_mpc).to(u.cm).value)
        l_ha_to_sfr = 5.3e-42
        f_ha_to_erg_s_cm2 = 1e-16
        pc2_to_kpc2 = 1e6
        pix_size_arcsec = 6.
        d_pc = np.float64(d_mpc.to(u.pc).value)
        pix_area_pc2 = (pix_size_arcsec * d_pc / 206265.)**2
        cos_i = np.cos(np.radians(gal_dict_i['incl']))

        res = dict()
        res['mstar'] = mstar * cos_i / pix_area_pc2

        out_name = fname_mstar(galname)

        with open(out_name, 'wb') as p:
            pickle.dump(res, p)


def get_array_reproj_stacked(galname, quantity, flux_ext=None):
    """
    Args:
        quantity (str) : one of the keys in /Users/ryan/venus/shared_data/califa/DR3-stack/%s/%s_result.pk
        flux_ext (int) : integer index for flux arrays
    Returns:
        quantity_map : array of the parameter you want (e.g. stellar mass, etc)
    """
    stack_result = pickle.load(
        open(
            '/Users/ryan/venus/shared_data/califa/DR3-stack/%s/%s_result.pk' %
            (galname, galname), 'rb'))
    pixel_mapping = pickle.load(
        open('/Users/ryan/Dropbox/mac/wise_w3_vs_co/pixel_mapping.pk',
             'rb'))[galname]

    pixel_map_co_proj = pixel_mapping['pixel_map_co_proj'].astype(int)
    res = np.full(pixel_map_co_proj.size,
                  np.nan)  # np.ones(pixel_map_co_proj.shape) * np.nan
    quantity_arr = stack_result[quantity]
    if flux_ext is not None:
        quantity_arr = stack_result[quantity][:, flux_ext]

    # r = reproj_califa_mapping(galname)
    # ha = maps_corr[1].data
    pmap = np.unique(pixel_mapping['pixel_map'])
    cmap = pixel_mapping['pixel_map_co_proj'].flatten()
    i = 0
    for p in pmap:
        res[cmap == p] = quantity_arr[i]
        i += 1

    res = res.reshape(pixel_map_co_proj.shape)

    # # for i in range(0, quantity_arr.size):
    # ignore = pixel_mapping['ignore']
    # j = 0
    # for i in np.unique(pixel_mapping['pixel_map']):
    #     if i not in ignore:
    #         quantity_map[pixel_map_co_proj == i] = quantity_arr[j]
    #     j += 1
    # quantity_map[quantity_map == 0] = np.nan
    return res


def reproj_mstar_stacked(galname):
    """
    Stellar mass surface density map from stacked spectral fitting.

    Returns:
        array of stellar mass surface density (Msun/pc^2),
        inclination-deprojected, in CO projection (6 arcsec pixels).
    """
    if len(
            glob.glob(
                '/Users/ryan/venus/shared_data/califa/DR3-stack/%s/%s_result.pk'
                % (galname, galname))) == 0:
        print("Nothing exists for this galaxy")
        return 0

    chi2 = get_array_reproj_stacked(galname, 'chi2')
    mstar_map = get_array_reproj_stacked(galname, 'M')

    gal_dict_i = gal_dict[galname]
    d_mpc = gal_dict_i['dist_Mpc'] * u.Mpc
    d_cm = np.float64((d_mpc).to(u.cm).value)
    l_ha_to_sfr = 5.3e-42
    f_ha_to_erg_s_cm2 = 1e-16
    pc2_to_kpc2 = 1e6
    pix_size_arcsec = 6.
    d_pc = np.float64(d_mpc.to(u.pc).value)
    pix_area_pc2 = (pix_size_arcsec * d_pc / 206265.)**2
    cos_i = np.cos(np.radians(gal_dict_i['incl']))

    res = dict()
    res['mstar'] = mstar_map * cos_i / pix_area_pc2

    out_name = fname_mstar_stacked(galname)
    with open(out_name, 'wb') as p:
        pickle.dump(res, p)


#
# def bpt_stacked(galname):
#     f_ha = get_array_reproj_stacked(galname, 'Ha6563', 2)
#     f_hb = get_array_reproj_stacked(galname, 'Hb4861', 2)
#     f_o3 = get_array_reproj_stacked(galname, 'OIII5007', 2)
#     f_n2 = get_array_reproj_stacked(galname, 'NII6583', 2)
#
#     f_ha_err = get_array_reproj_stacked(galname, 'Ha6563', 6)
#     f_hb_err = get_array_reproj_stacked(galname, 'Hb4861', 6)
#     f_o3_err = get_array_reproj_stacked(galname, 'OIII5007', 6)
#     f_n2_err = get_array_reproj_stacked(galname, 'NII6583', 6)
#
#     sn_cut = 3.  # require s/n > this number for all 4 lines
#     sn_ha = np.zeros(f_ha.shape)
#     sn_ha[np.nonzero(f_ha_err)] = f_ha[np.nonzero(f_ha_err)] / f_ha_err[
#         np.nonzero(f_ha_err)]
#
#     sn_hb = np.zeros(f_hb.shape)
#     sn_hb[np.nonzero(f_hb_err)] = f_hb[np.nonzero(f_hb_err)] / f_hb_err[
#         np.nonzero(f_hb_err)]
#
#     sn_o3 = np.zeros(f_o3.shape)
#     sn_o3[np.nonzero(f_o3_err)] = f_o3[np.nonzero(f_o3_err)] / f_o3_err[
#         np.nonzero(f_o3_err)]
#
#     sn_n2 = np.zeros(f_n2.shape)
#     sn_n2[np.nonzero(f_n2_err)] = f_n2[np.nonzero(f_n2_err)] / f_n2_err[
#         np.nonzero(f_n2_err)]
#
#     sngood = np.ones(f_ha.shape)
#     sngood[sn_ha < sn_cut] = 0.
#     sngood[sn_hb < sn_cut] = 0.
#     sngood[sn_o3 < sn_cut] = 0.
#     sngood[sn_n2 < sn_cut] = 0.
#
#     bpt_x = np.zeros(f_ha.shape)
#     bpt_y = np.zeros(f_ha.shape)
#
#     bpt_x[sngood != 0] = np.log10((f_n2[sngood != 0]) / (f_ha[sngood != 0]))
#     bpt_y[sngood != 0] = np.log10((f_o3[sngood != 0]) / (f_hb[sngood != 0]))
#
#     k03 = lambda bpt_x: 0.61 / (bpt_x - 0.05) + 1.3
#     kewley = lambda bpt_x: 0.61 / (bpt_x - 0.47) + 1.19
#     kewley_bpt_y = kewley(bpt_x)
#     yes_agn = np.zeros(f_ha.shape)
#     yes_agn[bpt_x >= 0.47] = 1
#     yes_agn[(kewley_bpt_y <= bpt_y) & (sngood != 0)] = 1
#     # yes_agn[sngood == 0] = 1.
#     return yes_agn
#


def bpt_stacked(galname):
    f_ha = get_array_reproj_stacked(galname, 'Ha6563', 2)
    f_hb = get_array_reproj_stacked(galname, 'Hb4861', 2)
    f_o3 = get_array_reproj_stacked(galname, 'OIII5007', 2)
    f_n2 = get_array_reproj_stacked(galname, 'NII6583', 2)
    ew_ha = get_array_reproj_stacked(galname, 'Ha6563', 3)

    f_ha_err = get_array_reproj_stacked(galname, 'Ha6563', 6)
    f_hb_err = get_array_reproj_stacked(galname, 'Hb4861', 6)
    f_o3_err = get_array_reproj_stacked(galname, 'OIII5007', 6)
    f_n2_err = get_array_reproj_stacked(galname, 'NII6583', 6)

    sn_cut = 3.0
    good = (f_n2 > 0) & (f_o3 > 0) & (f_ha > 0) & (f_hb > 0) & (
        ~np.isnan(ew_ha)) & (np.abs(f_n2 / f_n2_err) > sn_cut) & (np.abs(
            f_o3 / f_o3_err) > sn_cut) & (np.abs(
                f_ha / f_ha_err) > sn_cut) & (np.abs(f_hb / f_hb_err) > sn_cut)
    n2ha = np.full(f_n2.shape, np.nan)
    n2ha[good] = np.log10(f_n2[good]) - np.log10(f_ha[good])
    o3hb = np.full(f_o3.shape, np.nan)
    o3hb[good] = np.log10(f_o3[good]) - np.log10(f_hb[good])

    # Eq. 5 of 2001ApJ...556..121K
    kewley01 = lambda nii: 1.19 + 0.61 / (nii - 0.47)
    # Eq. 1 of 2003MNRAS.346.1055K
    kauffm03 = lambda nii: 1.30 + 0.61 / (nii - 0.05)
    # Eq. 3 of 2010MNRAS.403.1036C
    cidfer10 = lambda nii: 0.48 + 1.01 * nii

    bpt = np.full(n2ha.shape, np.nan)
    # Star forming: below Kauffmann line and EW > 6
    sf = (n2ha > -1.5) & (n2ha < -0.1) & (o3hb < kauffm03(n2ha)) & (abs(ew_ha)
                                                                    > 6.0)
    bpt[sf == True] = -1
    # Intermediate: below Kewley line and not star-forming
    inter = (~sf) & (n2ha < 0.3) & (o3hb < kewley01(n2ha))
    bpt[inter == True] = 0
    # LINER: above Kewley line and below Cid Fernandes line
    liner = (~sf) & (~inter) & (o3hb > -1) & (o3hb < cidfer10(n2ha))
    bpt[liner == True] = 1
    # Seyfert: above Kewley line and above Cid Fernandes line
    seyfert = (~sf) & (~inter) & (~liner) & (o3hb > -1)
    bpt[seyfert == True] = 2
    return bpt


# def calculate_o3n2(hb, o3, ha, n2, mask_zones = None, tau_V = None, correct = False):
#     if mask_zones is not None:
#         mask = mask_zones
#     else:
#         mask = np.zeros_like(Hb_obs, dtype = np.bool_)
#     Hb = np.ma.masked_array(Hb_obs, mask = mask)
#     O3 = np.ma.masked_array(O3_obs, mask = mask)
#     Ha = np.ma.masked_array(Ha_obs, mask = mask)
#     N2 = np.ma.masked_array(N2_obs, mask = mask)
#     if correct is True:
#         tau_V_m = np.ma.masked_array(tau_V, mask = mask)
#         from pystarlight.util import redenninglaws
#         q = redenninglaws.Cardelli_RedLaw([4861, 5007, 6563, 6583])
#         Hb *= np.ma.exp(q[0] * tau_V_m)
#         O3 *= np.ma.exp(q[1] * tau_V_m)
#         Ha *= np.ma.exp(q[2] * tau_V_m)
#         N2 *= np.ma.exp(q[3] * tau_V_m)
#     O3Hb = np.ma.log10(O3/Hb)
#     N2Ha = np.ma.log10(N2/Ha)
#     O3N2 = np.ma.log10(O3 * Ha / (N2 * Hb))
#     return O3, Hb, N2, Ha, O3Hb, N2Ha, O3N2


def metallicity_stacked_o3n2(galname):
    '''
    ** Not used anymore. See the [NII] function. **
    12+log(O/H) = 8.533 - 0.214 * O3N2,
    where fit parameters are from R.A. Marino+2013 (M13),
    and O3N2 = log10( [OIII]5007 / Hb * Ha / [NII]6583 ),
    first defined in Alloin+79.
    '''
    f_ha = get_array_reproj_stacked(galname, 'Ha6563', 2)
    f_hb = get_array_reproj_stacked(galname, 'Hb4861', 2)
    f_o3 = get_array_reproj_stacked(galname, 'OIII5007', 2)
    f_n2 = get_array_reproj_stacked(galname, 'NII6583', 2)
    ew_ha = get_array_reproj_stacked(galname, 'Ha6563', 3)

    f_ha_err = get_array_reproj_stacked(galname, 'Ha6563', 6)
    f_hb_err = get_array_reproj_stacked(galname, 'Hb4861', 6)
    f_o3_err = get_array_reproj_stacked(galname, 'OIII5007', 6)
    f_n2_err = get_array_reproj_stacked(galname, 'NII6583', 6)

    good = (f_n2 > 0) & (f_o3 > 0) & (f_ha > 0) & (f_hb >
                                                   0) & (~np.isnan(ew_ha))
    n2ha = np.full(f_n2.shape, np.nan)
    n2ha[good] = np.log10(f_n2[good]) - np.log10(f_ha[good])
    o3hb = np.full(f_o3.shape, np.nan)
    o3hb[good] = np.log10(f_o3[good]) - np.log10(f_hb[good])

    bpt = bpt_stacked(galname)

    o3n2 = np.log10(o3hb / n2ha)
    o3n2[bpt != -1] = np.nan

    m13 = 8.533 - 0.214 * o3n2

    res = dict()
    res['metallicity'] = m13
    out_name = fname_metallicity_stacked(galname)
    with open(out_name, 'wb') as p:
        pickle.dump(res, p)


def metallicity_stacked_n2(galname, save=True):
    '''
    12+log(O/H) = 9.12 + 0.73 * N2,
    where fit parameters are from Denicolo+02,
    and N2 = log10( [NII]6584 / Ha )
    '''
    f_ha = get_array_reproj_stacked(galname, 'Ha6563', 2)
    # f_hb = get_array_reproj_stacked(galname, 'Hb4861', 2)
    # f_o3 = get_array_reproj_stacked(galname, 'OIII5007', 2)
    f_n2 = get_array_reproj_stacked(galname, 'NII6583', 2)
    ew_ha = get_array_reproj_stacked(galname, 'Ha6563', 3)

    f_ha_err = get_array_reproj_stacked(galname, 'Ha6563', 6)
    # f_hb_err = get_array_reproj_stacked(galname, 'Hb4861', 6)
    # f_o3_err = get_array_reproj_stacked(galname, 'OIII5007', 6)
    f_n2_err = get_array_reproj_stacked(galname, 'NII6583', 6)

    good = (f_n2 > 0) & (f_ha > 0) & (~np.isnan(ew_ha))

    n2ha = np.full(f_n2.shape, np.nan)
    n2ha[good] = np.log10(f_n2[good]) - np.log10(f_ha[good])

    bpt = bpt_stacked(galname)

    # o3n2 = np.log10(o3hb / n2ha)
    n2ha[bpt != -1] = np.nan

    g12 = 9.12 + 0.73 * n2ha

    res = dict()
    res['metallicity'] = g12
    if save == True:
        out_name = fname_metallicity_stacked(galname)
        with open(out_name, 'wb') as p:
            pickle.dump(res, p)
    else:
        return g12


def alpha_co_metallicity(galname):
    '''
    12+log(O/H) = 9.12 + 0.73 * N2,
    where fit parameters are from Denicolo+02,
    and N2 = log10( [NII]6584 / Ha )
    '''
    g12 = metallicity_stacked_n2(galname, save=False)

    alpha_co = 10**(12. - 1.3 * g12)
    # Get rid of helium (pretty sure it is included in the Genzel paper)
    alpha_co /= 1.36
    alpha_co[np.isnan(g12)] = np.nan

    res = dict()
    res['alpha_co'] = alpha_co

    out_name = fname_alpha_co_stacked(galname)
    with open(out_name, 'wb') as p:
        pickle.dump(res, p)


def reproj_sfr_stacked(galname):
    """
    SFR surface density map from stacked spectral fitting.

    Returns:
        array of SFR surface density (Msun/yr/kpc^2),
        inclination-deprojected, in CO projection (6 arcsec pixels).
    """
    if len(
            glob.glob(
                '/Users/ryan/venus/shared_data/califa/DR3-stack/%s/%s_result.pk'
                % (galname, galname))) == 0:
        print("Nothing exists for this galaxy")
        return 0

    # chi2 = get_array_reproj_stacked(galname, 'chi2')
    # 2, 6
    bpt = bpt_stacked(galname)

    halpha_map = get_array_reproj_stacked(galname, 'Ha6563', 2)
    halpha_err = get_array_reproj_stacked(galname, 'Ha6563', 6)

    hbeta_map = get_array_reproj_stacked(galname, 'Hb4861', 2)
    hbeta_err = get_array_reproj_stacked(galname, 'Hb4861', 6)

    ext_corr = extcorr(halpha_map, hbeta_map)
    ext_corr[(hbeta_map / hbeta_err < 2) |
             (halpha_map / halpha_err < 2)] = np.nan

    halpha_ew_map = get_array_reproj_stacked(galname, 'Ha6563', 3)
    halpha_ew_err_map = get_array_reproj_stacked(galname, 'Ha6563', 7)
    ext_corr[bpt != -1] = np.nan
    # ext_corr[halpha_ew_map < 6] = np.nan

    pp = 10**(0.4 * 5.86 * np.log10(halpha_map / hbeta_map / 2.86))
    sigtemp = 2.303 * pp * 0.4 * 5.86 * 0.434 * np.sqrt(
        (halpha_err / halpha_map)**2 + (hbeta_err / hbeta_map)**2)

    halpha_err_corr = np.abs(
        halpha_map * pp) * np.sqrt((halpha_err / halpha_map)**2 +
                                   (sigtemp / pp)**2)

    halpha_map *= ext_corr
    # halpha_err *= ext_corr

    gal_dict_i = gal_dict[galname]
    d_mpc = gal_dict_i['dist_Mpc'] * u.Mpc
    d_cm = np.float64((d_mpc).to(u.cm).value)
    l_ha_to_sfr = 5.3e-42
    f_ha_to_erg_s_cm2 = 1e-16
    pc2_to_kpc2 = 1e6
    pix_size_arcsec = 6.
    d_kpc = np.float64(d_mpc.to(u.kpc).value)
    pix_area_kpc2 = (pix_size_arcsec * d_kpc / 206265.)**2
    cos_i = np.cos(np.radians(gal_dict_i['incl']))

    units_factor = cos_i * 4. * np.pi * d_cm**2 * l_ha_to_sfr * f_ha_to_erg_s_cm2 / pix_area_kpc2

    res = dict()
    # res['halpha_flux'] = halpha * units_factor
    # res['halpha_flux_err'] = halpha_err * units_factor
    sigma_sfr_err = halpha_err_corr * units_factor
    sigma_sfr = halpha_map * units_factor
    sigma_sfr_err[sigma_sfr > 1e7 * 7.9 / 5.3 * pix_area_kpc2 * 1e-6 /
                  cos_i] = np.nan
    sigma_sfr[sigma_sfr > 1e7 * 7.9 / 5.3 * pix_area_kpc2 * 1e-6 /
              cos_i] = np.nan

    res['sigma_sfr'] = sigma_sfr
    res['sigma_sfr_err'] = sigma_sfr_err

    out_name = fname_halpha_stacked(galname)

    with open(out_name, 'wb') as p:
        pickle.dump(res, p)


def loop_over_galaxies(fun):
    fname_surf_dens_all_dict = '/Users/ryan/Dropbox/mac/wise_w3_vs_co/surf_dens_all.pk'
    surf_dens_all = pickle.load(open(fname_surf_dens_all_dict, 'rb'))
    galnames_all = list(surf_dens_all.keys())

    for galname in galnames_all:
        if len(
                glob.glob(
                    '/Users/ryan/venus/shared_data/califa/DR3-stack/%s/%s_result.pk'
                    % (galname, galname))) == 0:
            print("Nothing exists for this galaxy")
            continue

        fun(galname)


if __name__ == '__main__':

    if sys.argv[1] == 'halpha':
        loop_over_galaxies(reproj_halpha)

    if sys.argv[1] == 'mstar':
        loop_over_galaxies(reproj_mstar)

    if sys.argv[1] == 'sfr_stacked':
        loop_over_galaxies(reproj_sfr_stacked)

    if sys.argv[1] == 'mstar_stacked':
        loop_over_galaxies(reproj_mstar_stacked)

    if sys.argv[1] == 'metallicity_stacked':
        loop_over_galaxies(metallicity_stacked_n2)

    if sys.argv[1] == 'alpha_co':
        loop_over_galaxies(alpha_co_metallicity)

    if sys.argv[1] == 'pixel_mapping':
        # fname_surf_dens_all_dict = '/Users/ryan/Dropbox/mac/wise_w3_vs_co/surf_dens_all.pk'
        fname_surf_dens_all_dict = '/Users/ryan/Dropbox/mac/wise_w3_vs_co/surf_dens_all_v2.pk'  # Without Aniano kernels, until that is fixed, and with more conservative s/n cuts in mom0 map making
        surf_dens_all = pickle.load(open(fname_surf_dens_all_dict, 'rb'))
        galnames_all = list(surf_dens_all.keys())
        res = dict()
        i = 0
        n = len(galnames_all)
        for galname in galnames_all:
            print(galname)
            print("%i / %i" % (i, n))
            rr = reproj_califa_mapping(galname)
            if rr != 0:
                res[galname] = rr
            i += 1

        with open('/Users/ryan/Dropbox/mac/wise_w3_vs_co/pixel_mapping.pk',
                  'wb') as p:
            pickle.dump(res, p)
