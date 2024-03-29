import copy
import itertools
import os

from astropy import wcs
from astropy.io import fits
from tqdm import tqdm
import astropy.units as u
import numpy as np
import scipy.interpolate as si
import scipy.spatial as sspatial


from .plot import PlotPointing
from .utils import SpiceFilename


def _get_spatial_wcs(header):
    """ Get subset WCS for the spatial coordinates of a SPICE HDU.

    Parameters
    ==========
    header : fits.Header
        SPICE L2 FITS header.

    Returns
    =======
    w : wcs.WCS
        Spatial WCSa.
    """
    w = wcs.WCS(header)
    w.wcs.pc[3, 0] = 0  # remove PC4_1 to decouple t from x
    return w.sub(2)  # drop wavelength and time axes


def _assert_hdus_have_the_same_spatial_coordinates(hdul):
    """ Assert that SPICE image HDUs have the same spatial coordinates

    Parameters
    ==========
    hdul : fits.HDUList
        HDU list with images extensions.

    Raises
    ======
    ValueError :
        When the spatial coordinates of the HDUs differ.

    """
    w_ref = _get_spatial_wcs(hdul[0].header)
    for hdu in hdul:
        if hdu.is_image and hdu.name != 'WCSDVARR':
            w = _get_spatial_wcs(hdu.header)
            if w.to_header() != w_ref.to_header():
                print(w_ref)
                print(w)
                raise ValueError('FITS extensions have different spatial '
                                 'coordinates')


def get_coordinates(header):
    """ Get SPICE coordinates (without distortion)
    
    Parameters
    ==========
    header : fits.Header
        Header of the FITS extension for which to compute the coordinates.

    Returns
    =======
    coordinates : np.ndarray
        Array of shape (2, ny, nx), containing the helioprojective coordinates
        `T_x` and `T_y` for each pixel of the field of view.
    """
    w = wcs.WCS(header)
    iy, ix = np.indices((header['NAXIS2'], header['NAXIS1']))
    iD = np.zeros_like(ix)  # dummy D indices
    it = np.zeros_like(ix)  # dummy t indices
    Tx, Ty, _, _ = w.pixel_to_world(ix, iy, iD, it)
    Tx = Tx.to('arcsec').value
    Ty = Ty.to('arcsec').value
    pi = u.Quantity(np.pi, 'rad').to('arcsec').value
    Tx = (Tx + pi) % (2*pi) - pi
    Ty = (Ty + pi) % (2*pi) - pi
    return np.array([Tx, Ty])


def add_distortion_to_coordinates(coordinates, hdul):
    """ Add distortion to coordinates

    Parameters
    ==========
    coordinates : array_like
        Array of shape (2, ny, nx), containing the helioprojective coordinates
        `T_x` and `T_y` for each pixel of the field of view.
    hdul : fits.HDUList
        HDU list with two extensions images ``WCSDVARR`` (1 and 2) containing
        the pointing distortion correction for `T_x` and `T_y`.

    Returns
    =======
    coordinates : np.ndarray
        Array of shape (2, ny, nx), containing the corrected coordinates.
    """
    Tx, Ty = coordinates
    Tx_corr = hdul['WCSDVARR', 1].data
    Ty_corr = hdul['WCSDVARR', 2].data
    return np.array([Tx - Tx_corr, Ty - Ty_corr])


def get_interpolation_points(hdul):
    """ Get the origin and target points for the interpolation

    Parameters
    ==========
    hdul : fits.HDUList
        HDU list with the image extensions for which to compute the
        coordinates, plus two image extensions ``WCSDVARR`` (1 and 2)
        containing the pointing distortion correction for `T_x` and `T_y`.

    Returns
    =======
    points : array_like
        Array of shape (ny*nx, 2) containing the coordinates of the data
        points.
    new_points : array_like
        Array of shape (ny*nx, 2) containing the points at which to interpolate
        the data.
    """
    _assert_hdus_have_the_same_spatial_coordinates(hdul)

    coord_grid = get_coordinates(hdul[0].header)
    coord_real = add_distortion_to_coordinates(coord_grid, hdul)

    # Interpolation initial and target points
    points = coord_real.reshape(2, -1).T
    new_points = np.moveaxis(coord_grid, 0, -1)
    return points, new_points


def remap_spice_hdu(hdu, points, new_points, sum_wvl=False):
    """ Remap a SPICE spectral cube to corrected coordinates

    Parameters
    ==========
    hdu : fits.PrimaryHDU, fits.ImageHDU, fits.BinTableHDU
        SPICE L2 FITS HDU to remap. (If the HDU is not of 'image' type, return
        it without modification.)
    points : array_like or sspatial.Delaunay
        Array of shape (ny*nx, 2) containing the coordinates of the data
        points, or precomputed Delaunay triangulation.
    new_points : array_like
        Array of shape (ny*nx, 2) containing the points at which to interpolate
        the data.
    new_points : array_like
    sum_wvl : bool
        If True, sum along wavelength axis to generate a quicklook image.

    Returns
    =======
    hdu : fits.PrimaryHDU, fits.ImageHDU, fits.BinTableHDU
        Aligned SPICE 'L2r' HDU, matching the type of the input HDU.
    """
    if not hdu.is_image or hdu.name == 'WCSDVARR':
        return hdu

    new_hdu = hdu.copy()
    if sum_wvl:
        # Integrated intensity
        img = np.nansum(hdu.data, axis=1)  # Sum over wavelengths
        img = np.squeeze(img)  # Collapse 1-depth axis (t or X)
        interp = si.LinearNDInterpolator(points, img.flatten())
        new_img = interp(new_points)
        new_hdu.data = new_img.reshape(1, 1, *new_img.shape)
    else:
        # Remap to new coordinates within time and/or wvl slices
        nt, nD, _, _ = hdu.data.shape
        itD = itertools.product(range(nt), range(nD))
        for it, iD in tqdm(itD, desc=f'Remapping {hdu.name}', total=nt * nD):
            img = hdu.data[it, iD]
            interp = si.LinearNDInterpolator(points, img.flatten())
            new_img = interp(new_points)
            new_hdu.data[it, iD] = new_img

    new_hdu.update_header()
    new_hdu.header.add_history('jitter_correction.py')
    new_hdu.add_datasum()
    new_hdu.add_checksum()
    return new_hdu


def correct_jitter(
        filename, output_dir, overwrite=False,
        plot_results=False, sum_wvl=False, windows=None,
        ):
    """ Correct the pointing in a SPICE level 2 FITS

    Parameters
    ==========
    filename : str
        Path to SPICE L2 FITS.
    output_dir : str
        Directory where the corrected FITS and plots are saved.
    overwrite : bool
        Overwrite output FITS if it already exists.
    plot_results : bool
        Generate plots to visualize the results.
    sum_wvl : bool
        If True, sum along wavelength axis to generate a quicklook image.
    windows : list of str or None (default: None)
        Windows to keep. If None, keep all windows.

    The aligned fits are saved in <output_dir> under the name
    <solo_L2r_spice-....fits> when sum_wvl is False, or
    <solo_L2r_quicklook_spice-....fits> when sum_wvl is True.

    Returns
    =======
    output_fits : str
        Path to the saved L2r FITS.
    """

    os.makedirs(output_dir, exist_ok=True)
    # filename operations
    filename = SpiceFilename(filename)
    filename_output = copy.copy(filename)
    filename_output['path'] = output_dir
    if sum_wvl:
        filename_output['level'] = 'L2r_quicklook'
    else:
        filename_output['level'] = 'L2r'

    if os.path.isfile(filename_output.to_str()) and not overwrite:
        print(f'Aligned file exists: {filename_output.to_str()}, exiting')
        return filename_output.to_str()

    # open FITS
    hdulist = fits.open(filename.to_str())
    points, new_points = get_interpolation_points(hdulist)
    points_trig = sspatial.Delaunay(points)

    # interpolate data
    new_hdulist = fits.HDUList(hdus=[])
    if windows is None:
        windows = [hdu.name for hdu in hdulist
                   if hdu.is_image and (hdu.name != 'WCSDVARR')]
    windows = windows + [hdu.name for hdu in hdulist
                         if not hdu.is_image or (hdu.name == 'WCSDVARR')]
    for win in windows:
        hdu = hdulist[win]
        new_hdu = remap_spice_hdu(
            hdu, points_trig, new_points, sum_wvl=sum_wvl)
        new_hdulist.append(new_hdu)

    # save data
    new_hdulist.writeto(filename_output.to_str(), overwrite=True)

    # generate plots
    if plot_results:
        p = PlotPointing()
        spice_name = os.path.basename(hdulist.filename())
        p.plot_hdulist(
            hdulist, spice_name,
            filename_output.to_str(ext='_plot_not_aligned.pdf')
            )
        p.plot_hdulist(
            new_hdulist, spice_name,
            filename_output.to_str(ext='_plot_aligned.pdf')
            )

    return filename_output.to_str()
