import math
import os

from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits, votable
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
import numpy as np
import numpy.core.records as rec
from astropy.table import Table

_allowed_weights = ['square', 'linear', 'none']


class IslandRange(object):
    def __init__(self, isle_id):
        self.isle_id = isle_id



def point_in_ellipse(origin, point, a, b, pa_rad, verbose=False):
    """
    Identify if the point is inside the ellipse.

    :param origin A SkyCoord defining the centre of the ellipse.
    :param point A SkyCoord defining the point to be checked.
    :param a The semi-major axis in arcsec of the ellipse
    :param b The semi-minor axis in arcsec of the ellipse
    :param pa_rad The position angle of the ellipse. This is the angle of the major axis measured in radians East of 
                   North (or CCW from the y axis).
    """
    # Convert point to be in plane of the ellipse, accounting for distortions at high declinations
    p_ra_dist = (point.icrs.ra.degree - origin.icrs.ra.degree)* math.cos(origin.icrs.dec.rad)
    p_dec_dist = point.icrs.dec.degree - origin.icrs.dec.degree

    # Calculate the angle and radius of the test opoint relative to the centre of the ellipse
    # Note that we reverse the ra direction to reflect the CCW direction
    radius = math.sqrt(p_ra_dist**2 + p_dec_dist**2)
    diff_angle = (math.pi/2 + pa_rad) if p_dec_dist == 0 else math.atan(p_ra_dist / p_dec_dist) - pa_rad

    # Obtain the point position in terms of the ellipse major and minor axes
    minor = radius * math.sin(diff_angle)
    major = radius * math.cos(diff_angle)
    if verbose:
        print ('point relative to ellipse centre angle:{} deg radius:{:.4f}" maj:{:.2f}" min:{:.2f}"'.format(math.degrees(diff_angle), radius*3600, 
                major*3600, minor*3600))
    
    a_deg = a / 3600.0
    b_deg = b / 3600.0

    # Calc distance from origin relative to a and b
    dist = math.sqrt((major / a_deg) ** 2 + (minor / b_deg) ** 2)
    if verbose:
        print("Point %s is %f from ellipse %f, %f, %f at %s." % (point, dist, a, b, math.degrees(pa_rad), origin))
    return round(dist,3) <= 1.0


def get_weighting_array(data, velocities, continuum_start_vel, continuum_end_vel, weighting='square'):
    """
    Calculate the mean of the continuum values. This is based on precalculated regions where there is no gas expected.
    :param data: A cubelet to be analysed, should be a 3D array of flux values.
    :param velocities: A numpy array of velocity values in m/s
    :param continuum_start_vel: The lower bound of the continuum velocity range (in m/s)
    :param continuum_end_vel: The upper bound of the continuum velocity range (in m/s)
    :param weighting: The weighting scheme to use, one of square, linear, none
    :return: A 2D array of weighting values for the
    """

    if (continuum_start_vel > np.max(velocities)) or (continuum_end_vel < np.min(velocities)):
        raise Exception("Continuum range {} to {} is outside of the data velocity range {} to {}".format(
            continuum_start_vel, continuum_end_vel, np.min(velocities), np.max(velocities)))

    continuum_range = np.where(continuum_start_vel < velocities)
    if len(continuum_range) == 0:
        return np.zeros(data.shape[1:2])

    if weighting not in _allowed_weights:
        raise Exception("Weighting must by one of ", ', '.join(_allowed_weights))

    bin_start = continuum_range[0][0]
    continuum_range = np.where(velocities < continuum_end_vel)
    bin_end = continuum_range[0][-1]

    # print("Using bins %d to %d (velocity range %d to %d) out of %d" % (
    #    bin_start, bin_end, continuum_start_vel, continuum_end_vel, len(velocities)))
    # print(data.shape)
    continuum_sample = np.array(data[bin_start:bin_end, :, :])
    continuum_sample[continuum_sample<0]=0

    # print ("...gave sample of", continuum_sample)
    mean_cont = np.nanmean(continuum_sample, axis=0)
    mean_cont[mean_cont<0]=0
    if weighting == 'square':
        mean_sq = mean_cont ** 2
        sum_sq = np.nansum(mean_sq)
        weights = mean_sq / sum_sq
    elif weighting == 'linear':
        weights = mean_cont / np.nansum(mean_cont)
    else:
        # No weights, just trim to ellipse
        weights = mean_cont
        weights[weights>0]=1
        weights = weights/np.sum(weights)
    # print("Got weighting of {} from {} and {}".format(weighting, mean_sq, sum_sq))
    return weights


def get_integrated_spectrum(image, w, src, velocities, continuum_start_vel, continuum_end_vel, radius=None, 
                            plot_weight_path=None, weighting='square', target_name=None):
    """
    Calculate the integrated spectrum of the component.
    :param image: The image's data array
    :param w: The image's world coordinate system definition
    :param src: The details of the component (or if an array, components) being processed, must have ra, dec, a, b, pa and comp_name keys
    :param velocities: A numpy array of velocity values in m/s
    :param continuum_start_vel: The lower bound of the continuum velocity range (in m/s)
    :param continuum_end_vel: The upper bound of the continuum velocity range (in m/s)
    :param radius: The radius of the box around the source centre where data will be checked for membership of the 
        source ellipse. Default is to use the semi-major axis of the source.
    :param plot_weight_path: The path to which diagnostic plots are output. Default is not to output plots.
    :param weighting: The weighting scheme to use, one of square, linear, none
    :return: An array of average flux/pixel across the component at each velocity step
    """

    if weighting not in _allowed_weights:
        raise Exception("Weighting must by one of ", ', '.join(_allowed_weights))

    all_sources = Table(names=('comp_name', 'ra', 'dec', 'a', 'b', 'pa'), dtype=(np.str_, np.float, np.float, np.float, np.float, np.float))
    if isinstance(src, list):
        for row in src:
            all_sources.add_row((row['comp_name'], row['ra'], row['dec'], row['a'], row['b'], row['pa']))
            #all_sources = Table(names=('comp_name', 'ra', 'dec', 'a', 'b', 'pa'), dtype=(np.str_, np.float, np.float, np.float, np.float, np.float),
            #    data=(src['comp_name'], src['ra'], src['dec'], src['a'], src['b'], src['pa']))
    else:
        print (src)
        all_sources.add_row((src['comp_name'], src['ra'], src['dec'], src['a'], src['b'], src['pa']))

    if plot_weight_path:
        if isinstance(src, list):
            print ("getting spectrum for sources " + ', '.join([item['comp_name'] for item in src]))
        else:
            print ("getting spectrum for source " + str(src))

    if target_name is None:
        target_name = all_sources[0]['comp_name']

    has_stokes = len(image.shape) > 3
    #pix = w.wcs_world2pix(all_sources['ra'], all_sources['dec'], 0, 0, 1) if has_stokes else w.wcs_world2pix(all_sources['ra'], all_sources['dec'], 0, 1)
    pix = w.wcs_world2pix(all_sources['ra'], all_sources['dec'], 0, 0, 1) if has_stokes else w.wcs_world2pix(all_sources['ra'], all_sources['dec'], 0, 1)
    x_coord = np.round(pix[0]).astype(int) - 1  # 266
    y_coord = np.round(pix[1]).astype(int) - 1  # 197
    #if not radius:
    #    radius = np.ceil(all_sources['a'])
    #print("Translated %.4f, %.4f to %d, %d" % (
    #    src['ra'], src['dec'], x_coord, y_coord))
    #print (w)
    y_min = 0
    y_max = image.shape[-2]-1
    x_min = 0
    x_max = image.shape[-1]-1
    #data = np.copy(image[0, :, y_min:y_max+1, x_min:x_max+1]) if has_stokes else np.copy(image[:, y_min:y_max+1, x_min:x_max+1])
    data = np.copy(image[0, :, :, :]) if has_stokes else np.copy(image[:, :, :])
    if plot_weight_path:
        # non wcs plot
        fig, ax = plt.subplots(1, 1, figsize=(9, 3))
        summed_data = np.nansum(data, axis=0)
        ax.imshow(summed_data, origin='lower')
        plt.title(target_name) 
        fname = plot_weight_path + '/'+ target_name + '_data.png'
        print ('Plotting data to ' + fname) 
        plt.savefig(fname, bbox_inches='tight')
        plt.close()

        # wcs plot
        plt.subplot(projection=w.celestial)
        center_chan = data[data.shape[0] // 2,:,:]
        plt.imshow(center_chan, origin='lower')
        #if has_stokes:
        #    plt.imshow(image[0,10,:,:], origin='lower')
        #else:
        #    plt.imshow(image[10,:,:], origin='lower')
        plt.grid(color='white', ls='solid')
        fname = plot_weight_path + '/'+ target_name + '_image_wcs.png'
        print ('Plotting image wcs to ' + fname) 
        plt.savefig(fname, bbox_inches='tight')
        plt.close()


    origin = np.asarray(SkyCoord(all_sources['ra'], all_sources['dec'], frame='icrs', unit="deg"))
    pa_rad = np.radians(all_sources['pa'])
    total_pixels = (y_max-y_min +1) * (x_max-x_min +1)
    #total_pixels = data.shape[1] * data.shape[2]
    outside_pixels = 0
    for x in range(x_min, x_max+1):
        for y in range(y_min, y_max+1):
            eq_pos = w.wcs_pix2world(x, y, 0, 0, 0) if has_stokes else w.wcs_pix2world(x, y, 0, 0)
            point = SkyCoord(eq_pos[0], eq_pos[1], frame='icrs', unit="deg")
            in_ellipse = False
            for idx, a_src in enumerate(all_sources):
                in_ellipse |= point_in_ellipse(origin[idx], point, a_src['a'], a_src['b'], pa_rad[idx])
            if not in_ellipse:
                data[:, y-y_min, x-x_min] = 0
                outside_pixels += 1
            #print (point.ra, point.dec, x, y, in_ellipse)

    # print("Found {} pixels out of {} inside the component {} at {} {}".format(total_pixels - outside_pixels, total_pixels,
    #                                                                   src['comp_name'],
    #                                                                   point.galactic.l.degree,
    #                                                                   point.galactic.b.degree))
    weights = get_weighting_array(data, velocities, continuum_start_vel, continuum_end_vel, weighting=weighting)
    integrated = np.nansum(data * weights, axis=(1, 2))
    inside_pixels = total_pixels - outside_pixels
    if inside_pixels <= 0:
        print ("Error: No data for component!")
    else:
        integrated /= inside_pixels

    if plot_weight_path:
        fig, ax = plt.subplots(1, 1, figsize=(9, 3))
        pos = ax.imshow(weights, origin='lower')
        fig.colorbar(pos, ax=ax)
        plt.title(target_name)
        fname = plot_weight_path + '/'+ target_name + '_weights.png'
        print ('Plotting weights to ' + fname) 
        for idx, a_src in enumerate(all_sources):
            print ('Ellipse ra={} dec={} pa={:.03f} deg {:.03f}pi rad'.format(a_src['ra'], a_src['dec'], a_src['pa'], pa_rad[idx]/math.pi))
        plt.savefig(fname, bbox_inches='tight')
        plt.close()

    return integrated


def find_edges(fluxes, num_edge_chan):
    """
    Seek from the edges to find where the data starts for this set of fluxes.
    This accounts for an optional number of channels in the data which have no
    data recorded.
    :param fluxes: The array of fluxes to be checked.
    :param num_edge_chan: The number of edge channels with data to be skipped
    :return: The index of the first and last cell to have data.
    """

    l_edge = 0
    r_edge = len(fluxes)-1

    while l_edge < len(fluxes)-1 and fluxes[l_edge] == 0:
        l_edge += 1

    while r_edge > 0 and fluxes[r_edge] == 0:
        r_edge -= 1

    return l_edge + num_edge_chan, r_edge - num_edge_chan
