# Python script to produce spectra for a series of sub-cubes.
# The spectra are assumed to be absoprtion spectra and rated according to their suitability for absorption detection.
# Emission spectra extraction is based on work by Claire Murray

# Author James Dempsey
# Date 19 Feb 2020

import argparse
import csv 
import glob
import itertools
import math
import os
import time
import warnings
import zipfile

import astropy.units as u
from astropy.coordinates import SkyCoord, FK4
from astropy.io import fits, votable
from astropy.io.ascii import Csv
from astropy.io.votable.tree import Param,Info
from astropy.io.votable import from_table, parse_single_table, writeto
from astropy.table import QTable, Table, Column
from astropy.utils.exceptions import AstropyWarning
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors
import numpy as np
import numpy.core.records as rec

from spectral_cube import SpectralCube
import radio_beam
import aplpy


from RadioAbsTools import cube_tools, spectrum_tools


class AbsRun:
    def __init__(self):
        self.length = 0
        self.start_idx = 0
        self.start_vel = 0
        self.end_vel = 0
        self.max_sigma = 0
    
    def __init__(self, length, start_idx, start_vel, end_vel, max_sigma):
        self.length = length
        self.start_idx = start_idx
        self.start_vel = start_vel
        self.end_vel = end_vel
        self.max_sigma = max_sigma

    def __str__(self):
        return "AbsRun: " + str(self.__dict__)
    
    def __repr__(self):
        return "<AbsRun len:%s start:%s>" % (self.length, self.start_idx)



def parseargs():
    """
    Parse the command line arguments
    :return: An args map with the parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Produce absorption spectra from a set of subcubes")
    parser.add_argument("-s", "--sbid", help="The id of the ASKAP scheduling block to be processed",
                        type=int, required=True)
    parser.add_argument("-c", "--catalogue", help="The catalgoue of source positions and characteristics",
                        required=True)
    parser.add_argument("-p", "--parent", help="The parent folder for the processing, will default to sbnnn/ where nnn is the sbid.",
                        required=False)
    #parser.add_argument("-e", "--emission", help="The HI emission cube to source emission data around each source.",
    #                    required=False)
    parser.add_argument("--skip_abs", help="Use the existing absorption spectra",
                        action='store_true', required=False)
    parser.add_argument("--continuum_start", help="The lowest velocity (in km/s) of the conitnuum region.", type=int, 
                        default=-100, required=False)
    parser.add_argument("--continuum_end", help="The lowest velocity (in km/s) of the continuum region.", type=int, 
                        default=-60, required=False)
    parser.add_argument("--no_zoom", help="Don't zoom in the combined spectra",
                        action='store_true', required=False)

    args = parser.parse_args()
    return args


def read_targets(targets_file):
    ids=[]
    comp_names=[]
    ras=[]
    decs=[]
    beams=[]
    i=1
    with open(targets_file, 'r') as csvfile:
        tgt_reader = csv.reader(csvfile, delimiter=',',
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for row in tgt_reader:
            if (tgt_reader.line_num == 1):
                # skip header
                continue

            #print (row)
            ids.append(i)
            comp_names.append(row[1])
            ras.append(float(row[2]))
            decs.append(float(row[3]))
            beams.append(row[4:])

            i+=1

            
    table = Table()
    table.add_column(Column(name='id', data=ids))
    table.add_column(Column(name='comp_name', data=comp_names))
    table.add_column(Column(name='ra', data=ras))
    table.add_column(Column(name='dec', data=decs))
    table.add_column(Column(name='beams', data=beams))

    return table


def rename_columns(table):
    names = np.asarray(table.colnames)
    for name in names:
        if name.startswith('col_'):
            table.rename_column(name, name[4:])

def plot_mom0(fname, comp_name, out_folder, src_ra, src_dec, src_maj, src_min, src_pa):
    cube = SpectralCube.read(fname)
    cube.beam_threshold = 0.13
    m0 = cube.moment0().to(u.Jy*u.km/(u.beam*u.s))
    m0.write('moment-0.fits', overwrite=True)

    fig = aplpy.FITSFigure('moment-0.fits')
    fig.show_grayscale()#pmax=100)
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text('Brightness (Jy km / s / beam)')
    fig.add_beam()
    fig.set_theme('publication')

    # Plot ellipse for each source
    a_deg = src_maj*u.arcsec.to(u.deg)*2
    b_deg = src_min*u.arcsec.to(u.deg)*2
    fig.show_ellipses(src_ra, src_dec, a_deg, b_deg, angle=src_pa-90, edgecolor='red') #, 

    figname = '{}/{}_mom0'.format(out_folder, comp_name)
    #fig.savefig(figname+'.pdf')
    fig.savefig(figname+'.png')
    fig.close()

def get_source(file_list, target, selavy_table, folder=None, scaling_factor=1.0):
    
    fname = ''
    if folder:
        fname = '{}{}_sl.fits'.format(folder, target['comp_name'])
        if not fname in file_list:
            return None
    
    comp_cat_rows = selavy_table[selavy_table['component_name']==target['comp_name']]
    comp_cat_row = None
    for row in comp_cat_rows:
        if comp_cat_row is None or row['flux_peak'] > comp_cat_row['flux_peak']:
            comp_cat_row = row
    src = {'ra':comp_cat_row['ra_deg_cont'], 'dec':comp_cat_row['dec_deg_cont'], 
           'a':comp_cat_row['maj_axis']/2*scaling_factor, 'b':comp_cat_row['min_axis']/2*scaling_factor, 'pa': comp_cat_row['pos_ang'],
          'comp_name': target['comp_name'], 'fname': fname, 'flux_peak': comp_cat_row['flux_peak']}
    return src

def save_spectrum(velocity, opacity, flux, em_mean, em_std, filename, 
                   sigma_tau):
    """
    Write the spectrum (velocity, flux and opacity) to a votable format file.

    :param spectrum: The spectrum to be output.
    :param opacity:  The opacity to be output.
    :param filename:  The filename to be created
    :param longitude: The galactic longitude of the target object
    :param latitude: The galactic latitude of the target object
    """
    table = Table(meta={'name': filename, 'id': 'opacity'})
    table.add_column(Column(name='velocity', data=velocity, unit='m/s', description='LSRK velocity'))
    table.add_column(Column(name='opacity', data=opacity, description='Absorption as exp(-tau)'))
    table.add_column(Column(name='flux', data=flux, unit='Jy', description='Flux per beam'))
    table.add_column(Column(name='sigma_opacity', data=sigma_tau, description='Noise in the absorption profile, relative to exp(-tau)'))
    table.add_column(Column(name='em_mean', data=em_mean, unit='K', description='Mean brightness temperature around the source'))
    table.add_column(Column(name='em_std', data=em_std, unit='K', description='Noise level in the brightness temperature'))

    votable = from_table(table)
    writeto(votable, filename)


def extract_spectrum(fname, src, continuum_start_vel, continuum_end_vel, figures_folder, num_edge_chan=10):

    hdulist = fits.open(fname)
    image = hdulist[0].data
    header = hdulist[0].header
    w = WCS(header)
    index = np.arange(header['NAXIS3'])
    velocities = w.wcs_pix2world(10,10,index[:],0,0)[2]
    print (image.shape)
    print (continuum_start_vel, continuum_end_vel, velocities[0], velocities[-1])
    radius=math.ceil(min(max(6, src['a']), (image.shape[-1]-1)//2))

    img_slice = cube_tools.get_integrated_spectrum(image, w, src, velocities, continuum_start_vel, continuum_end_vel, 
        radius=radius, plot_weight_path=figures_folder)

    l_edge, r_edge = cube_tools.find_edges(img_slice, num_edge_chan)

    spectrum_array = rec.fromarrays(
        [np.arange(img_slice.size)[l_edge:r_edge],
         velocities[l_edge:r_edge],
         img_slice[l_edge:r_edge]],
        names='plane,velocity,flux')

    del image
    del header
    hdulist.close()

    return spectrum_array


def highlight_features(ax, ranges):
    if ranges is not None:
        ylim = ax.get_ylim()
        colors = list(mcolors.TABLEAU_COLORS.keys())
        #y = np.arange(0.0, 2, 0.01)
        for idx, rng in enumerate(ranges):
            color_idx = (idx+1) % len(colors)
            ax.fill_betweenx(ylim, rng[0], rng[1], color=colors[color_idx], alpha=0.2, zorder=-1)


def plot_combined_spectrum(velocity, em_mean, em_std, opacity, sigma_opacity, filename, title, vel_range=None, ranges=None):
    fig, axs = plt.subplots(2, 1, figsize=(8, 10))
    
    axs[0].set_title(title, fontsize=16)

    axs[0].axhline(1, color='r')
    axs[0].plot(velocity, opacity, zorder=4, lw=1)
    axs[0].fill_between(velocity, 1-sigma_opacity, 1+sigma_opacity, color='grey', alpha=0.4, zorder=1)
    axs[0].plot(velocity, 1-(3*sigma_opacity), zorder=1, lw=1, ls=":")
    axs[0].set_ylabel(r'$e^{(-\tau)}$')
    axs[0].grid(True)
    highlight_features(axs[0], ranges)

    axs[1].plot(velocity, em_mean, color='firebrick', zorder=4)
    axs[1].fill_between(velocity, em_mean-em_std, em_mean+em_std, color='grey', alpha=0.4, zorder=1)
    axs[1].set_ylabel(r'$T_B (K)$')
    axs[1].set_xlabel(r'Velocity relative to LSR (km/s)')
    axs[1].grid(True)
    highlight_features(axs[1], ranges)
    
    if vel_range is None:
        vel_range = (np.min(velocity), np.max(velocity))
    axs[0].set_xlim(vel_range)
    axs[1].set_xlim(vel_range)

    plt.savefig(filename, bbox_inches='tight')
    #plt.show()
    plt.close()
    return

def output_spectrum(spec_folder, spectrum_array, opacity, sigma_opacity, em_mean, em_std, comp_name, src_id, continuum_start_vel, continuum_end_vel):
    filename = '{}/{}_spec'.format(spec_folder, comp_name)

    save_spectrum(spectrum_array.velocity, opacity, spectrum_array.flux, em_mean, em_std, filename+'.vot', sigma_opacity)
    spectrum_tools.plot_absorption_spectrum(spectrum_array.velocity, opacity, filename+'.png', comp_name, continuum_start_vel, continuum_end_vel, sigma_opacity, range=None)
    spectrum_tools.plot_absorption_spectrum(spectrum_array.velocity, opacity, filename+'_zoom.png', comp_name, continuum_start_vel, continuum_end_vel, sigma_opacity, range=(75,275))


def match_emission_to_absorption(em_mean, em_std, em_velocity, abs_velocity):
    em_vel_pos = (em_velocity[1] - em_velocity[0]) > 0
    abs_vel_pos = (abs_velocity[1] - abs_velocity[0]) > 0
    if abs_vel_pos != em_vel_pos:
        #print("Reversing emission data")
        em_mean = np.flip(em_mean, 0)
        em_std = np.flip(em_std, 0)
        em_velocity = np.flip(em_velocity, 0)
    # tODO: Match the emission velocities against the absoprtion ones. Pad or trim the emisiosn as needed to get a consistent array size with the absorption velocities
    # Can we assume these are already in the same scale? Maybe not if we intend to do higher spectral resolution spectra for some sources
    em_mean_interp = np.interp(abs_velocity, em_velocity, em_mean, left=0, right=0)
    em_std_interp = np.interp(abs_velocity, em_velocity, em_std, left=0, right=0)
    return em_mean_interp, em_std_interp


def calc_sigma_tau(cont_sd, em_mean, opacity):
    """
    Calculate the noise in the absorption profile at each velocity step. Where emission data is available, this is
    based on the increased antenna temperature due to the received emission.

    :param cont_sd: The measured noise in the continuum region of the spectrum in absorption units.
    :param em_mean: The mean emission brightness temperature in K
    :param opacity: The optical depth spectrum, used only for the shape of the data
    :return: A numpy array containing the noise level in the optical depth data at each velocity step.
    """

    # ATCA figures tsys = 44.7 eta=0.5
    tsys = 50.0 # Preliminary PAF figures from ASKAP Pilot Survey Plan, Hotan et al 2018
    antenna_eff = 0.7 # 0.5
    if len(em_mean) > 0:
        floor = np.zeros(em_mean.shape)
        sigma_tau = cont_sd * ((tsys + np.fmax(floor, antenna_eff*em_mean)) / tsys)
    else:
        sigma_tau = np.full(opacity.shape, cont_sd)
    return sigma_tau


def check_noise(opacity, sigma_opacity):
    """
    Check if the noise estimate of the spectrum is poor. A poor estimate will have a higher standard deviation than the
    noise estimate even when any potential absoprtion is clipped out. Only low (less than 10) signal to noise spectrum will be flagged as 
    having poor noise.

    :param opacity: The optical depth spectrum
    :param sigma_opacity: The noise lecel at each velocity step
    :return true if the noise estimate is poor, false if it is good
    """

    opacity_ar = np.asarray(opacity)
    noise = np.array(sigma_opacity)

    signal_to_noise = (1-opacity_ar) / noise
    max_sn = np.max(signal_to_noise)

    clipped = np.array(opacity_ar, copy=True)
    clipped[clipped<1-sigma_opacity] = 1
    clipped_mean = np.mean(clipped)
    clipped_std = np.std(clipped)
    clipped_ratio = clipped_std/np.min(noise)
    return (clipped_ratio > 1) & (max_sn < 10)


def extract_all_spectra(targets, file_list, cutouts_folder, selavy_table, figures_folder,spectra_folder, sbid, max_spectra = 5000, cont_range=(-100,-60)):
    print('Processing {} cutouts into spectra, input:{}'.format(len(targets), cutouts_folder))

    i = 0
    for tgt in targets:
        comp_name = tgt['comp_name']
        src = get_source(file_list, tgt, selavy_table, folder=cutouts_folder, scaling_factor=0.8)
        if not src:
            print('Skipping missing src #{} {}'.format(tgt['id'], comp_name))
            continue
        i+=1
        if i > max_spectra:
            print ("Reaching maximum spectra for this run")
            break
            
        print('\nExtracting spectrum for src #{} {} maj:{} min:{} pa:{:.3f}'.format(tgt['id'], comp_name, src['a'], src['b'], src['pa']))

        plot_mom0(src['fname'], comp_name, figures_folder, src['ra'], src['dec'], 
                src['a'], src['b'], src['pa'])
        
        # Extract the spectrum
        #continuum_start_vel = -100*u.km.to(u.m)
        #continuum_end_vel = -60*u.km.to(u.m)
        continuum_start_vel = cont_range[0]*u.km.to(u.m)
        continuum_end_vel = cont_range[1]*u.km.to(u.m)
        spectrum = extract_spectrum(src['fname'], src, continuum_start_vel, continuum_end_vel, figures_folder)
        
        mean_cont, sd_cont = spectrum_tools.get_mean_continuum(spectrum.velocity, spectrum.flux, continuum_start_vel, continuum_end_vel)
        opacity = spectrum.flux/mean_cont
        sigma_opacity = sd_cont*np.ones(spectrum.velocity.shape)

        # Read the primary beam emission spectrum
        filename = '{0}/{1}_pb_emission.vot'.format(spectra_folder, comp_name)
        emission_votable = votable.parse(filename, pedantic=False)
        emission_tab = emission_votable.get_first_table().to_table()
        pb_em_mean, pb_em_std = match_emission_to_absorption(emission_tab['pb_em_mean'], emission_tab['pb_em_std'], emission_tab['velocity'], spectrum.velocity)

        # Calculate the noise envelope
        sigma_opacity = calc_sigma_tau(sd_cont, pb_em_mean, opacity)

        # Read the emission spectrum
        filename = '{0}/{1}_emission.vot'.format(spectra_folder, comp_name)
        if os.path.exists(filename):
            emission_votable = votable.parse(filename, pedantic=False)
            emission_tab = emission_votable.get_first_table().to_table()
            em_mean, em_std = match_emission_to_absorption(emission_tab['em_mean'], emission_tab['em_std'], emission_tab['velocity'], spectrum.velocity)
        else:
            em_mean = em_std = np.zeros(spectrum.velocity.shape)

        output_spectrum(spectra_folder, spectrum, opacity, sigma_opacity, em_mean, em_std, comp_name, tgt['id'], continuum_start_vel, continuum_end_vel)

    return None


def find_runs(spec_table, min_sigma=3, min_len=2):
    runs = []
    curr_idx = 0
    single_min_sigma = min_sigma[0] if isinstance(min_sigma, list) else min_sigma
    run_min_sigma = min_sigma[1] if isinstance(min_sigma, list) else min_sigma

    single_sigma_filter = 1-spec_table['opacity'] > single_min_sigma*spec_table['sigma_opacity']
    run_sigma_filter = 1-spec_table['opacity'] > run_min_sigma*spec_table['sigma_opacity']
    sigma = (1-spec_table['opacity']) / spec_table['sigma_opacity']
    #sigma3 = 1-spec_table['opacity'] > min_sigma*spec_table['sigma_opacity']
    for bit, group in itertools.groupby(run_sigma_filter):
        result = list(group)
        length = sum(result)
        #print (bit, group, length, result)
        if length >= min_len:
            # Check at least one channel passes the single channel threshold
            if np.sum(single_sigma_filter[curr_idx:curr_idx+len(result)]) > 0:
                vel_start = spec_table[curr_idx]['velocity']/1000
                vel_end = spec_table[curr_idx+length-1]['velocity']/1000
                max_sigma = np.max(sigma[curr_idx:curr_idx+length])
                runs.append(AbsRun(length, curr_idx, vel_start, vel_end, max_sigma))
                #print ("Found run of {} at {} from {:.2f} to {:.2f} km/s".format(length, curr_idx, vel_start, vel_end))
        curr_idx += len(result)
    return runs


def merge_runs(runs_a, runs_b):
    fullset = list(runs_a)
    for abs_run in runs_b:
        #print (abs_run)
        present = False
        for existing in runs_a:
            #print ("..", existing)
            if existing.start_vel < abs_run.end_vel and existing.end_vel > abs_run.start_vel:
                # features overlap
                present = True
        if not present:
            fullset.append(abs_run)
    return fullset

def assess_spectra(targets, file_list, selavy_table, figures_folder, spectra_folder, sbid, cont_range=(-100,-60)):
    src_table = QTable(names=('id', 'comp_name', 'ra', 'dec', 'rating', 'flux_peak', 'mean_cont', 'sd_cont', 'opacity_range', 
            'max_s_max_n', 'max_noise', 'num_chan_noise', 'min_opacity', 'vel_min_opacity', 'has_mw_abs', 'has_other_abs', 
            'semi_maj_axis', 'semi_min_axis', 'pa', 'n_h', 'noise_flag'),
            dtype=('int', 'U32', 'float64', 'float64', 'str', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 
            'float64', 'float64', 'float64', 'bool', 'bool', 'float64', 'float64', 'float64', 'float64', None),
            meta={'name': 'ASKAP Spectra for sbid {}'.format(sbid)})

    abs_table = QTable(names=('id', 'comp_name', 'abs_name', 'ra', 'dec', 'rating', 'flux_peak', 'mean_cont', 'sd_cont', 'opacity_range', 
            'max_s_max_n', 'max_noise', 'num_chan_noise', 'semi_maj_axis', 'semi_min_axis', 'pa', 
            'start_vel', 'end_vel', 'length', 'min_opacity', 'max_sigma'),
            dtype=('int', 'U32', 'U32', 'float64', 'float64', 'str', 'float64', 'float64', 'float64', 'float64', 
            'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 
            'float64', 'float64', 'int', 'float64', 'float64'),
            meta={'name': 'ASKAP Absorption detections for sbid {}'.format(sbid)})

    print('Assessing {} spectra'.format(len(targets)))

    i = 0
    for tgt in targets:
        comp_name = tgt['comp_name']
        src = get_source(file_list, tgt, selavy_table, scaling_factor=0.8)
        if not src:
            print('Skipping missing src #{} {}'.format(tgt['id'], tgt['comp_name']))
            continue

        abs_spec_filename = '{}/{}_spec.vot'.format(spectra_folder, comp_name)
        if not os.path.exists(abs_spec_filename):
            print ('No absorption spectrum found for {} at {}'.format(comp_name, abs_spec_filename))
            continue
        abs_spec_votable = parse_single_table(abs_spec_filename)
        abs_spec = abs_spec_votable.to_table()

        opacity = abs_spec['opacity']
        sigma_opacity = abs_spec['sigma_opacity'] 
        min_opacity = np.min(opacity)
        vel_min_opacity = abs_spec['velocity'][np.argmin(opacity)]
        max_opacity = np.max(opacity)
        num_noise = (opacity > 1+sigma_opacity).sum()

        #continuum_start_vel = -100*u.km.to(u.m)
        #continuum_end_vel = -60*u.km.to(u.m)
        continuum_start_vel = cont_range[0]*u.km.to(u.m)
        continuum_end_vel = cont_range[1]*u.km.to(u.m)
        mean_cont, sd_cont = spectrum_tools.get_mean_continuum(abs_spec['velocity'], abs_spec['flux'], continuum_start_vel, continuum_end_vel)

        poor_noise_flag = check_noise(opacity, sigma_opacity)

        rating, opacity_range, max_s_max_n = spectrum_tools.rate_spectrum(opacity, np.mean(sigma_opacity))

        has_mw_abs = False
        has_other_abs = False
        runs = find_runs(abs_spec, min_sigma=[3, 2.8], min_len=2)
        if (len(runs) > 0):
            for idx, absrow in enumerate(runs):
                abs_name = '{}_{}'.format(tgt['comp_name'], int(absrow.start_vel))
                if absrow.start_vel < 50:
                    has_mw_abs = True
                if absrow.start_vel >= 50:
                    has_other_abs = True
                abs_table.add_row([tgt['id'], tgt['comp_name'], abs_name, src['ra']*u.deg, src['dec']*u.deg, 
                    rating, src['flux_peak']*u.Jy, mean_cont*u.Jy, sd_cont, opacity_range, max_s_max_n, max_opacity, num_noise, 
                    src['a'], src['b'], src['pa'], absrow.start_vel*u.km/u.s, absrow.end_vel*u.km/u.s, absrow.length, 
                    np.nanmin(opacity[absrow.start_idx:absrow.start_idx+absrow.length]), absrow.max_sigma])

        
        src_table.add_row([tgt['id'], tgt['comp_name'], src['ra']*u.deg, src['dec']*u.deg, 
                        rating, src['flux_peak']*u.Jy, mean_cont*u.Jy, sd_cont, opacity_range, max_s_max_n, max_opacity, num_noise, 
                        min_opacity, vel_min_opacity, has_mw_abs, has_other_abs, src['a'], src['b'], src['pa'], 0, poor_noise_flag])

    for rating in 'ABCDEF':
        count = np.sum(src_table['rating'] == rating)
        print ("  Rating {}: {:.0f} sources".format(rating, count))

    hits = np.sum(src_table['has_mw_abs'])
    print ("  Has MW absorption: {:.0f} sources or {:.1f}%".format(hits, 100*hits/len(targets)))
    hits = np.sum(src_table['has_other_abs'])
    print ("  Has other absorption: {:.0f} sources or {:.1f}%".format(hits, 100*hits/len(targets)))
    print ("  Found {} absorption instances".format(len(abs_table)))

    return src_table, abs_table


def add_column_density(spectra_table):
    # read in mom0 map
    gass_mom0 = fits.open('hi_zea_ms_mom0.fits')
    gass_wcs = WCS(gass_mom0[0].header)
    gass_no_nan_data = np.nan_to_num(gass_mom0[0].data)
    gass_nh_data = gass_no_nan_data * 1.82 * 10**18 /1e3

    # Calculate all column density values
    pix_pos = gass_wcs.wcs_world2pix(spectra_table['ra'], spectra_table['dec'], 0)
    if pix_pos and np.min(pix_pos) >= 0 and np.max(pix_pos[0]) < gass_nh_data.shape[1] and np.max(pix_pos[1]) < gass_nh_data.shape[0]:
        n_h_vals = gass_nh_data[pix_pos[1].astype(int), pix_pos[0].astype(int)]
        spectra_table['n_h'] = n_h_vals
    else:
        print ("No column density data available for this region")

def add_col_metadata(vo_table, col_name, description, units=None, ucd=None, datatype=None):
    col = vo_table.get_first_table().get_field_by_id(col_name)
    col.description = description
    if units:
        col.unit = units
    if ucd:
        col.ucd = ucd
    if datatype:
        col.datatype = datatype


def add_spectra_column_metadata(spectra_vo_table):
    add_col_metadata(spectra_vo_table, 'id', 'Unique identifier of the source.')
    add_col_metadata(spectra_vo_table, 'comp_name', 'Background source component name.', ucd='meta.id;meta.main')
    add_col_metadata(spectra_vo_table, 'ra', 'J2000 right ascension in decimal degrees.', units='deg', ucd='pos.eq.ra;meta.main')
    add_col_metadata(spectra_vo_table, 'dec', 'J2000 declination in decimal degrees.', units='deg', ucd='pos.eq.dec;meta.main')
    add_col_metadata(spectra_vo_table, 'rating', 'Quality rating of the absorption spectrum.')
    add_col_metadata(spectra_vo_table, 'flux_peak', 'Peak flux density of the background source.', units='mJy/beam', ucd='phot.flux.density;stat.max;em.radio;stat.fit')
    add_col_metadata(spectra_vo_table, 'mean_cont', 'Mean continuum level per channel of the spectrum.', units='mJy/beam', ucd='phot.flux.density;stat.mean;em.radio;stat.fit')
    add_col_metadata(spectra_vo_table, 'sd_cont', '1-sigma noise level of the spectrum in opacity units. Does not include emission noise.')
    add_col_metadata(spectra_vo_table, 'opacity_range', 'The range in the absorption spectrum from highest emission to lowest absorption.')
    add_col_metadata(spectra_vo_table, 'max_s_max_n', 'Ratio of maximum absorption versus maximum emission noise in the spectrum.')
    add_col_metadata(spectra_vo_table, 'max_noise', 'The maximum emission noise in the spectrum, in opacity units.')
    add_col_metadata(spectra_vo_table, 'num_chan_noise', 'The number of channels in the spectrum that have emission above the 1-sigma noise level.')
    add_col_metadata(spectra_vo_table, 'min_opacity', 'The minimum opacity in the spectrum, in opacity units, representing the peak absorption.')
    add_col_metadata(spectra_vo_table, 'vel_min_opacity', 'The velocity at which the peak absorption is found.', units='km/s')
    add_col_metadata(spectra_vo_table, 'has_mw_abs', 'Flag to indicate that the spectrum has absorption in the Milky Way velocity range (<50 km/s).', datatype='boolean')
    add_col_metadata(spectra_vo_table, 'has_other_abs', 'Flag to indicate that the spectrum has absorption away from the Milky Way velocity range (>50 km/s).', datatype='boolean')
    add_col_metadata(spectra_vo_table, 'semi_maj_axis', 'HWHM major axis before deconvolution.', units='arcsec', ucd='phys.angSize.smajAxis;em.radio;stat.fit')
    add_col_metadata(spectra_vo_table, 'semi_min_axis', 'HWHM minor axis before deconvolution.', units='arcsec', ucd='phys.angSize.sminAxis;em.radio;stat.fit')
    add_col_metadata(spectra_vo_table, 'pa', 'Position angle before deconvolution measured as degrees CCW of North.', units='deg', ucd='phys.angSize;pos.posAng;em.radio;stat.fit')
    add_col_metadata(spectra_vo_table, 'n_h', 'Column density towards the background source excluding Milky Way velocities, from the GASS survey.', units='cm-2', ucd='phys.columnDensity')
    add_col_metadata(spectra_vo_table, 'noise_flag', 'Flag to indicate that the noise for the spectrum may be underestimated.', datatype='boolean')


def add_absorption_column_metadata(abs_vo_table):
    add_col_metadata(abs_vo_table, 'id', 'Unique identifier of the absorption feature.')
    add_col_metadata(abs_vo_table, 'comp_name', 'Background source component name.', ucd='meta.id')
    add_col_metadata(abs_vo_table, 'abs_name', 'Absorption feature name.', ucd='meta.id;meta.main')
    add_col_metadata(abs_vo_table, 'ra', 'J2000 right ascension in decimal degrees.', units='deg', ucd='pos.eq.ra;meta.main')
    add_col_metadata(abs_vo_table, 'dec', 'J2000 declination in decimal degrees.', units='deg', ucd='pos.eq.dec;meta.main')
    add_col_metadata(abs_vo_table, 'rating', 'Quality rating of the absorption spectrum.')
    add_col_metadata(abs_vo_table, 'flux_peak', 'Peak flux density of the background source.', units='mJy/beam', ucd='phot.flux.density;stat.max;em.radio;stat.fit')
    add_col_metadata(abs_vo_table, 'mean_cont', 'Mean continuum level per channel of the spectrum.', units='mJy/beam', ucd='phot.flux.density;stat.mean;em.radio;stat.fit')
    add_col_metadata(abs_vo_table, 'sd_cont', '1-sigma noise level of the spectrum in opacity units. Does not include emission noise.')
    add_col_metadata(abs_vo_table, 'opacity_range', 'The range in the absorption spectrum from highest emission to lowest absorption.')
    add_col_metadata(abs_vo_table, 'max_s_max_n', 'Ratio of maximum absorption versus maximum emission noise in the spectrum.')
    add_col_metadata(abs_vo_table, 'max_noise', 'The maximum emission noise in the spectrum, in opacity units.')
    add_col_metadata(abs_vo_table, 'num_chan_noise', 'The number of channels in the spectrum that have emission above the 1-sigma noise level.')
    add_col_metadata(abs_vo_table, 'semi_maj_axis', 'HWHM major axis before deconvolution.', units='arcsec', ucd='phys.angSize.smajAxis;em.radio;stat.fit')
    add_col_metadata(abs_vo_table, 'semi_min_axis', 'HWHM minor axis before deconvolution.', units='arcsec', ucd='phys.angSize.sminAxis;em.radio;stat.fit')
    add_col_metadata(abs_vo_table, 'pa', 'Position angle before deconvolution measured as degrees CCW of North.', units='deg', ucd='phys.angSize;pos.posAng;em.radio;stat.fit')
    add_col_metadata(abs_vo_table, 'start_vel', 'The lowest velocity covered by this absorption feature.', units='km/s', ucd='phys.velocity;stat.min')
    add_col_metadata(abs_vo_table, 'end_vel', 'The highest velocity covered by this absorption feature.', units='km/s', ucd='phys.velocity;stat.max')
    add_col_metadata(abs_vo_table, 'length', 'The number of channels included in this absorption feature.', units='chan')
    add_col_metadata(abs_vo_table, 'min_opacity', 'The minimum opacity in the feature, in opacity units, representing the peak absorption.')
    add_col_metadata(abs_vo_table, 'max_sigma', 'The maximum significance of the absorption feature, representing the peak absorption.', ucd='stat.snr')


def prep_folders(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.mkdir(folder)
            print ("Created " + folder)


def recenter(ax, wcs, x, y, radius=None, width=None, height=None):
    '''
    Center the image on a given position and with a given radius.

    Either the radius or width/heigh arguments should be specified.

    Parameters
    ----------

    x, y : float
        Coordinates to center on

    radius : float, optional
        Radius of the region to view. This produces a square plot.

    width : float, optional
        Width of the region to view. This should be given in
        conjunction with the height argument.

    height : float, optional
        Height of the region to view. This should be given in
        conjunction with the width argument.
    '''

    xpix, ypix = wcs.wcs_world2pix(x, y, 0)
    
    pix_scale = proj_plane_pixel_scales(wcs)
    sx, sy = pix_scale[1], pix_scale[0]

    if radius:
        dx_pix = radius / sx
        dy_pix = radius / sy
    elif width and height:
        dx_pix = width / sx * 0.5
        dy_pix = height / sy * 0.5
    else:
        raise Exception("Need to specify either radius= or width= and height= arguments")

    if (xpix + dx_pix < -0.5 or
        xpix - dx_pix > wcs.array_shape[1] - 0.5 or
        ypix + dy_pix < -0.5 or
            ypix - dy_pix > wcs.array_shape[1]):

        raise Exception("Zoom region falls outside the image")

    ax.set_xlim(xpix - dx_pix, xpix + dx_pix)
    ax.set_ylim(ypix - dy_pix, ypix + dy_pix)


def add_source_box(ax, sources):
    min_ra = np.min(sources['ra'])
    max_ra = np.max(sources['ra'])
    min_dec = np.min(sources['dec'])
    max_dec = np.max(sources['dec'])
    print (min_ra, max_ra, min_dec, max_dec)
    width = max_ra-min_ra
    height = max_dec-min_dec
    centre_ra = min_ra + width/2
    centre_dec = min_dec + height/2
    print (centre_ra, centre_dec, width, height)

    r = Rectangle((min_ra, min_dec), width, height, edgecolor='green', facecolor='none',
              transform=ax.get_transform('fk5'), lw=3)
    ax.add_patch(r)


def get_field_centre(sources):
    min_ra = np.min(sources['ra'])
    max_ra = np.max(sources['ra'])
    min_dec = np.min(sources['dec'])
    max_dec = np.max(sources['dec'])
    print (min_ra, max_ra, min_dec, max_dec)
    width = max_ra-min_ra
    height = max_dec-min_dec
    centre_ra = min_ra + width/2
    centre_dec = min_dec + height/2
    field_centre = SkyCoord(ra=centre_ra*u.deg, dec=centre_dec*u.deg, frame=FK4)
    print (field_centre)
    return field_centre


def add_sources(ax, sources,detection, best, is_mw):

    for src in sources:
        colour_name = 'darkgreen'
        facecolor = 'darkgreen'
        sigma = (1-src['min_opacity'])/src['sd_cont']
    
        marker = '.'
        if src['has_other_abs'] and not (is_mw and src['has_mw_abs']):
            marker = 'o'
            facecolor = 'none'
            colour_name = 'orange'

        elif src['has_mw_abs']:
            marker = 's'
            facecolor = 'none'
            colour_name = 'yellow'
            
        elif src['rating'] <= 'B':
            marker = '+'
            facecolor = 'darkgreen'

        ax.scatter([src['ra']], [src['dec']], transform=ax.get_transform('world'), 
               marker=marker, edgecolor=colour_name, facecolor=facecolor)


def plot_background_map(fig, background):

    # moment zero map
    mom0 = fits.open(background)

    no_nan_data = np.nan_to_num(mom0[0].data)
    nh_data = no_nan_data * 1.82 * 10**18 /1e3
    vmin=np.percentile(nh_data, 0.25)
    vmax=np.percentile(nh_data, 99.75)

    # Create an ImageNormalize object
    vmid = vmin - (vmax - vmin) / 30.
    asinh_a = (vmid - vmin) / (vmax - vmin)
    norm_kwargs = {'asinh_a': abs(asinh_a)}
    norm = simple_norm(nh_data, 'asinh', min_cut=vmin, max_cut=vmax, clip=False,
                                    **norm_kwargs)

    wcs = WCS(mom0[0].header)
    ax = fig.add_subplot(111, projection=wcs)
    im = ax.imshow(nh_data, cmap='gist_yarg', vmin=vmin, vmax=vmax, alpha=0.6, norm=norm)

    lon = ax.coords[0]
    lat = ax.coords[1]
    #lon.set_ticks(number=12)
    lon.set_ticks_position('rb')
    lon.set_ticklabel_position('rb')
    lon.set_axislabel_position('b')

    #lat.set_ticks(number=12)
    lat.set_ticks_position('tl')
    lat.set_ticklabel_position('tl')
    lat.set_axislabel_position('l')

    # Add axes labels
    ax.set_xlabel("Right Ascension (hours)", fontsize=16)
    ax.set_ylabel("Declination (degrees)", fontsize=16)

    ax.grid()
    return ax, wcs


def plot_source_loc_map(spectra_table, figures_folder, is_mw, background='hi_zea_ms_mom0.fits'):

    print('\nPlotting {} source locations over background of {}.'.format(len(spectra_table), background))

    fig = plt.figure(figsize=(10.5, 9))
    ax, wcs = plot_background_map(fig, background)

    field_centre = get_field_centre(spectra_table)
    try:
        recenter(ax, wcs, field_centre.ra.value, field_centre.dec.value, width=8.25, height=7.25)  # degrees
    except Exception as ex:
        print(ex)
        return

    # Display the moment map image
    #plt.colorbar(im,fraction=0.046, pad=0.04)
    add_sources(ax, spectra_table, 3, 5, is_mw)

    prefix = figures_folder + "/source_loc"
    print ("  Output", prefix+".png")
    fig.savefig(prefix + ".png", bbox_inches='tight')
    fig.savefig(prefix + ".pdf", bbox_inches='tight')


def plot_field_loc_map(spectra_table, figures_folder, background='hi_zea_ms_mom0.fits'):

    print('\nPlotting field locations over background of {}.'.format(background))

    fig = plt.figure(figsize=(10.5, 9))
    ax, wcs = plot_background_map(fig, background)

    # Display the moment map image
    #plt.colorbar(im,fraction=0.046, pad=0.04)
    add_source_box(ax, spectra_table)

    prefix = figures_folder + "/field_loc"
    print ("  Output", prefix+".png")
    fig.savefig(prefix + ".png", bbox_inches='tight')
    fig.savefig(prefix + ".pdf", bbox_inches='tight')


def output_emission_spectra(spectra_table, spectra_folder):
    for idx, source in enumerate(spectra_table):
        comp_name = source['comp_name']
        filename = '{0}/{1}_emission.vot'.format(spectra_folder, comp_name)
        emission_votable = votable.parse(filename, pedantic=False)
        emission_tab = emission_votable.get_first_table().to_table()

        velocities = emission_tab['velocity']
        tb_mean = emission_tab['em_mean']
        tb_std = emission_tab['em_std']

        if np.sum(tb_mean) == 0:
            print ("Source {} has all no emisision data".format(comp_name))

        abs_spec_filename = '{}/{}_spec.vot'.format(spectra_folder, comp_name)
        abs_spec_votable = parse_single_table(abs_spec_filename)
        abs_spec = abs_spec_votable.to_table()

        title = 'Source #{} {}'.format(source['id'], comp_name)
        filename = '{}/{}_combined.png'.format(spectra_folder, comp_name)
        plot_combined_spectrum(velocities/1000, tb_mean, tb_std, 
            abs_spec['velocity']/1000, abs_spec['opacity'], abs_spec['sigma_opacity'], 
            filename, title)


def plot_all_spectra(spectra_table, abs_table, spectra_folder, no_zoom):
    for idx, source in enumerate(spectra_table):
        comp_name = source['comp_name']

        abs_spec_filename = '{}/{}_spec.vot'.format(spectra_folder, comp_name)
        abs_spec_votable = parse_single_table(abs_spec_filename)
        abs_spec = abs_spec_votable.to_table()
        abs_velocity = abs_spec['velocity']

        tgt_abs = abs_table[abs_table['comp_name']==comp_name]
        start_vel = tgt_abs['start_vel'].data
        end_vel = tgt_abs['end_vel'].data
        ranges = np.stack((start_vel,end_vel), axis=-1)
        vel_ranges = None if no_zoom else (75,350)

        title = 'Source #{} {}'.format(source['id'], comp_name)
        filename = '{}/{}_combined.png'.format(spectra_folder, comp_name)
        plot_combined_spectrum(abs_velocity/1000, abs_spec['em_mean'], abs_spec['em_std'], 
            abs_spec['opacity'], abs_spec['sigma_opacity'], 
            filename, title, ranges=ranges, vel_range=vel_ranges)


def output_reg_file(filename, sources):
    with open(filename, 'w') as writer:
        for src in sources:
            smaj = src['semi_maj_axis']*u.arcsec.to(u.deg)
            smin = src['semi_min_axis']*u.arcsec.to(u.deg)
            shape = 'j2000;ellipse {} {} {} {} {} # text="{}"\n'.format(
                src['ra'], src['dec'], smaj, smin, src['pa'], src['comp_name'])
            writer.write(shape)

def export_ds9_regions(spectra_table, parent_folder, detect_field='has_mw_abs'):

    print('\nOutputting DS9 region files.')

    detections = spectra_table[spectra_table[detect_field]]
    filename = parent_folder+'detections.reg'
    print (" Outputting {} detections to {}".format(len(detections), filename))
    output_reg_file(filename, detections)

    non_detections = spectra_table[~spectra_table[detect_field]]
    filename = parent_folder+'non_detections.reg'
    print (" Outputting {} non detections to {}".format(len(non_detections), filename))
    output_reg_file(filename, non_detections)


def check_milky_way(selavy_table):
    loc = SkyCoord(ra=selavy_table['ra_deg_cont'], dec=selavy_table['dec_deg_cont'], frame='icrs')
    is_mw = np.min(loc.galactic.b.value) > -20
    return is_mw

def main():
    warnings.simplefilter('ignore', category=AstropyWarning)

    args = parseargs()

    start = time.time()
    print("#### Started ASKAP spectra extraction for sbid {} at {} ####".format(args.sbid,
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))))

    parent_folder = 'sb{}/'.format(args.sbid)
    if args.parent:
        parent_folder = args.parent + '/'
    if not os.path.exists(parent_folder):
        print("Error: Folder {} does not exist.".format(parent_folder))
        return 1

    cutout_folder = parent_folder + 'cutouts/'
    figures_folder = parent_folder + 'figures/'
    spectra_folder = parent_folder + 'spectra/'
    prep_folders([spectra_folder, figures_folder])

    file_list = glob.glob(cutout_folder + 'J*_sl.fits')
    file_list.sort()

    targets = read_targets('{}/targets_{}.csv'.format(parent_folder, args.sbid))
    
    # Read and filter catalogue
    src_votable = votable.parse(args.catalogue, pedantic=False)
    selavy_table = src_votable.get_first_table().to_table()
    rename_columns(selavy_table)
    continuum_range = (args.continuum_start, args.continuum_end)
    is_milky_way = check_milky_way(selavy_table)
    print ('is milky way?', is_milky_way)

    if not args.skip_abs:    
        # Extract absorption spectra (including noise and emission)
        spectra = extract_all_spectra(targets, file_list, cutout_folder, selavy_table, figures_folder, spectra_folder, 
            args.sbid, cont_range=continuum_range)
    else:
        print ("**Skipping spectra extraction - reusing existing spectra")

    # Assess spectra - rating, consecutive significant channels (flag if mw or other abs)
    spectra_table, abs_table = assess_spectra(targets, file_list, selavy_table, figures_folder, spectra_folder, 
        args.sbid,  cont_range=continuum_range)
    add_column_density(spectra_table)

    # Save spectra catalogue
    spectra_vo_table = from_table(spectra_table)
    add_spectra_column_metadata(spectra_vo_table)
    writeto(spectra_vo_table, parent_folder+'askap_spectra.vot')

    # Save absorption catalogue
    abs_vo_table = from_table(abs_table)
    add_absorption_column_metadata(abs_vo_table)
    writeto(abs_vo_table, parent_folder+'askap_absorption.vot')

    # Produce consolidated plots
    mom0_file = 'hi_zea_all_mom0.fits' if is_milky_way else 'hi_zea_ms_mom0.fits'
    plot_source_loc_map(spectra_table, figures_folder, is_milky_way, background=mom0_file)
    plot_field_loc_map(spectra_table, figures_folder, background=mom0_file)
    plot_all_spectra(spectra_table, abs_table, spectra_folder, args.no_zoom)

    # DS9 files
    export_ds9_regions(spectra_table, parent_folder)

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Extracted %d spectra and %d absorption instances in %.02f s' %
          (len(spectra_table), len(abs_table), end - start))

    return 0



if __name__ == '__main__':
    exit(main())
    
