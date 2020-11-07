# Python script to produce spectra for a series of sub-cubes.
# The spectra are assumed to be absoprtion spectra and rated according to their suitability for absorption detection.

# Author James Dempsey
# Date 19 Feb 2020

import argparse
import csv 
import glob
import math
import os
import time
import warnings

import astropy.units as u
from astropy.coordinates import SkyCoord, FK4
from astropy.io import fits, votable
from astropy.io.ascii import Csv
from astropy.io.votable.tree import Param,Info
from astropy.io.votable import from_table, writeto
from astropy.table import QTable, Table, Column
from astropy.utils.exceptions import AstropyWarning
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import numpy.core.records as rec

from spectral_cube import SpectralCube
import radio_beam
import aplpy
import astropy.units as u


from RadioAbsTools import cube_tools, spectrum_tools


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

def get_source(file_list, target, folder, selavy_table):
    
    fname = '{}{}_sl.fits'.format(folder, target['comp_name'])
    if not fname in file_list:
        return None
    
    comp_cat_row= selavy_table[selavy_table['component_name']==target['comp_name']][0]
    src = {'ra':comp_cat_row['ra_deg_cont'], 'dec':comp_cat_row['dec_deg_cont'], 
           'a':comp_cat_row['maj_axis']/2, 'b':comp_cat_row['min_axis']/2, 'pa': comp_cat_row['pos_ang'],
          'comp_name': target['comp_name'], 'fname': fname, 'flux_peak': comp_cat_row['flux_peak']}
    return src
    
def save_spectrum(velocity, opacity, flux, filename, 
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

    votable = from_table(table)
    writeto(votable, filename)


def extract_spectrum(fname, src, continuum_start_vel, continuum_end_vel, figures_folder, num_edge_chan=10):

    hdulist = fits.open(fname)
    image = hdulist[0].data
    header = hdulist[0].header
    w = WCS(header)
    index = np.arange(header['NAXIS3'])
    velocities = w.wcs_pix2world(10,10,index[:],0,0)[2]
    radius=math.ceil(max(6, src['a']))

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

def output_spectrum(spec_folder, spectrum_array, opacity, sigma_opacity, comp_name, continuum_start_vel, continuum_end_vel):
    filename = '{}/{}_spec'.format(spec_folder, comp_name)

    save_spectrum(spectrum_array.velocity, opacity, spectrum_array.flux, filename+'.vot', sigma_opacity)
    spectrum_tools.plot_absorption_spectrum(spectrum_array.velocity, opacity, filename+'.png', comp_name, continuum_start_vel, continuum_end_vel, sigma_opacity, range=None)
    spectrum_tools.plot_absorption_spectrum(spectrum_array.velocity, opacity, filename+'_zoom.png', comp_name, continuum_start_vel, continuum_end_vel, sigma_opacity, range=(75,275))



def extract_all_spectra(targets, file_list, cutouts_folder, selavy_table, figures_folder,spectra_folder, max_spectra = 5000):
    src_table = QTable(names=('id', 'comp_name', 'ra', 'dec', 'rating', 'flux_peak', 'mean_cont', 'sd_cont', 'opacity_range', 
            'max_s_max_n', 'min_opacity', 'semi_maj_axis', 'semi_min_axis', 'pa'),
            dtype=('int', 'U32', 'float64', 'float64', 'str', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 
            'float64', 'float64', 'float64'),
            meta={'name': 'ASKAP SMC Spectra'})

    print('Processing {} cutouts into spectra, input:{}'.format(len(targets), cutouts_folder))

    i = 0
    for tgt in targets:
        src = get_source(file_list, tgt, cutouts_folder, selavy_table)
        if not src:
            print('Skipping missing src #{} {}'.format(tgt['id'], tgt['comp_name']))
            continue
        i+=1
        if i > max_spectra:
            print ("Reaching maximum spectra for this run")
            break
            
        print('\nExtracting spectrum for src #{} {} maj:{} min:{} pa:{:.3f}'.format(tgt['id'], tgt['comp_name'], src['a'], src['b'], src['pa']))

        plot_mom0(src['fname'], tgt['comp_name'], figures_folder, src['ra'], src['dec'], 
                src['a'], src['b'], src['pa'])
        
        continuum_start_vel = -100*u.km.to(u.m)
        continuum_end_vel = -60*u.km.to(u.m)
        spectrum = extract_spectrum(src['fname'], src, continuum_start_vel, continuum_end_vel, figures_folder)
        
        mean_cont, sd_cont = spectrum_tools.get_mean_continuum(spectrum.velocity, spectrum.flux, continuum_start_vel, continuum_end_vel)
        opacity = spectrum.flux/mean_cont
        sigma_opacity = sd_cont*np.ones(spectrum.velocity.shape)
        min_opacity = np.min(opacity)
        rating, opacity_range, max_s_max_n = spectrum_tools.rate_spectrum(opacity, sd_cont)

        output_spectrum(spectra_folder, spectrum, opacity, sigma_opacity, src['comp_name'], continuum_start_vel, continuum_end_vel)
        
        src_table.add_row([tgt['id'], tgt['comp_name'], src['ra']*u.deg, src['dec']*u.deg, 
                        rating, src['flux_peak']*u.Jy, mean_cont*u.Jy, sd_cont, opacity_range, max_s_max_n, min_opacity, 
                        src['a'], src['b'], src['pa']])

    return src_table


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
              transform=ax.get_transform('fk5'))
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


def add_sources(ax, sources,detection, best):

    for src in sources:
        colour_name = 'darkgray'
        facecolor = 'darkgray'
        sigma = (1-src['min_opacity'])/src['sd_cont']
        marker = '.'
        if sigma > best:
            marker = 'o'
            facecolor = 'none'
            colour_name = 'red'

        elif sigma > detection:
            marker = '^'
            facecolor = 'none'
            colour_name = 'blue'
            
        elif 'rating' not in src or src['rating'] <= 'B':
            marker = '+'
            facecolor = 'darkgreen'

        ax.scatter([src['ra']], [src['dec']], transform=ax.get_transform('world'), 
               marker=marker, edgecolor=colour_name, facecolor=facecolor)#,
               #s=300, edgecolor=colour_name, facecolor='none', lw=2)


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
    ax.grid()
    return ax, wcs


def plot_source_loc_map(spectra_table, figures_folder, background='hi_zea_ms_mom0.fits'):

    print('\nPlotting {} source locations over background of {}.'.format(len(spectra_table), background))

    fig = plt.figure(figsize=(10.5, 9))
    ax, wcs = plot_background_map(fig, background)

    field_centre = get_field_centre(spectra_table)
    recenter(ax, wcs, field_centre.ra.value, field_centre.dec.value, width=8.25, height=7.25)  # degrees

    # Display the moment map image
    #plt.colorbar(im,fraction=0.046, pad=0.04)
    add_sources(ax, spectra_table, 3, 5)

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

    targets = read_targets('targets_{}.csv'.format(args.sbid))
    
    # Read and filter catalogue
    #src_votable = votable.parse('AS037_Continuum_Component_Catalogue_8906_100.votable', pedantic=False)
    src_votable = votable.parse(args.catalogue, pedantic=False)
    selavy_table = src_votable.get_first_table().to_table()
    
    spectra_table = extract_all_spectra(targets, file_list, cutout_folder, selavy_table, figures_folder,spectra_folder)

    spectra_vo_table = from_table(spectra_table)
    writeto(spectra_vo_table, parent_folder+'askap_spectra.vot')

    plot_source_loc_map(spectra_table, figures_folder)
    plot_field_loc_map(spectra_table, figures_folder)

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Extracted %d spectra in %.02f s' %
          (len(spectra_table), end - start))

    return 0



if __name__ == '__main__':
    exit(main())
    
