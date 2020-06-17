# Python script to produce spectra for a series of sub-cubes.
# The spectra are assumed to be absoprtion spectra and rated according to their suitability for absorption detection.

# Author James Dempsey
# Date 19 Feb 2020

import argparse
import csv 
import glob
import os
import time
import warnings

import astropy.units as u
from astropy.table import Table, Column
from astropy.io import ascii
from astropy.io.ascii import Csv
from astropy.io.votable.tree import Param,Info
from astropy.io.votable import from_table, writeto
from astropy.io import fits, votable
from astropy.wcs import WCS
from astropy.table import QTable
from astropy.utils.exceptions import AstropyWarning
import matplotlib.pyplot as plt
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


def output_spectra(velocity, opacity, flux, filename, 
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
    table.add_column(Column(name='velocity', data=velocity, unit='m/s'))
    table.add_column(Column(name='opacity', data=opacity))
    table.add_column(Column(name='flux', data=flux, unit='Jy', description='Flux per beam'))
    table.add_column(Column(name='sigma_opacity', data=sigma_tau, description='Noise in the absorption profile'))

    votable = from_table(table)
    writeto(votable, filename)


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
    a_deg = src_maj*u.arcsec.to(u.deg)
    b_deg = src_min*u.arcsec.to(u.deg)
    fig.show_ellipses(src_ra, src_dec, a_deg, b_deg, angle=src_pa, edgecolor='red') #, 

    figname = '{}/{}_mom0.'.format(out_folder, comp_name)
    fig.savefig(figname+'.pdf')
    fig.savefig(figname+'.png')
    fig.close()

def get_source(file_list, target, folder, selavy_table):
    
    fname = '{}{}_sl.fits'.format(folder, target['comp_name'])
    if not fname in file_list:
        return None
    
    comp_cat_row= selavy_table[selavy_table['component_name']==target['comp_name']][0]
    src = {'ra':comp_cat_row['ra_deg_cont'], 'dec':comp_cat_row['dec_deg_cont'], 
           'a':comp_cat_row['maj_axis'], 'b':comp_cat_row['min_axis'], 'pa': comp_cat_row['pos_ang'],
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
    table.add_column(Column(name='velocity', data=velocity, unit='m/s'))
    table.add_column(Column(name='opacity', data=opacity))
    table.add_column(Column(name='flux', data=flux, unit='Jy', description='Flux per beam'))
    table.add_column(Column(name='sigma_opacity', data=sigma_tau, description='Noise in the absorption profile'))

    votable = from_table(table)
    writeto(votable, filename)


def extract_spectrum(fname, src, continuum_start_vel, continuum_end_vel, num_edge_chan=10):

    hdulist = fits.open(fname)
    image = hdulist[0].data
    header = hdulist[0].header
    w = WCS(header)
    index = np.arange(header['NAXIS3'])
    velocities = w.wcs_pix2world(10,10,index[:],0,0)[2]

    img_slice = cube_tools.get_integrated_spectrum(image, w, src, velocities, continuum_start_vel, continuum_end_vel, radius=6)

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
    src_table = QTable(names=('id', 'comp_name', 'ra', 'dec', 'rating', 'flux_peak', 'mean_cont', 'sd_cont', 'opacity_range', 'max_s_max_n', 'min_opacity'),
            dtype=('int', 'U32', 'float64', 'float64', 'str', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64'),
                meta={'name': 'ASKAP SMC Spectra'})

    print('Processing {} cutouts into spectra'.format(len(targets)))

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
            
        print('\nExtracting spectrum for src #{} {}'.format(tgt['id'], tgt['comp_name']))

        plot_mom0(src['fname'], tgt['comp_name'], figures_folder, src['ra'], src['dec'], 
                src['a'], src['b'], src['pa'])
        
        continuum_start_vel = -100*u.km.to(u.m)
        continuum_end_vel = -60*u.km.to(u.m)
        spectrum = extract_spectrum(src['fname'], src, continuum_start_vel, continuum_end_vel)
        
        mean_cont, sd_cont = spectrum_tools.get_mean_continuum(spectrum.velocity, spectrum.flux, continuum_start_vel, continuum_end_vel)
        opacity = spectrum.flux/mean_cont
        sigma_opacity = sd_cont*np.ones(spectrum.velocity.shape)
        min_opacity = np.min(opacity)
        rating, opacity_range, max_s_max_n = spectrum_tools.rate_spectrum(opacity, sd_cont)

        output_spectrum(spectra_folder, spectrum, opacity, sigma_opacity, src['comp_name'], continuum_start_vel, continuum_end_vel)
        
        src_table.add_row([tgt['id'], tgt['comp_name'], src['ra']*u.deg, src['dec']*u.deg, 
                        rating, src['flux_peak']*u.Jy, mean_cont*u.Jy, sd_cont, opacity_range, max_s_max_n, min_opacity])

    return src_table


def prep_folders(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.mkdir(folder)
            print ("Created " + folder)
               

def main():
    warnings.simplefilter('ignore', category=AstropyWarning)

    args = parseargs()

    start = time.time()
    print("#### Started ASKAP spectra extraction for sbid {} at {} ####".format(args.sbid,
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))))

    parent_folder = 'sb{}/'.format(args.sbid)
    if args.parent:
        parent_folder = args.parent
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

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Extracted %d spectra in %.02f s' %
          (len(spectra_table), end - start))

    return 0



if __name__ == '__main__':
    exit(main())
    
