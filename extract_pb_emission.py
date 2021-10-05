# Python script to produce emission spectra for a series of sources.
# The script works with GASS data to produce primary beam emission.
# Emission spectra extraction is based on work by Claire Murray

# Author James Dempsey
# Date 15 Dec 2020

import argparse
import csv 
import glob
import os
import time
import warnings

import astropy.units as u
from astropy.io import fits, votable
from astropy.io.votable.tree import Param,Info
from astropy.io.votable import from_table, parse_single_table, writeto
from astropy.table import Table, Column
from astropy.utils.exceptions import AstropyWarning
import numpy as np

from spectral_cube import SpectralCube


def parseargs():
    """
    Parse the command line arguments
    :return: An args map with the parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Produce emission spectra for a set of sources")
    parser.add_argument("-s", "--sbid", help="The id of the ASKAP scheduling block to be processed",
                        type=int, required=True)
    parser.add_argument("-c", "--catalogue", help="The catalgoue of source positions and characteristics",
                        required=False)
    parser.add_argument("-e", "--emission_folder", help="The folder in which GASS ZEA data can be found.",
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


def extract_targets(selavy_table):
    table = Table()
    table.add_column(Column(name='id', data=range(len(selavy_table))))
    table.add_column(Column(name='comp_name', data=selavy_table['component_name']))
    table.add_column(Column(name='ra', data=selavy_table['ra_deg_cont']))
    table.add_column(Column(name='dec', data=selavy_table['dec_deg_cont']))
    #table.add_column(Column(name='beams', data=beams))

    return table


def rename_columns(table):
    names = np.asarray(table.colnames)
    for name in names:
        if name.startswith('col_'):
            table.rename_column(name, name[4:])


def prep_folders(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.mkdir(folder)
            print ("Created " + folder)


def output_emission_spectrum(source, comp_name, velocity, em_mean, em_std, filename):
    title = 'Primary beam emission for source #{} {}'.format(source['id'], comp_name)
    em_out_tab = Table(meta={'name': title})
    em_out_tab.add_column(Column(name='velocity', data=velocity, unit='m/s', description='LSRK velocity'))
    em_out_tab.add_column(Column(name='pb_em_mean', data=em_mean, unit='K', description='Mean brightness temperature in the primary beam'))
    em_out_tab.add_column(Column(name='pb_em_std', data=em_std, unit='K', description='Noise level in the brightness temperature in the primary beam'))
    
    votable = from_table(em_out_tab)
    votable.params.append(Param(votable, id='id', value=source['id'], datatype='int'))
    votable.params.append(Param(votable, id='comp_name', value=comp_name, datatype='char', arraysize='*'))
    votable.params.append(Param(votable, id='ra', value=source['ra'], unit='deg', datatype='double'))
    votable.params.append(Param(votable, id='dec', value=source['dec'], unit='deg', datatype='double'))
    writeto(votable, filename)


def extract_emission_around_source(slab, pos, radius_outer, radius_inner):
    xp = np.int(pos[0])
    yp = np.int(pos[1])

    # Only pull spectra whose positions fall on the footprint of the cube (this should be all, right?)
    if (xp < radius_outer) or (xp > slab.shape[2]-radius_outer) or (yp < radius_outer) or (yp > slab.shape[1]-radius_outer):
        print ("Skipping")
        empty_result = np.zeros(slab.shape[0])
        return empty_result, empty_result

    # Define pixel coordinates of a grid surrounding each target
    center = (xp, yp)
    imin = center[0] - radius_outer
    imax = center[0] + radius_outer + 1
    jmin = center[1] - radius_outer
    jmax = center[1] + radius_outer + 1

    # loop through and pile in all spectra which fall in the annulus, based on their distance
    # from the central pixel
    #print (imin, imax, jmin, jmax)
    sub_specx = []
    for k in np.arange(imin, imax):
        for j in np.arange(jmin, jmax):
            kj = np.array([k, j])
            dist = np.linalg.norm(kj - np.array(center))
            if dist > radius_inner and dist <= radius_outer:
                spec = slab[:, kj[1], kj[0]]
                sub_specx = sub_specx + [spec]
    
    #print (sub_specx)
    # Compute the mean over all spectra
    tb_mean = np.nanmedian(sub_specx, axis=0)
    # Estimate the uncertainty per channel via the standard deviation over all spectra
    tb_std = np.nanstd(sub_specx, axis=0)
    #print ('mean=',tb_mean)

    return tb_mean, tb_std


def calc_channels(gass_files):
    prev_max = None
    channel_count = 0
    slab_min_idx = []
    for filename in gass_files:
        print (filename)
        cube = SpectralCube.read(filename)
        mid_x = cube.shape[-1]//2
        mid_y = cube.shape[-2]//2

        spectral_axis = cube.world[:, mid_x, mid_y][0].ravel()
        #print (spectral_axis)
        min_idx = 0
        if prev_max is not None:
            min_idx = np.argmax((spectral_axis-1*(u.m/u.s)) > prev_max)
        print ("Slicing {}:-1 which is {:.3f} to {:.3f}".format(min_idx, spectral_axis[min_idx], spectral_axis[-1]))
        prev_max = spectral_axis[-1]
        channel_count += len(spectral_axis)-min_idx
        slab_min_idx.append(min_idx)    
    
    return channel_count, slab_min_idx


def calc_pixcoords(gass_files, selavy_table):
    cube = SpectralCube.read(gass_files[0])
    pix_coord = cube.wcs.celestial.wcs_world2pix(selavy_table['ra'],selavy_table['dec'], 0)
    return pix_coord


def extract_emission_spectra(gass_files, pixcoord, channel_count, slab_min_idx):
    start = time.time()
    prev = start
    velocities = np.zeros(channel_count)
    num_spectra = len(pixcoord[0])
    all_em_spectra = np.zeros((num_spectra, channel_count))
    all_em_noise = np.zeros((num_spectra, channel_count))
    radius_outer = 7
    radius_inner = 1

    max_chan = 0
    for idx, filename in enumerate(gass_files):
        print (filename, slab_min_idx[idx])
        cube = SpectralCube.read(filename)
        mid_x = cube.shape[-1]//2
        mid_y = cube.shape[-2]//2
        spectral_axis = cube.world[:, mid_x, mid_y][0].ravel()
        min_slab_chan = slab_min_idx[idx]
        min_chan = max_chan
        max_chan = min_chan+len(spectral_axis)-min_slab_chan
        velocities[min_chan:max_chan] = spectral_axis[min_slab_chan:].to(u.km/u.s)

        with fits.open(filename) as hdul:
            slab = hdul[0].data
            post_slab = time.time()
            for j in range(len(pixcoord[0])):
                pos = [pixcoord[0][j],pixcoord[1][j]]
                #if j > 2:
                #    break
                tb_mean, tb_std = extract_emission_around_source(slab, pos, radius_outer, radius_inner)
                all_em_spectra[j][min_chan:max_chan] = tb_mean[min_slab_chan:]
                all_em_noise[j][min_chan:max_chan] = tb_std[min_slab_chan:]
                        
        checkpoint = time.time()
        print ("Reading file {}, took {:.2f} s".format(filename, post_slab-prev))
        print ("Scanning GASS slice for {} sources, took {:.2f} s".format(len(pixcoord[0]), checkpoint-post_slab))
        prev = checkpoint

    end = time.time()
    print("  ## Finished emission spectra extract at {}, took {:.2f} s ##".format(
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)), end-start))
    return all_em_spectra, all_em_noise, velocities


def output_emission_spectra(spectra_table, tb_mean_all, tb_std_all, velocities, spectra_folder):
    print (velocities)
    for idx, source in enumerate(spectra_table):
        tb_mean = tb_mean_all[idx]
        tb_std = tb_std_all[idx]
        comp_name = source['comp_name']
        if np.sum(tb_mean) == 0:
            print ("Source {} has all no emission data".format(comp_name))
        filename = '{}/{}_pb_emission.vot'.format(spectra_folder, comp_name)
        output_emission_spectrum(source, comp_name, velocities, tb_mean, tb_std, filename)

def main():
    warnings.simplefilter('ignore', category=AstropyWarning)

    args = parseargs()

    start = time.time()
    print("#### Started GASS emission spectra extraction for sbid {} at {} ####".format(args.sbid,
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))))

    parent_folder = 'sb{}/'.format(args.sbid)
    if args.parent:
        parent_folder = args.parent + '/'
    if not os.path.exists(parent_folder):
        print("Error: Folder {} does not exist.".format(parent_folder))
        return 1
    if not os.path.exists(args.emission_folder):
        print("Error: Folder {} does not exist.".format(args.emission_folder))
        return 1

    figures_folder = parent_folder + 'figures/'
    spectra_folder = parent_folder + 'spectra/'
    prep_folders([spectra_folder])

    targets_fname = '{}/targets_{}.csv'.format(parent_folder, args.sbid)
    if os.path.exists(targets_fname):
        targets = read_targets(targets_fname)
    else:
        src_votable = votable.parse(args.catalogue, pedantic=False)
        selavy_table = src_votable.get_first_table().to_table()
        rename_columns(selavy_table)
        targets = extract_targets(selavy_table)

    gass_files = glob.glob(args.emission_folder + '/gass_*.zea.fits')
    gass_files = sorted(gass_files)

    # Extract emission spectra
    channel_count, slab_min_idx = calc_channels(gass_files)
    pixcoord = calc_pixcoords(gass_files, targets)
    #tb_mean_all, tb_std_all, velocities = extract_emission_spectra(args.emission, targets)
    tb_mean_all, tb_std_all, velocities = extract_emission_spectra(gass_files, pixcoord, channel_count, slab_min_idx)
    output_emission_spectra(targets, tb_mean_all, tb_std_all, velocities *1000, spectra_folder)

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Extracted %d spectra in %.02f s' %
          (len(targets), end - start))

    return 0



if __name__ == '__main__':
    exit(main())
    
