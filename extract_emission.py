# Python script to produce emission spectra for a series of sources.
# The script works wityh ASKAP data to produce synthesised beam emission.
# Emission spectra extraction is based on work by Claire Murray

# Author James Dempsey
# Date 7 Dec 2020

import argparse
import csv 
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
    parser.add_argument("-e", "--emission", help="The HI emission cube to source emission data around each source.",
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


def prep_folders(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.mkdir(folder)
            print ("Created " + folder)


def output_emission_spectrum(source, comp_name, velocity, em_mean, em_std, filename):
    title = 'Emission for source #{} {}'.format(source['id'], comp_name)
    em_out_tab = Table(meta={'name': title})
    em_out_tab.add_column(Column(name='velocity', data=velocity, unit='m/s', description='LSRK velocity'))
    em_out_tab.add_column(Column(name='em_mean', data=em_mean, unit='K', description='Mean brightness temperature'))
    em_out_tab.add_column(Column(name='em_std', data=em_std, unit='K', description='Noise level in the brightness temperature'))
    
    votable = from_table(em_out_tab)
    votable.params.append(Param(votable, id='id', value=source['id'], datatype='int'))
    votable.params.append(Param(votable, id='comp_name', value=comp_name, datatype='char', arraysize='*'))
    votable.params.append(Param(votable, id='ra', value=source['ra'], unit='deg', datatype='double'))
    votable.params.append(Param(votable, id='dec', value=source['dec'], unit='deg', datatype='double'))
    writeto(votable, filename)


def extract_channel_slab(filename, chan_start, chan_end):
    cube = SpectralCube.read(filename)
    vel_cube = cube.with_spectral_unit(u.m/u.s, velocity_convention='radio')
    slab = vel_cube[chan_start:chan_end,:, :].with_spectral_unit(u.km/u.s)

    return slab


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
    print (imin, imax, jmin, jmax)
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
    tb_mean = np.nanmean(sub_specx, axis=0)
    # Estimate the uncertainty per channel via the standard deviation over all spectra
    tb_std = np.nanstd(sub_specx, axis=0)
    #print ('mean=',tb_mean)

    return tb_mean, tb_std


def extract_emission_around_source_by_plane(slab, pos, radius_outer, radius_inner):
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
    print (imin, imax, jmin, jmax)
    sub_specx = []
    num_channels = slab.shape[0]
    cutout_center = (radius_outer, radius_outer)
    for plane in range(num_channels):
        cutout = slab[plane, jmin:jmax, imin:imax]
        print (cutout)
        idx = 0
        for k in np.arange(imin, imax):
            for j in np.arange(jmin, jmax):
                kj = np.array([k-imin, j-jmin])
                dist = np.linalg.norm(kj - np.array(cutout_center))
                if dist > radius_inner and dist <= radius_outer:
                    value = slab[plane, kj[1], kj[0]]
                    if plane == 0:
                        spec = np.zeros((num_channels))
                        sub_specx = sub_specx + [spec]
                    spec = sub_specx[idx]
                    spec[plane] = value
                    idx += 1
    
    #print (sub_specx)
    # Compute the mean over all spectra
    tb_mean = np.nanmean(sub_specx, axis=0)
    # Estimate the uncertainty per channel via the standard deviation over all spectra
    tb_std = np.nanstd(sub_specx, axis=0)
    #print ('mean=',tb_mean)

    return tb_mean, tb_std


def extract_emission_spectra(cube, spectra_table, slab_size=160):
    
    # Read the cube
    spec_cube = SpectralCube.read(cube, mode='readonly', memmap=False)
    vel_cube = spec_cube.with_spectral_unit(u.m/u.s, velocity_convention='radio')
    wcs = vel_cube.wcs.celestial
    spec_len =vel_cube.shape[0]
    header = fits.getheader(cube)
    velocities = vel_cube.spectral_axis

    # Identify the target pixels for each spectrum
    pixcoord = wcs.wcs_world2pix(spectra_table['ra'],spectra_table['dec'], 0)

    radius_outer = 8
    radius_inner = 4

    # Extract the spectra
    start = time.time()
    print("  ## Started emission spectra extract at {} ##".format(
          (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))))
    prev = start
    tb_mean_all = []
    tb_std_all = []
    for s in spectra_table:
        tb_mean_all.append(np.zeros(spec_len))
        tb_std_all.append(np.zeros(spec_len))
    
    # Extract using slabs
    unit = None
    prev = time.time()
    for i in range(0,spec_len,slab_size):
        max_idx = min(i+slab_size, spec_len)
        #slab = extract_channel_slab(cube, i, max_idx)
        slab = vel_cube[i:max_idx,:, :].with_spectral_unit(u.km/u.s)

        print (slab)
        unit = slab.unit
        post_slab = time.time()

        for j in range(len(pixcoord[0])):
            pos = [pixcoord[0][j],pixcoord[1][j]]
            tb_mean, tb_std = extract_emission_around_source(slab, pos, radius_outer, radius_inner)
            tb_mean_all[j][i:max_idx] = tb_mean
            tb_std_all[j][i:max_idx] = tb_std
        
        checkpoint = time.time()
        print ("Reading slab of channels {} to {}, took {:.2f} s".format(i, max_idx-1, post_slab-prev))
        print ("Scanning slab of channels {} to {} for {} sources, took {:.2f} s".format(i, max_idx-1, len(pixcoord[0]), checkpoint-post_slab))
        prev = checkpoint

    end = time.time()
    print("  ## Finished emission spectra extract at {}, took {:.2f} s ##".format(
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)), end-start))
    return tb_mean_all, tb_std_all, velocities.value


def output_emission_spectra(spectra_table, tb_mean_all, tb_std_all, velocities, spectra_folder):
    print (velocities)
    for idx, source in enumerate(spectra_table):
        tb_mean = tb_mean_all[idx]
        tb_std = tb_std_all[idx]
        comp_name = source['comp_name']
        if np.sum(tb_mean) == 0:
            print ("Source {} has all no emission data".format(comp_name))
        filename = '{}/{}_emission.vot'.format(spectra_folder, comp_name)
        output_emission_spectrum(source, comp_name, velocities, tb_mean, tb_std, filename)

def main():
    warnings.simplefilter('ignore', category=AstropyWarning)

    args = parseargs()

    start = time.time()
    print("#### Started ASKAP emission spectra extraction for sbid {} at {} ####".format(args.sbid,
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))))

    parent_folder = 'sb{}/'.format(args.sbid)
    if args.parent:
        parent_folder = args.parent + '/'
    if not os.path.exists(parent_folder):
        print("Error: Folder {} does not exist.".format(parent_folder))
        return 1
    if not os.path.exists(args.emission):
        print("Error: File {} does not exist.".format(args.emission))
        return 1

    spectra_folder = parent_folder + 'spectra/'
    prep_folders([spectra_folder])

    targets = read_targets('targets_{}.csv'.format(args.sbid))
    
    # Extract emission spectra
    tb_mean_all, tb_std_all, velocities = extract_emission_spectra(args.emission, targets)
    output_emission_spectra(targets, tb_mean_all, tb_std_all, velocities, spectra_folder)

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Extracted %d spectra in %.02f s' %
          (len(targets), end - start))

    return 0



if __name__ == '__main__':
    exit(main())
    
