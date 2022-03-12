#!/usr/bin/env python -u

# Script to prepare a GASKAP absorption cutout run

# Author James Dempsey
# Date 6 Mar 2021

import argparse
import csv
import difflib
import glob
import math
import os
import shutil
import subprocess
import sys
import time
import warnings
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.io import ascii, votable
from astropy.table import Table, Column
from astropy import units as u
from astropy.utils.exceptions import AstropyWarning
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np

class CommandFailedError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

def parseargs():
    """
    Parse the command line arguments
    :return: An args map with the parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Script to prepare a GASKAP absorption cutout run for an ASKAP scheduing block")
    parser.add_argument("-s", "--sbid", help="The id of the ASKAP scheduling block to be processed",
                        type=str, required=True)
    parser.add_argument("-m", "--ms", help="The path to the measurement sets to be used for the cutouts. If omitted the pre-exiting beam details will be used",
                        type=str, required=False)
    parser.add_argument("-c", "--catalogue", help="The path to the Selevy compnent catalogue for the scheduling block",
                        type=str, required=True)
    parser.add_argument("-p", "--ms_pat", help="The pattern to be used to match the folders of the measurement sets top be used.",
                        type=str, default='*.ms')


    parser.add_argument("--output_folder", help="The output folder which will contain the job configuration files",
                        default=None)
    parser.add_argument("--status_folder", help="The status folder which will contain the job completion or failed files",
                        default='status')

    parser.add_argument("--pbs", help="Run the jobs via PBS qsub command", default=False,
                        action='store_true')
    parser.add_argument("-l", "--log_folder", help="The folder which will contain the stdout and stderr files from the jobs",
                        default='logs')
    args = parser.parse_args()
    return args


def run_os_cmd(cmd, failOnErr=True):
    """
    Run an operating system command ensuring that it finishes successfully.
    If the comand fails, the program will exit.
    :param cmd: The command to be run
    :return: None
    """
    print(">", cmd)
    sys.stdout.flush()
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0:
            message = "Command '"+cmd+"' failed with code " + str(retcode)
            print(message, file=sys.stderr)
            if failOnErr:
                raise CommandFailedError(message)
    except OSError as e:
        message = "Command '" + cmd + "' failed " + e
        print(message, file=sys.stderr)
        if failOnErr:
            raise CommandFailedError(message)
    return None

def rename_columns(table):
    names = np.asarray(table.colnames)
    for name in names:
        if name.startswith('col_'):
            table.rename_column(name, name[4:])

def slice_strings(a,start,end):
    if end is None:
        if start > 0:
            raise('end must be present when start is positive')
        b = a.view((str,1)).reshape(len(a),-1)[:,start:]
        return np.frombuffer(b.tostring(),dtype=(str,start*-1))
    
    b = a.view((str,1)).reshape(len(a),-1)[:,start:end]
    return np.frombuffer(b.tostring(),dtype=(str,end-start))

def get_beams_near_src(target_loc, beam_locs, beams, max_sep = 0.8*(1*u.deg)):
    """
    Find the beams within a certain angular distance of the target location.
    
    Parameters
    ----------
    target_loc: SkyCoord
        The sky location of the target.
    beam_locs: SkyCoord[]
        The locations of each beam.
    beams: str[]
        Array of beam names
    max_sep: dimension
        Maximum distance a beam centre can be from the target to be included (optional).
    
    Returns
    -------
    List of beam names within the requested distance and a list of their distances from the target.
    """
    beam_sep = beam_locs.separation(target_loc)
    beams_covering_target = beam_sep < max_sep
    src_beams = beams[beams_covering_target]
    src_beam_sep = beam_sep[beams_covering_target]
    return src_beams, src_beam_sep


def find_mismatches(file_list):
    """
    Scan a list of file names and build a list of regions of the file names (including paths) that are not the same 
    for all files.
    
    Parameters
    ----------
    file_list: list(str)
        The list of the file paths to be compared
    
    Returns
    -------
    Array of regions in the paths that do not match, each region will be a pair of positions representing the start and 
    end+1 of the range.
    """
    # Scan the list of files and record the positions where each file name differs from the previous file
    mismatches=set()
    prev_value = None
    for file_path in file_list:
        if prev_value != None:
            s = difflib.SequenceMatcher(None, prev_value, file_path)
            for opcode in s.get_opcodes():
                if opcode[0] != 'equal':
                    mismatches.add((opcode[1], opcode[2]))
        prev_value = file_path

    # Identify the unique regions of difference, storing each region once only and ignoring contained sub-regions
    unique_mismatches = []
    for i, region in enumerate(sorted(mismatches)):
        contained = False
        for j, other in enumerate(sorted(mismatches)):
            if i != j:
                if region[0] >= other[0] and region[1] <= other[1]:
                    #print (region,"inside", other)
                    contained = True
        if not contained:
            unique_mismatches.append(region)

    #print (sorted(unique_mismatches))
    return unique_mismatches


def build_ms_pattern(ms_loc, ms_pat):
    """
    Build up a single pattern for filenames of the beam measurement sets with placeholders for the
    beam number and interleave.
    """
    # Create a list the names of the measurement set files
    search_path = '{}/*/{}'.format(ms_loc, ms_pat)
    print ("Find measurement sets matching", search_path)
    file_list = glob.glob(search_path)
    if len(file_list) == 0:
        search_path = '{}/{}'.format(ms_loc, ms_pat)
        print ("Find measurement sets matching", search_path)
        file_list = glob.glob(search_path)

    # Find which parts of the file names vary between files
    mismatches=find_mismatches(file_list)

    # Idenitfy if the mismatches are for the beam number or interleave name 
    # and build up definitions of the varying regions
    regions = []
    for start,end in sorted(mismatches):
        if start==end:
            continue
        if file_list[0][start:end].isnumeric():
            print("beam {} to {} is {}".format(start, end, file_list[0][start:end]))
            regions.append(("{1}", start, end))
        else:
            print("interleave {} to {} is {}".format(start, end, file_list[0][start:end]))
            regions.append(("{0}", start, end))
    regions.reverse()

    # Build a pattern for the measurement sets by replacing the varying regions in a sample path with placeholders
    pattern = str(file_list[0])
    for region in regions:
        pattern = pattern[0:region[1]]+region[0]+pattern[region[2]:]
    return pattern


def record_data_loc(sbid, ms_pat, data_loc_fname):
    index = None
    if os.path.exists(data_loc_fname):
        data_loc = ascii.read(data_loc_fname, guess=False, delimiter=',')    
        for idx, row in enumerate(data_loc):
            if str(row['sbid']) == str(sbid):
                index = idx
    else:
        data_loc = Table(names=('sbid','pattern'), dtype=('i4', 'S500'))
        #data_loc.add_column()
    
    data_loc_new = Table(data_loc)
    all_pat = data_loc['pattern'].data.tolist()
    if index is None:
        data_loc_new.add_row([str(sbid), ''])
        all_pat.append(ms_pat)
    else:
        all_pat[index] = ms_pat
    data_loc_new['pattern'] = Column(data=all_pat, name='pattern')
    if os.path.exists(data_loc_fname):
        shutil.copyfile(data_loc_fname, data_loc_fname+'.bak')
    ascii.write(data_loc_new, output=data_loc_fname, format='csv', overwrite=True)
    print ("Recorded ms pattern of sbid {} as {}".format(sbid, ms_pat))
    

def get_target_list(catalogue, flux_peak_min=15):
    # Read and filter catalogue
    src_votable = votable.parse_single_table(catalogue, pedantic=False)
    table = src_votable.to_table()
    rename_columns(table)
    targets = table[table['flux_peak']>flux_peak_min]

    # Filter out non-unique component names
    names, counts = np.unique(targets['component_name'], return_counts=True)
    duplicates = names[counts > 1]
    for comp_name in duplicates:
        indexes = targets['component_name'] == comp_name
        max_peak_flux = np.max(targets['flux_peak'][indexes])
        to_remove = (targets['component_name'] == comp_name) & (targets['flux_peak'] < max_peak_flux)
        idx_to_rem = np.arange(len(targets))[to_remove]
        print ('Ignoring weaker duplicate sources', targets['component_id'][idx_to_rem].data )
        targets.remove_rows(idx_to_rem)

    print ("Found {} targets from {} sources".format(len(targets), len(table)))
    return targets


def generate_beam_listing(sb_folder, sbid):
    me = Path(__file__)
    script = Path(me.parent, 'list_beams.sh')

    run_os_cmd('{} {}'.format(script, sbid))
    beam_listing = '{0}/beam_listing_SB{1}.csv'.format(sb_folder, sbid)
    return beam_listing


def find_beam_locs(beam_listing):
    # Build the list of beam locations
    beams = ascii.read(beam_listing, format='no_header', guess=False, delimiter=',')
    names = ('col1','col2', 'col3', 'col4')
    new_names = ('filename','field', 'ra_rad', 'dec_rad')
    beams.rename_columns(names, new_names)

    mismatches=find_mismatches(beams['filename'])
    for region in mismatches:
        if region[1]-region[0] ==1:
            fn_il_loc = region
        else:
            fn_beam_loc = region
    beam_id = slice_strings(beams['filename'], fn_beam_loc[0], fn_beam_loc[1])
    #print (beam_id)
    #interleave_id = slice_strings(beams['interleave'], -1, None)# beams['interleave'][:][-1:]

    mismatches=find_mismatches(beams['field'])
    il_loc = list(mismatches)[0]
    interleave_id = slice_strings(beams['field'], il_loc[0], il_loc[1])# beams['interleave'][:][-1:]

    file_interleave = slice_strings(beams['filename'], fn_il_loc[0], fn_il_loc[1])

    ids = np.stack((beam_id, interleave_id), axis=-1)
    unique_ids, unique_idx = np.unique(ids, axis=0, return_index=True)

    beams['beam_id'] = beam_id
    beams['interleave'] = interleave_id
    beams['file_interleave'] = file_interleave
    unique_beams = beams[beams['interleave'] == beams['file_interleave']]

    u_beam_locs = SkyCoord(ra=unique_beams['ra_rad']*u.rad, dec=unique_beams['dec_rad']*u.rad, frame='icrs')
    return unique_beams, u_beam_locs


def plot_beams_and_targets(sb_folder, targets, beams, beam_locs):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    im = ax.scatter(targets['ra_deg_cont'], targets['dec_deg_cont'], c=targets['flux_peak'], cmap='RdBu_r', vmin=15, vmax=100)
    for il in 'ABC':
        filt = beams['interleave'] == il
        ax.scatter(beam_locs.ra.value[filt], beam_locs.dec.value[filt], marker='+')
        # Plot footprints for one interleave
        if il == 'A':
            for loc in beam_locs[filt]:
                c=Ellipse((loc.ra.value, loc.dec.value), width=0.8*2/math.cos(loc.dec.rad), height=0.8, fill=False, ls=':', zorder=-1)
                ax.add_artist(c)

    cb = fig.colorbar(im)
    cb.set_label('Peak flux (mJy/beam)')
    ax.grid()
    plt.gca().invert_xaxis()
    plt.xlabel('Right Ascension')
    plt.ylabel('Declination')
    plt.savefig(sb_folder+'/sources_beams.png', bbox_inches='tight')


def get_image_params_table(dest_folder, sbid, targets, beams, beam_locs, max_sep_close=0.55*u.deg, max_sep_far=0.8*u.deg):
    comp_names = []
    comp_ra = []
    comp_dec = []
    included_beam_nums = []
    included_beam_interleaves = []
    included_beam_ids = []
    included_beam_sep = []

    for tgt in targets:
        target_loc = SkyCoord(ra=tgt['ra_deg_cont']*u.degree, dec=tgt['dec_deg_cont']*u.degree, frame='icrs')
        src_beams, src_beam_sep = get_beams_near_src(target_loc, beam_locs, beams, max_sep=max_sep_close)
        if len(src_beams) == 0:
            #print ("No beams close to {}, checking out to {}".format(tgt['component_name'], max_sep_far))
            src_beams, src_beam_sep = get_beams_near_src(target_loc, beam_locs, beams, max_sep=max_sep_far)

        for i in range(len(src_beams)):
            comp_names.append(tgt['component_name'])
            comp_ra.append(tgt['ra_deg_cont'])
            comp_dec.append(tgt['dec_deg_cont'])
            included_beam_nums.append(src_beams['beam_id'].data[i])
            included_beam_interleaves.append(src_beams['interleave'].data[i])
            included_beam_ids.append(src_beams['beam_id'].data[i]+src_beams['interleave'].data[i])
            included_beam_sep.append(src_beam_sep.to(u.deg).value[i])

    image_params = Table()
    image_params['component_name'] = comp_names
    image_params['comp_ra'] = comp_ra
    image_params['comp_dec'] = comp_dec
    image_params['beam_nums'] = included_beam_nums
    image_params['beam_interleaves'] = included_beam_interleaves
    image_params['beam_ids'] = included_beam_ids
    image_params['beam_sep'] = included_beam_sep

    image_params_vot = votable.from_table(image_params)
    filename = "{}/sb{}_srcs_image_params.vot".format(dest_folder, sbid)
    votable.writeto(image_params_vot, filename)

    print ("Produced VO table file {} with {} target beam combos.".format(filename, len(image_params)))

    return image_params


def report_beam_usage(image_params):
    ar = np.array(image_params['beam_ids'])
    for i in range(36):
        for interleave in ('A', 'B', 'C'):
            key = '{:02d}{}'.format(i, interleave)
            count = len(ar[ar==key])
            if count == 0:
                print ("Warning: Beam {} is unused!".format(key))
    #print ("Total instances of beam usage {}".format(len(ar)))
    print ('Mean targets per beam {:.2f}'.format(len(ar)/(36*3)))
    mean_beams_per_source = len(ar) / len(np.unique(image_params['component_name']))
    print ('Mean beams per target {:.2f}'.format(mean_beams_per_source))


def create_targets_csv(dest_folder, sbid, targets, image_params):
    # Create the csv file of targets and included beams
    csv_filename = '{}/targets_{}.csv'.format(dest_folder, sbid)
    with open(csv_filename, 'w') as csvfile:
        tgtwriter = csv.writer(csvfile, delimiter=',',
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)
        tgtwriter.writerow(['index', 'component_name', 'ra', 'dec', 'beams'])
        i = 1
        for tgt in targets:
            comp_name = tgt['component_name']
            row = []
            row.append(str(i))
            row.append(comp_name)
            row.append(tgt['ra_deg_cont'])
            row.append(tgt['dec_deg_cont'])
            for tgt_beam in image_params:
                if comp_name == tgt_beam['component_name']:
                    row.append(tgt_beam['beam_ids'])
            tgtwriter.writerow(row)
            i+= 1

    mean_beams_per_source = len(image_params) / len(np.unique(image_params['component_name']))
    print ('Produced csv file {} with {} targets (mean {:.2f} beams/target).'.format(csv_filename, i-1, mean_beams_per_source))


def main():
    # Parse command line options
    args = parseargs()

    start = time.time()
    print("#### Started preparing GASKAP absorption cutout run for sbid {} at {} ####".format
          (args.sbid, time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))))

    # Check locations
    if args.ms and not os.path.exists(args.ms):
        print("Could not find measurement sets at", args.ms)
        return 1
    if not os.path.exists(args.catalogue):
        print("Could not find component catalogue at", args.catalogue)
        return 1
    sb_folder = 'sb{}'.format(args.sbid)
    if not os.path.exists(sb_folder):
        os.makedirs(sb_folder)
    beam_listing = '{0}/beam_listing_SB{1}.csv'.format(sb_folder, args.sbid)
    if not args.ms:
        if not os.path.exists(beam_listing):
            print("No measurement sets supplied and beam listing {} does not existm", beam_listing)
            return 1
        print ('Using pre-existing beam information at', beam_listing)

    # Setup the data_loc file
    warnings.simplefilter('ignore', category=AstropyWarning)
    if args.ms:
        ms_patterns = build_ms_pattern(args.ms, args.ms_pat)
        record_data_loc(args.sbid, ms_patterns, 'data_loc.csv')
    

    # Build the beam listing
    if args.ms:
        beam_listing = generate_beam_listing(sb_folder, args.sbid)
    beams, beam_locs = find_beam_locs(beam_listing)

    # Build the target list
    targets = get_target_list(args.catalogue)
    plot_beams_and_targets(sb_folder, targets, beams, beam_locs)
    num_targets = len(targets)
    image_params = get_image_params_table(sb_folder, args.sbid, targets, beams, beam_locs)
    report_beam_usage(image_params)
    create_targets_csv(sb_folder, args.sbid, targets, image_params)

    # Report
    end = time.time()
    print('#### Processing completed at {} ####'.format(
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))))
    print('Processed {0} targets in {1:.0f} sec'.format(
          num_targets, (end - start)))

    return 0


if __name__ == '__main__':
    exit(main())