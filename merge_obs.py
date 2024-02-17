# Python script to merge absorption data from multiple observations. This will produce a consolidated set of absorption 
# spectra organised by field. Where a source has been observed more than once the cubelets for the observations will be 
# stacked (weighted by noise) to produce a single cubelet and a consolidated spectrum will be retrieved from the 
# combined subcube.
#
# Note that to run this you will need to have access to all subcubes for the absorption spectra to be combined.

# Author James Dempsey
# Date 09 Jul 2023


import argparse
import glob
import os
import re
import shutil
import time
import warnings

from datetime import datetime, timezone

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from astropy.coordinates import SkyCoord, FK5, search_around_sky
from astropy.io import ascii, fits
from astropy.io.votable import from_table, parse_single_table, writeto
from astropy.table import QTable, Table, Column, hstack, vstack
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS
from regions import Regions
from scipy import interpolate

import extract_spectra


DATA_RELEASE_VERSION=0.3

class DuplicateSource:
    def __init__(self):
        self.field = ''
        self.comp_name = 0
        self.sbids = []
        self.comp_names = []
        self.sd_conts = []
    
    def __init__(self, field, comp_name, sbids, comp_names, sd_conts):
        self.field = field
        self.comp_name = comp_name
        self.sbids = sbids
        self.comp_names = comp_names
        self.sd_conts = sd_conts

    def __str__(self):
        return "Duplicate: " + str(self.__dict__)
    
    def __repr__(self):
        return "<Duplicate field:%s comp_name:%s>" % (self.field, self.comp_name)


def parseargs():
    """
    Parse the command line arguments
    :return: An args map with the parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Merge GASKAP absorption observations.")
    parser.add_argument("input_folder", help="The folder where the source absoprtion data can be found.")
    parser.add_argument("output_folder", help="The folder where the consolidated spectra will be written.")
    parser.add_argument("field_list", help="The file listing the GASKAP fields. This should be a DS9 region file " + 
                        "with text naming the fields.")
    parser.add_argument("-s", "--sbids", help="The file where the list of sbids can be found. This should be a csv file " + 
                        "with an sbid column.")
    parser.add_argument("-r", "--release_version", help="The file where the list of sbids can be found. This should be a csv file " + 
                        "with an sbid column.", default=DATA_RELEASE_VERSION, required=False)
    parser.add_argument("--averaging", help="How many channels should be averaged together, default is no averaging", type=int,
                        default=4, required=False)
            

    args = parser.parse_args()
    return args

def find_poloygon_centre(region):
    # TODO: this is not robust and should be replaced with a library function
    vertices = region.vertices
    ras = []
    decs = []
    for coord in vertices:
        ras.append(coord.ra.deg)
        decs.append(coord.dec.deg)
    centre_ra = np.mean(ras)
    centre_dec = np.mean(decs)
    #print (centre_ra, centre_dec)
    return centre_ra*u.deg, centre_dec*u.deg

def read_fields(field_list):
    regions = Regions.read(field_list, format='ds9')
    #print (regions[1].visual)
    #print (regions[1].meta)
    #print (regions[1].serialize(format='ds9'))
    #print (regions[0].center)

    fld_ids = []
    fld_names = []
    fld_cras = []
    fld_cdecs = []
    fld_regions = []
    fields = []
    for idx in range(0,len(regions),2):
        text = regions[idx+1].meta['label'].split(' ')
        c_ra, c_dec = find_poloygon_centre(regions[idx])
        fld_ids.append(text[0])
        fld_names.append(text[1])
        fld_cras.append(c_ra)
        fld_cdecs.append(c_dec)
        fld_regions.append(regions[idx])
    fields = QTable([fld_ids, fld_names, fld_cras, fld_cdecs, fld_regions], 
                    names=('id', 'name', 'centre_ra', 'centre_dec', 'region'),)
    print (fields[0])
    print ("Read {} GASKAP fields".format(len(fields)))

    return fields


def prep_folders(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.mkdir(folder)
            print ("Created " + folder)

def log_config(args):
    print ('Configuration:')
    print (' {: >20} : {}'.format('input_folder', args.input_folder))
    print (' {: >20} : {}'.format('output_folder', args.output_folder))
    print (' {: >20} : {}'.format('field_list', args.field_list))
    print (' {: >20} : {}'.format('sbids', args.sbids))
    print (' {: >20} : {}'.format('data release version', args.release_version))
    print (' {: >20} : {}'.format('averaging', args.averaging))
    print ('')

def check_config(args):
    if not os.path.exists(args.input_folder):
        print("Error: Folder {} does not exist.".format(args.input_folder))
        return False
    if not os.path.exists(args.field_list):
        print("Error: Field file {} does not exist.".format(args.field_list))
        return False
    if args.sbids and not os.path.exists(args.sbids):
        print("Error: Sbid file {} does not exist.".format(args.sbids))
        return False
    if not os.path.exists(args.output_folder):
        print("Error: Folder {} does not exist.".format(args.output_folder))
        return False
    return True
    
def find_sbids(input_folder, allowed_sbids_file):
    pattern = '{}/sb*'.format(input_folder)
    folders = glob.glob(pattern)
    all_sbids = []
    sbid_pat = re.compile("^.*/sb([0-9]+)$")
    for folder in folders:
        match = re.search(sbid_pat, folder)
        if match:
            sbid = match.group(1)
            if not os.path.exists(os.path.join(folder, "cutouts")):
                print ("Skipping folder {} as it does not have cutouts".format(folder))
                continue
            spec_pat = '{}/**/spectra'.format(folder, sbid)
            if len(glob.glob(spec_pat)) == 0:
                print ("Skipping folder {} as it does not have spectra".format(folder))
                continue
            cat_pat = '{}/averaged/gaskap_sb{}_abs_spectra.vot'.format(folder, sbid)
            if len(glob.glob(cat_pat)) == 0:
                print ("Skipping folder {} as it does not have a spectra catalogue".format(folder))
                continue
            if not os.path.exists(os.path.join(folder, "spectra")):
                print ("WARNING: Folder {} does not have full resolution spectra".format(folder))
            if not os.path.exists(os.path.join(folder, "averaged", "spectra")):
                print ("WARNING: Folder {} does not have averaged spectra".format(folder))

            all_sbids.append(sbid)

    # Filter by sbids allowed llist
    if allowed_sbids_file:
        allowed_sbids_table = ascii.read(allowed_sbids_file, format="csv")
        allowed_sbids =  allowed_sbids_table[allowed_sbids_table.colnames[0]]
        found_sbids = all_sbids
        all_sbids = np.intersect1d(allowed_sbids, all_sbids)
        if (len(found_sbids) > len(all_sbids)):
            print ("Found the following sbids: ", found_sbids)
            print("Skipped {} sbids that were not in the allowed list.".format(len(found_sbids)-len(all_sbids)))
        
    print ("Processing sbids:", all_sbids)
    return all_sbids


def prep_field_folders(fields, output_folder):
    folders = []
    for field in fields:
        field_folder = os.path.join(output_folder, field)
        folders.append(field_folder)
        folders.append(os.path.join(field_folder, 'cutouts'))
        folders.append(os.path.join(field_folder, 'spectra'))
        folders.append(os.path.join(field_folder, 'figures'))
        avg_folder = os.path.join(field_folder, 'averaged')
        folders.append(avg_folder)
        folders.append(os.path.join(avg_folder, 'spectra'))
        folders.append(os.path.join(avg_folder, 'figures'))
    
    prep_folders(folders)


def copy_source_data(src_folder, dest_folder, comp_name):
    file_pat = "{}/*{}*".format(src_folder, comp_name)
    file_list = glob.glob(file_pat)
    for f in file_list:
        shutil.copy2(f, dest_folder)
    return len(file_list)


def first_src_better(src1, src2):
    if src1['rating'] != src2['rating']:
        return src1['rating'] < src2['rating']
    else:
        return src1['sd_cont'] <= src2['sd_cont']
    
def plot_field_distribution(uniq_fields, spec_avg_dat, sbid, output_folder):
    #colours = np.full((len(spec_dat), 3), sns.color_palette()[0])
    fig, ax = plt.subplots(1,1)
    for idx, field in enumerate(uniq_fields):
        fld_dat = spec_avg_dat[spec_avg_dat['field'] == field]
        plt.scatter(fld_dat['ra']/15, fld_dat['dec'], color=sns.color_palette()[idx], label=field)
    ax.invert_xaxis()
    ax.grid(linestyle=':')
    ax.set_xlabel("Right Ascension (hours)")
    ax.set_ylabel("Declination (degrees)")
    ax.legend(bbox_to_anchor=(1.05, 1),
                         loc='upper left', borderaxespad=0.)
    plt.title("Field distribution for sbid {}".format(sbid))
    plt.tight_layout()
    plt.savefig('{}/sb{}_fields.png'.format(output_folder, sbid))


def add_unique_spectra_to_global_catalogue(spec_dat, idxes_new, sbid, abs_dat, global_spec_catalogue, global_abs_catalogue):
    ids_to_eliminate = spec_dat['id'][idxes_new]
    dupe_mask = np.isin(spec_dat['id'], ids_to_eliminate)
    spec_to_add = spec_dat[~dupe_mask]
    #all_sbids = [[int(sbid)] for i in range(len(spec_to_add))]
    all_sbids = np.full((len(spec_to_add)), sbid, dtype=object)
    spec_to_add['all_sbids'] = all_sbids
    
    dupe_abs_mask = np.isin(abs_dat['src_id'], ids_to_eliminate)
    abs_to_add = abs_dat[~dupe_abs_mask]

    global_spec_catalogue = vstack([global_spec_catalogue, spec_to_add])
    global_abs_catalogue = vstack([global_abs_catalogue, abs_to_add])
    return global_spec_catalogue, global_abs_catalogue


def distribute_obs_spectra(sbid, input_folder, output_folder, fields, field_centres, global_avg_spec_catalogue, 
                           global_avg_abs_catalogue, global_native_spec_catalogue, global_native_abs_catalogue, 
                           duplicate_srcs):
    # Assumption: Both the native and averaged catalogues have the same number of entries and the same component names
    # (i.e. they are baed on the same continuum catalogue) but they do not need to be in the same order.

    # Read averaged spectra catalogue
    sb_path = os.path.join(input_folder, "sb"+sbid)
    avg_path = os.path.join(sb_path, "averaged")
    print (avg_path)
    spec_cat_votable = parse_single_table(os.path.join(avg_path, 'gaskap_sb{}_abs_spectra.vot'.format(sbid)))
    spec_avg_dat = spec_cat_votable.to_table()
    abs_cat_votable = parse_single_table(os.path.join(avg_path, 'gaskap_sb{}_absorption.vot'.format(sbid)))
    abs_avg_dat = abs_cat_votable.to_table()
    #spec_dat['sbid'] = np.full((len(spec_dat)), str(sbid), dtype=object)
    #print (spec_avg_dat.columns)
    #print (spec_dat)

    # Read native spectra catalogue
    spec_cat_votable = parse_single_table(os.path.join(sb_path, 'gaskap_sb{}_abs_spectra.vot'.format(sbid)))
    spec_native_dat = spec_cat_votable.to_table()
    abs_cat_votable = parse_single_table(os.path.join(sb_path, 'gaskap_sb{}_absorption.vot'.format(sbid)))
    abs_native_dat = abs_cat_votable.to_table()
    if len(spec_native_dat) != len(spec_avg_dat):
        raise ValueError("SBID " + sbid + " had different numbers of native and averaged spectra")

    # Determine fields for each spectrum
    src_pos = SkyCoord(ra=spec_avg_dat['ra'], dec=spec_avg_dat['dec'])
    idx, d2d, d3d = src_pos.match_to_catalog_sky(field_centres)
    spec_avg_dat['field'] = fields['name'][idx]
    all_sbids = np.full((len(spec_avg_dat)), sbid, dtype=object)
    spec_avg_dat['all_sbids'] = all_sbids

    src_pos = SkyCoord(ra=spec_native_dat['ra'], dec=spec_native_dat['dec'])
    idx, d2d, d3d = src_pos.match_to_catalog_sky(field_centres)
    spec_native_dat['field'] = fields['name'][idx]
    all_sbids = np.full((len(spec_native_dat)), sbid, dtype=object)
    spec_native_dat['all_sbids'] = all_sbids

    plot_field_distribution(np.unique(spec_avg_dat['field']), spec_avg_dat, sbid, output_folder)

    # Determine which sources are unique in the global set so far and update the global set
    if global_avg_spec_catalogue is None:
        global_avg_spec_catalogue = spec_avg_dat
        spec_to_add = spec_avg_dat
        #all_sbids = [[int(sbid)] for i in range(len(global_spec_catalogue))]
        all_sbids = np.full((len(spec_avg_dat)), sbid, dtype=object)
        global_avg_spec_catalogue['all_sbids'] = all_sbids
        #print (global_spec_catalogue['all_sbids'][0])
        #global_spec_catalogue['all_sbids'][0] = np.concatenate(( global_spec_catalogue['all_sbids'][0], np.array([12345])))
        #print (global_spec_catalogue)
        global_avg_abs_catalogue = abs_avg_dat

        global_native_spec_catalogue = spec_native_dat
        all_sbids = np.full((len(spec_native_dat)), sbid, dtype=object)
        global_native_spec_catalogue['all_sbids'] = all_sbids
        global_native_abs_catalogue = abs_native_dat
    else:
        # Identify sources already in catalogue
        pos_exist = SkyCoord(global_avg_spec_catalogue['ra'], global_avg_spec_catalogue['dec'])
        pos_new_avg = SkyCoord(spec_avg_dat['ra'], spec_avg_dat['dec'])
        idxes_new_avg, idxes_exist_avg, sep, _ = search_around_sky(pos_new_avg, pos_exist, 3*u.arcsec)
        spec_to_add = spec_avg_dat[idxes_new_avg]
        print ("Found {} averaged spectra that are already in the combined set.".format(len(idxes_new_avg)))
        print(idxes_new_avg)
        print(idxes_exist_avg)
        pos_exist = SkyCoord(global_native_spec_catalogue['ra'], global_native_spec_catalogue['dec'])
        pos_new_native = SkyCoord(spec_native_dat['ra'], spec_native_dat['dec'])
        idxes_new_native, idxes_exist_native, sep, _ = search_around_sky(pos_new_native, pos_exist, 3*u.arcsec)
        print ("Found {} native spectra that are already in the combined set.".format(len(idxes_new_native)))
        print(idxes_new_native)
        print(idxes_exist_native)

        # Record details of duplicate sources for later stacking
        for i, idx in enumerate(idxes_exist_avg):
            exist_avg_spec_record = global_avg_spec_catalogue[idx]
            new_avg_spec_record = spec_avg_dat[idxes_new_avg[i]]
            #new_native_spec_record = spec_native_dat[idxes_new_native[i]]
            new_native_idx = np.where(spec_native_dat['comp_name'] == new_avg_spec_record['comp_name'])[0][0]
            print (new_native_idx)
            new_native_spec_record = spec_native_dat[new_native_idx]
            #native_idx = idxes_exist_native[i]
            native_idx = np.where(global_native_spec_catalogue['comp_name'] == exist_avg_spec_record['comp_name'])[0][0]
            exist_native_spec_record = global_native_spec_catalogue[native_idx]
            #print(exist_spec_record['all_sbids'])
            #all_sbids = global_spec_catalogue['all_sbids'][idx].tolist()
            #all_sbids.append(sbid)
            global_avg_spec_catalogue['all_sbids'][idx] = global_avg_spec_catalogue['all_sbids'][idx] + "," + sbid
            global_native_spec_catalogue['all_sbids'][native_idx] = global_native_spec_catalogue['all_sbids'][native_idx] + "," + sbid
            # We make the decision on which spectrum is best based on the avergaed spectrum only
            exist_src_better = first_src_better(exist_avg_spec_record, new_avg_spec_record)
            first_comp_name = exist_avg_spec_record['comp_name']
            if first_comp_name in duplicate_srcs:
                duplicate_rec = duplicate_srcs[first_comp_name]
                if not exist_src_better:
                    duplicate_rec.comp_name = new_avg_spec_record['comp_name']
                    del duplicate_srcs[first_comp_name]
                    duplicate_srcs[new_avg_spec_record['comp_name']] = duplicate_rec

            else:
                first_comp_name = exist_avg_spec_record['comp_name'] if exist_src_better else new_avg_spec_record['comp_name']
                duplicate_rec = DuplicateSource(exist_avg_spec_record['field'], first_comp_name, 
                                          [exist_avg_spec_record['sbid']], [exist_avg_spec_record['comp_name']], 
                                          [exist_avg_spec_record['sd_cont']])
                duplicate_srcs[first_comp_name] = duplicate_rec
            duplicate_rec.sbids.append(new_avg_spec_record['sbid'])
            duplicate_rec.comp_names.append(new_avg_spec_record['comp_name'])
            duplicate_rec.sd_conts.append(new_avg_spec_record['sd_cont'])

            if not exist_src_better:
                # Replace row in existing table
                print ("Replaced {} - {} with {} - {}".format(exist_avg_spec_record['sbid'], exist_avg_spec_record['comp_name'],
                                                              new_avg_spec_record['sbid'], new_avg_spec_record['comp_name']))
                print ("     and {} - {} with {} - {}".format(exist_native_spec_record['sbid'], exist_native_spec_record['comp_name'],
                                                              new_native_spec_record['sbid'], new_native_spec_record['comp_name']))
                new_avg_spec_record['all_sbids'] = global_avg_spec_catalogue['all_sbids'][idx]
                global_avg_spec_catalogue[idx] = new_avg_spec_record
                new_native_spec_record['all_sbids'] = global_native_spec_catalogue['all_sbids'][native_idx]
                global_native_spec_catalogue[native_idx] = new_native_spec_record
    
        # Update global catalogue 
        global_avg_spec_catalogue, global_avg_abs_catalogue = add_unique_spectra_to_global_catalogue(spec_avg_dat, idxes_new_avg, 
                                                                                      sbid, abs_avg_dat, global_avg_spec_catalogue, 
                                                                                      global_avg_abs_catalogue)
        global_native_spec_catalogue, global_native_abs_catalogue = add_unique_spectra_to_global_catalogue(spec_native_dat, idxes_new_native, 
                                                                                      sbid, abs_native_dat, global_native_spec_catalogue, 
                                                                                      global_native_abs_catalogue)

    # Copy cubelets and spectra for non-duplicate sources to fields
    if True:
        print ("Distributing {} unique spectra for sbid {}".format(len(spec_to_add), sbid))
        num_files = 0
        prep_field_folders(np.unique(spec_avg_dat['field']), output_folder)
        for s in spec_to_add:
            for folder in ('cutouts', 'spectra', 'figures', 'averaged/spectra', 'averaged/figures'):
                src_folder = os.path.join(sb_path, folder)
                dest_folder = os.path.join(output_folder, s['field'], folder)
                num_files += copy_source_data(src_folder, dest_folder, s['comp_name'])
        print ("Distributed {} files".format(num_files))

    return global_avg_spec_catalogue, global_avg_abs_catalogue, global_native_spec_catalogue, global_native_abs_catalogue


def get_velocity_axis(hdu):
    wcs1 = WCS(hdu[0].header)
    spec_idx1 = np.arange(0,wcs1.array_shape[-3])
    velocity1 = wcs1.spectral.all_pix2world(spec_idx1, 0)[0]
    return velocity1


def stack_cutouts(src_to_merge, input_folder, output_folder):
    velocity_comb = None
    hdus = []
    weights = []

    # Read all of the cutouts, calculate the weights, calculate the intersection of the velocity axes of all input 
    # cutouts. We order by increasing noise to avoid interpolation of the best spectrum
    noise_order = np.argsort(src_to_merge.sd_conts)
    print ("Got noise order for {} of {}".format(src_to_merge.comp_name, np.array(src_to_merge.sd_conts)[noise_order]))
    for idx, sbid in enumerate(np.array(src_to_merge.sbids)[noise_order]):
        src_idx = noise_order[idx]
        comp_name = src_to_merge.comp_names[src_idx]
        sb_path = os.path.join(input_folder, "sb{}".format(sbid))
        cutout_path = os.path.join(sb_path, 'cutouts', '{}_sl.fits'.format(comp_name))
        print (" {}", cutout_path)
        single_hdu = fits.open(cutout_path)
        hdus.append(single_hdu)
        weights.append(1/src_to_merge.sd_conts[src_idx]**2) # Weight by inverse square of the noise
        curr_velocity = get_velocity_axis(single_hdu)
        if velocity_comb is None:
            velocity_comb = curr_velocity
            #print (velocity_comb)
        else:
            vel_intersect_mask = (velocity_comb >= np.min(curr_velocity)) & (velocity_comb <= np.max(curr_velocity)) #| (np.isclose(velocity_comb,  np.min(velocity1))) | (np.isclose(velocity_comb,  np.max(velocity1)))
            velocity_comb = velocity_comb[vel_intersect_mask]
            #print (velocity_comb)

    # Normalise the weights
    weights = np.array(weights)
    weights = weights/np.nansum(weights)

    # Start with the data from the best cube trimmed to the combined velocity axis and weighted
    hdu1 = hdus[0]
    fitsdata1 = hdu1[0].data
    velocity1 = get_velocity_axis(hdu1)
    vel_intersect_mask = (velocity1 >= np.min(velocity_comb)) & (velocity1 <= np.max(velocity_comb))
    fitsdata_comb = fitsdata1[:,vel_intersect_mask,:,:]*weights[0]

    # Add the weighted data from the rest of the cubes
    for idx in range(1, len(src_to_merge.sd_conts)):
        hdu2 = hdus[idx]
        fitsdata2 = hdu2[0].data
        velocity2 = get_velocity_axis(hdu2)
        for i in range(fitsdata2.shape[-1]):
            for j in range(fitsdata2.shape[-2]):
                f = interpolate.interp1d(velocity2, fitsdata2[0, :, j, i], fill_value=0, bounds_error=False, kind='cubic')
                regridded_spec = f(velocity_comb)
                fitsdata_comb[0,:,j,i]+=regridded_spec*weights[idx]

    # Note the data is already normalised as the weights are normalised
    #print ("Renormalising data by ", np.nansum(weights))
    #fitsdata_comb /= np.nansum(weights)

    # Prepare the FITS header for the stacked cutout
    comp_name = src_to_merge.comp_names[0]
    field_path = os.path.join(output_folder, src_to_merge.field)
    cutout_path = os.path.join(field_path, 'cutouts', '{}_sl.fits'.format(comp_name))
    wcs1 = WCS(hdu1[0].header)
    spec_col = wcs1.world_axis_physical_types.index('spect.dopplerVeloc.radio')
    wcs1.wcs.crpix[spec_col] = 1
    wcs1.wcs.crval[spec_col] = velocity_comb[0]
    header=wcs1.to_header()
    cutout_header = hdu1[0].header.copy()
    for key in header:
        if key not in ('COMMENT', 'HISTORY'):
            cutout_header.set(key, value=header[key], comment=header.comments[key])

    cutout_header['OBJECT'] = (comp_name, 'Name of the target GASKAP absorption component')
    cutout_header.insert('OBJECT', ('FIELD', src_to_merge.field, 'GASKAP field identifier'), after=True)
    cutout_header.add_history('Merged by GASKAP HI Absorption Pipeline merge_obs.py', after=-1)
    cutout_header.add_history(' on {:%Y-%m-%dT%H:%M:%S}'.format(datetime.now(timezone.utc)), after=-1)
    cutout_header.add_comment('Included ASKAP sbids are:', after=-1)
    for idx, sbid in enumerate(src_to_merge.sbids):
        cutout_header.add_comment('sbid {} comp {} sd_cont {:.4f}'.format(sbid, src_to_merge.comp_names[idx], src_to_merge.sd_conts[idx]), after=-1)

    # Save the stacked cutout
    new_hdu = fits.PrimaryHDU(data=fitsdata_comb, header=cutout_header)
    new_hdul = fits.HDUList(new_hdu)
    for old_hdu in hdu1[1:]:
        if type(old_hdu) is fits.BinTableHDU and 'EXTNAME' in old_hdu.header and str(old_hdu.header['EXTNAME']).strip() == 'BEAMS':
            # Trim the beams table to match the trimmed spectral axis
            old_hdr = old_hdu.header
            old_data = old_hdu.data[vel_intersect_mask]
            old_hdr['NAXIS2'] = old_data.shape[0]
            new_hdul.append(fits.BinTableHDU(header=old_hdr, data=old_data))
        else:
            old_hdr = old_hdu.header
            new_hdul.append(old_hdu)
    new_hdul.writeto(cutout_path, overwrite=True)
    print ("Wrote stacked cube to", cutout_path)
    return cutout_path


def check_milky_way(ra, dec):
    loc = SkyCoord(ra=ra, dec=dec, frame='icrs')
    is_mw = np.abs(np.min(loc.galactic.b.value)) < 20
    return is_mw


def add_to_table(target, to_be_added):
    for row in to_be_added:
        target.add_row(row)

def assess_spectrum(tgt, src_def, abs_spec, is_milky_way, cont_range, comp_cat_row):
    src_table, abs_table = extract_spectra.define_spectra_tables("")
    extract_spectra.assess_single_spectrum(tgt, src_def, abs_spec, src_table, abs_table, 0, is_milky_way, False, cont_range)
    src_table.add_column(Column(name='field', data=comp_cat_row['field'][0]))
    src_table.add_column(Column(name='all_sbids', data=comp_cat_row['all_sbids'][0], dtype=object))
    return src_table, abs_table


def merge_spectra(src_to_merge, global_avg_spec_catalogue, global_avg_abs_catalogue, global_native_spec_catalogue, 
                  global_native_abs_catalogue, input_folder, output_folder, 
                  averaging, cont_range=(-110,-70), weighting='square'):
    print ('   ---    ')
    print ("Merging {} cutouts for {}".format(len(src_to_merge.sbids), src_to_merge.comp_name))

    # Stack the cutouts
    cutout_path = stack_cutouts(src_to_merge, input_folder, output_folder)

    # Save merged cubelets to fields
    comp_cat_row = global_avg_spec_catalogue[global_avg_spec_catalogue['comp_name'] == src_to_merge.comp_name]
    #print ("ra type is ", type(comp_cat_row['ra']))
    src_def = {'ra':comp_cat_row['ra'].quantity[0], 'dec':comp_cat_row['dec'].quantity[0], 
           'a':comp_cat_row['semi_maj_axis'].quantity[0], 'b':comp_cat_row['semi_min_axis'].quantity[0], 'pa': comp_cat_row['pa'].quantity[0],
          'comp_name': src_to_merge.comp_name, 'fname': cutout_path, 'flux_peak': comp_cat_row['flux_peak'].quantity[0],
           'flux_int': comp_cat_row['flux_int'].quantity[0], 'component_id': comp_cat_row['component_id'][0]}
    tgt = {'ra':comp_cat_row['ra'].quantity[0], 'dec':comp_cat_row['dec'].quantity[0],'id':0, 'comp_name': src_to_merge.comp_name}
    is_milky_way = check_milky_way(tgt['ra'], tgt['dec'])
    
    figures_folder = os.path.join(output_folder, src_to_merge.field, 'figures')
    spectra_folder = os.path.join(output_folder, src_to_merge.field, 'spectra')
    avg_figures_folder = os.path.join(output_folder, src_to_merge.field, 'averaged', 'figures')
    avg_spectra_folder = os.path.join(output_folder, src_to_merge.field, 'averaged', 'spectra')
    #print (figures_folder)

    extract_spectra.plot_mom0(cutout_path, src_to_merge.comp_name, figures_folder, [src_def])

    # Recategorise the merged spectra
    continuum_start_vel = cont_range[0]*u.km.to(u.m)
    continuum_end_vel = cont_range[1]*u.km.to(u.m)
    native_spectrum = extract_spectra.extract_spectrum(cutout_path, src_def, continuum_start_vel, continuum_end_vel, 
                                                figures_folder, weighting)
    extract_spectra.process_spectrum(native_spectrum, tgt, spectra_folder, 
                                                       src_to_merge.comp_name, 
                                                       continuum_start_vel, continuum_end_vel, "")
    native_abs_spec_filename = '{}/{}_spec.vot'.format(spectra_folder, src_to_merge.comp_name)
    native_abs_spec_votable = parse_single_table(native_abs_spec_filename)
    native_abs_spec = native_abs_spec_votable.to_table()
        
    avg_spectrum = extract_spectra.average_spectrum(native_spectrum, averaging)

    src_to_merge.stacked_sd_cont = extract_spectra.process_spectrum(avg_spectrum, tgt, avg_spectra_folder, 
                                                       src_to_merge.comp_name, 
                                                       continuum_start_vel, continuum_end_vel, "")

    # Read back in the processed spectrum 
    avg_abs_spec_filename = '{}/{}_spec.vot'.format(avg_spectra_folder, src_to_merge.comp_name)
    avg_abs_spec_votable = parse_single_table(avg_abs_spec_filename)
    avg_abs_spec = avg_abs_spec_votable.to_table()

    print ("field {} source {} stacked noise {:.6f}".format(src_to_merge.field, src_to_merge.comp_name, src_to_merge.stacked_sd_cont))
    # Assess the combined spectra
    native_src_table, native_abs_table = assess_spectrum(tgt, src_def, native_abs_spec, is_milky_way, cont_range, comp_cat_row)
    avg_src_table, avg_abs_table = assess_spectrum(tgt, src_def, avg_abs_spec, is_milky_way, cont_range, comp_cat_row)

    stacked_noise = avg_src_table['sd_cont'].data[0]
    #print ("Unstacked noise {:.4f} and stacked noise {:.4f}".format(comp_cat_row['sd_cont'][0],stacked_noise))
    if stacked_noise > comp_cat_row['sd_cont'][0]:
        stacked = False
        print ("WARN: Ignoring stacking for {} as the noise is worse. {:.4f} vs {:.4f}".format(src_to_merge.comp_name, comp_cat_row['sd_cont'][0], stacked_noise))
    else:
        stacked = True
        # Replace the old entries for this source with the combined spectrum entries
        src_filter = global_avg_spec_catalogue['comp_name'] == src_to_merge.comp_name
        global_avg_spec_catalogue = global_avg_spec_catalogue[~src_filter]
        abs_filter = global_avg_abs_catalogue['comp_name'] == src_to_merge.comp_name
        global_avg_abs_catalogue = global_avg_abs_catalogue[~abs_filter]
        print (global_avg_spec_catalogue.colnames)
        add_to_table(global_avg_spec_catalogue, avg_src_table)
        add_to_table(global_avg_abs_catalogue, avg_abs_table)

        src_filter = global_native_spec_catalogue['comp_name'] == src_to_merge.comp_name
        global_native_spec_catalogue = global_native_spec_catalogue[~src_filter]
        abs_filter = global_native_abs_catalogue['comp_name'] == src_to_merge.comp_name
        global_native_abs_catalogue = global_native_abs_catalogue[~abs_filter]
        print (global_native_spec_catalogue.colnames)
        add_to_table(global_native_spec_catalogue, native_src_table)
        add_to_table(global_native_abs_catalogue, native_abs_table)

        #global_spec_catalogue = vstack([global_spec_catalogue, src_table])
        #global_abs_catalogue = vstack([global_abs_catalogue, abs_table])
        print ("Replaced {} spectrum and {} abs detections with {} spectrum and {} abs detections for {}".format(
            np.sum(src_filter), np.sum(abs_filter), len(avg_src_table), len(avg_abs_table), src_to_merge.comp_name))
    #print ("Geronimo!!")
    #exit(1)
    return global_avg_spec_catalogue, global_avg_abs_catalogue, global_native_spec_catalogue, global_native_abs_catalogue, stacked

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
    add_col_metadata(spectra_vo_table, 'id', 'Unique identifier of the source.', ucd='meta.id')
    #add_col_metadata(spectra_vo_table, 'sbid', 'The id of the scheduling block id in which the source was observed.', ucd='meta.id')
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
    add_col_metadata(spectra_vo_table, 'peak_tau', 'The maximum tau value in the spectrum, representing the peak absorption.')
    add_col_metadata(spectra_vo_table, 'e_peak_tau', 'The uncertainty in tau at the peak absorption velocity.')
    add_col_metadata(spectra_vo_table, 'has_mw_abs', 'Flag to indicate that the spectrum has absorption in the Milky Way velocity range (<50 km/s).', datatype='boolean')
    add_col_metadata(spectra_vo_table, 'has_other_abs', 'Flag to indicate that the spectrum has absorption away from the Milky Way velocity range (>50 km/s).', datatype='boolean')
    add_col_metadata(spectra_vo_table, 'semi_maj_axis', 'HWHM major axis before deconvolution.', units='arcsec', ucd='phys.angSize.smajAxis;em.radio;stat.fit')
    add_col_metadata(spectra_vo_table, 'semi_min_axis', 'HWHM minor axis before deconvolution.', units='arcsec', ucd='phys.angSize.sminAxis;em.radio;stat.fit')
    add_col_metadata(spectra_vo_table, 'pa', 'Position angle before deconvolution measured as degrees CCW of North.', units='deg', ucd='phys.angSize;pos.posAng;em.radio;stat.fit')
    add_col_metadata(spectra_vo_table, 'n_h', 'Column density towards the background source excluding Milky Way velocities, from the GASS survey.', units='cm-2', ucd='phys.columnDensity')
    add_col_metadata(spectra_vo_table, 'noise_flag', 'Flag to indicate that the noise for the spectrum may be underestimated.', datatype='boolean')
    add_col_metadata(spectra_vo_table, 'continuum_slope', 'The slope of the continuum, shown as change in relative absorption per 10 km/s. Expected to be close to zero.')
    add_col_metadata(spectra_vo_table, 'all_sbids', 'The list of scheduling blocks in which this source was observed.')


def add_absorption_column_metadata(abs_vo_table):
    add_col_metadata(abs_vo_table, 'src_id', 'Identifier of the source.', ucd='meta.id')
    add_col_metadata(abs_vo_table, 'comp_name', 'Background source component name.', ucd='meta.id')
    add_col_metadata(abs_vo_table, 'abs_name', 'Absorption feature name.', ucd='meta.id;meta.main')
    add_col_metadata(abs_vo_table, 'ra', 'J2000 right ascension in decimal degrees.', units='deg', ucd='pos.eq.ra;meta.main')
    add_col_metadata(abs_vo_table, 'dec', 'J2000 declination in decimal degrees.', units='deg', ucd='pos.eq.dec;meta.main')
    add_col_metadata(abs_vo_table, 'rating', 'Quality rating of the absorption spectrum.')
    add_col_metadata(abs_vo_table, 'flux_peak', 'Peak flux density of the background source.', units='mJy/beam', ucd='phot.flux.density;stat.max;em.radio;stat.fit')
    add_col_metadata(abs_vo_table, 'mean_cont', 'Mean continuum level per channel of the spectrum.', units='mJy/beam', ucd='phot.flux.density;stat.mean;em.radio;stat.fit')
    add_col_metadata(abs_vo_table, 'sd_cont', '1-sigma noise level of the spectrum in optical depth units. Does not include emission noise.')
    add_col_metadata(abs_vo_table, 'opacity_range', 'The range in the absorption spectrum from highest emission to lowest absorption.')
    add_col_metadata(abs_vo_table, 'max_s_max_n', 'Ratio of maximum absorption versus maximum emission noise in the spectrum.')
    add_col_metadata(abs_vo_table, 'max_noise', 'The maximum emission noise in the spectrum, in optical depth units.')
    add_col_metadata(abs_vo_table, 'num_chan_noise', 'The number of channels in the spectrum that have emission above the 1-sigma noise level.')
    add_col_metadata(abs_vo_table, 'semi_maj_axis', 'HWHM major axis before deconvolution.', units='arcsec', ucd='phys.angSize.smajAxis;em.radio;stat.fit')
    add_col_metadata(abs_vo_table, 'semi_min_axis', 'HWHM minor axis before deconvolution.', units='arcsec', ucd='phys.angSize.sminAxis;em.radio;stat.fit')
    add_col_metadata(abs_vo_table, 'pa', 'Position angle before deconvolution measured as degrees CCW of North.', units='deg', ucd='phys.angSize;pos.posAng;em.radio;stat.fit')
    add_col_metadata(abs_vo_table, 'start_vel', 'The lowest velocity covered by this absorption feature.', units='km/s', ucd='phys.velocity;stat.min')
    add_col_metadata(abs_vo_table, 'end_vel', 'The highest velocity covered by this absorption feature.', units='km/s', ucd='phys.velocity;stat.max')
    add_col_metadata(abs_vo_table, 'length', 'The number of channels included in this absorption feature.', units='chan')
    add_col_metadata(abs_vo_table, 'min_optical_depth', 'The minimum optical depth in the feature, in optical depth units, representing the peak absorption.')
    add_col_metadata(abs_vo_table, 'e_min_optical_depth', 'The uncertainty in the optical depth at the peak absorption velocity.')
    add_col_metadata(abs_vo_table, 'peak_tau', 'The maximum tau value in the feature, representing the peak absorption.')
    add_col_metadata(abs_vo_table, 'e_peak_tau', 'The uncertainty in tau at the peak absorption velocity.')
    add_col_metadata(abs_vo_table, 'max_sigma', 'The maximum significance of the absorption feature, representing the peak absorption.', ucd='stat.snr')
    add_col_metadata(abs_vo_table, 'ew', 'The integral of absorption (tau) across the feature.', units='km/s', ucd='spect.line.eqWidth;em.radio')
    add_col_metadata(abs_vo_table, 'e_ew', 'The uncertainty in the integral of absorption (tau) across the feature.', units='km/s', ucd='stat.error;spect.line.eqWidth;em.radio')

def write_spectra_votable(spectra_table, filename):
    spectra_vo_table = from_table(spectra_table)
    add_spectra_column_metadata(spectra_vo_table)
    writeto(spectra_vo_table, filename)
    print ("\nWrote {} spectra to {}".format(len(spectra_table), filename))


def write_absorption_votable(abs_table, filename):
    abs_vo_table = from_table(abs_table)
    add_absorption_column_metadata(abs_vo_table)
    writeto(abs_vo_table, filename)
    print ("Wrote {} absorption features to {}".format(len(abs_table), filename))


def main():
    warnings.simplefilter('ignore', category=AstropyWarning)

    args = parseargs()

    start = time.time()
    print("#### Started merging GASKAP absorption observations at {} ####".format(
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))))

    log_config(args)
    prep_folders([args.output_folder])
    if not check_config(args):
        return 1
    
    fields = read_fields(args.field_list)
    field_centres = SkyCoord(ra=fields['centre_ra'], dec=fields['centre_dec'], frame=FK5)

    all_sbids = find_sbids(args.input_folder, args.sbids)
    if len(all_sbids) == 0:
        print ("ERROR: No scheduling blocks to process")
        return 1

    # Read spectra for each sbid
    global_avg_spec_catalogue = None
    global_avg_abs_catalogue = None
    global_native_spec_catalogue = None
    global_native_abs_catalogue = None
    duplicate_srcs = {}
    for sbid in all_sbids:
        global_avg_spec_catalogue, global_avg_abs_catalogue, global_native_spec_catalogue, global_native_abs_catalogue = distribute_obs_spectra(
                                        sbid, args.input_folder, args.output_folder, fields, field_centres, 
                                        global_avg_spec_catalogue, global_avg_abs_catalogue, 
                                        global_native_spec_catalogue, global_native_abs_catalogue, duplicate_srcs)
        
    # Merge spectra across multiple sbids
    print ("\n\nMerging {} sources with multiple observations".format(len(duplicate_srcs)) )
    if len(duplicate_srcs) > 0:
        print(duplicate_srcs)
        sample = list(duplicate_srcs.values())[0]
        print (sample.field, sample.comp_name, sample.sbids, sample.comp_names, sample.sd_conts)
    num_stacked = 0
    merged_stats = Table(names=('field','comp_name', 'sd_cont', 'num_spec', 'min_input_sd_cont', 'max_input_sd_cont'), 
                        dtype=('U32', 'U32', 'float64', 'int', 'float64', 'float64'))
    for src_to_merge in duplicate_srcs.values():
        global_avg_spec_catalogue, global_avg_abs_catalogue, global_native_spec_catalogue, global_native_abs_catalogue, stacked = merge_spectra(
            src_to_merge, global_avg_spec_catalogue, global_avg_abs_catalogue, 
            global_native_spec_catalogue, global_native_abs_catalogue, 
            args.input_folder, args.output_folder, args.averaging)
        merged_stats.add_row([src_to_merge.field, src_to_merge.comp_name, src_to_merge.stacked_sd_cont, len(src_to_merge.sd_conts),
                        np.min(src_to_merge.sd_conts), np.max(src_to_merge.sd_conts)])
        if stacked:
            num_stacked += 1
    ascii.write(merged_stats, output=args.output_folder+'/merged_stats.csv', format='csv', overwrite=True)

    # Save the global catalogues
    num_spectra = len(global_avg_spec_catalogue)
    global_avg_spec_catalogue.remove_column('sbid')
    filename = os.path.join(args.output_folder, 'gaskap_abs_spectra_v{}.vot'.format(args.release_version))
    write_spectra_votable(global_avg_spec_catalogue, filename)
    filename = os.path.join(args.output_folder, 'gaskap_absorption_v{}.vot'.format(args.release_version))
    write_absorption_votable(global_avg_abs_catalogue, filename)
    global_native_spec_catalogue.remove_column('sbid')
    filename = os.path.join(args.output_folder, 'gaskap_native_abs_spectra_v{}.vot'.format(args.release_version))
    write_spectra_votable(global_native_spec_catalogue, filename)
    filename = os.path.join(args.output_folder, 'gaskap_native_absorption_v{}.vot'.format(args.release_version))
    write_absorption_votable(global_native_abs_catalogue, filename)

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Merged %d sbids containing %d spectra including %d stacked spectra (from %d) in %.02f s' %
          (len(all_sbids), num_spectra, num_stacked, len(duplicate_srcs), end - start))

    return 0





if __name__ == '__main__':
    exit(main())
