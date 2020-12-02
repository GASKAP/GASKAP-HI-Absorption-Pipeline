# GenerateSpectraPages.py

import argparse
import os
import time

from astropy.coordinates import SkyCoord
from astropy.io import ascii, votable
from astropy.io.votable import from_table, writeto
from astropy.table import Column, Table
import astropy.units as u
import numpy as np


def parseargs():
    """
    Parse the command line arguments
    :return: An args map with the parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Produce preview pages for a set of spectra")
    parser.add_argument("-s", "--sbid", help="The id of the ASKAP scheduling block",
                        type=int, required=True)
    parser.add_argument("-g", "--good", help="The sigma threshold for spectra to be included in the detections.html page",
                        type=float, default=3.0)
    parser.add_argument("-b", "--best", help="The sigma threshold for spectra to be included in the best.html page",
                        type=float, default=5.0)
    parser.add_argument("-p", "--parent", help="The parent folder for the processing, will default to sbnnn/ where nnn is the sbid.",
                        required=False)
    args = parser.parse_args()
    return args

def output_header(f, title):
    f.write('<!doctype html>\n<html lang="en">\n<head>\n<title>{}</title>'.format(title))
    with open('style.html') as style:
        f.write(style.read())
    f.write('\n</head>\n<body>')
    f.write('\n<div class="container-fluid">')
    f.write('\n<h1 align="middle">{}</h1>'.format(title))
    return

def output_location_plots(f):
    f.write('\n<div class="row px-3" id="maps">')
    f.write('\n<div class="col-md-auto"><h2 class="d-inline font-weight-light text-center text-lg-left mt-4 mb-0">Location</hs></div>')
    f.write('\n<div class="col-md-auto">')
    f.write('\n<a href="figures/field_loc.png" class="d-block mb-4 h-100"  data-lightbox="maps">')
    f.write('\n<img class="img-fluid img-thumbnail" style="height: 180px" src="figures/field_loc.png" alt="Map of the location of the field.">')
    f.write('\n</a>\n</div>')
    f.write('\n<div class="col-md-auto">')
    f.write('\n<a href="figures/source_loc.png" class="d-block mb-4 h-100"  data-lightbox="maps">')
    f.write('\n<img class="img-fluid img-thumbnail" style="height: 180px" src="figures/source_loc.png" alt="Map of the location of the sources.">')
    f.write('\n</a>\n</div>')
    f.write('\n</div>')

def output_block_title(f, rating, first, count):
    if not first:
        f.write('\n\n</div><br/>\n')
    spec = 'spectrum' if count == 1 else 'spectra'
    title = '{} Rating {} {}'.format(count, rating, spec) if rating else '{} Missed {} (with closest source)'.format(count, spec)
    f.write('\n<div>')
    f.write('\n<div class="col-9 d-inline"><h2 class="d-inline font-weight-light text-center text-lg-left mt-4 mb-0">{}</h2></div>'.format(title))
    f.write('\n<div class="col-3 pull-right d-inline"><a class="btn btn-primary" data-toggle="collapse" href="#spectra{0}" role="button" aria-expanded="false" aria-controls="spectra{0}" style="font-size: x-small;">Hide/Show</a></div>'.format(rating))
    f.write('\n</div>')
    f.write('\n<div class="row text-center text-lg-left collapse show" id="spectra{}">'.format(rating))

    
def output_img(f, comp_name, rating):
    f.write('\n<div class="col-lg-3 col-md-4 col-6 px-2">')
    f.write('\n<a href="spectra/{0}_spec.png" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(comp_name, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="spectra/{0}_spec_zoom.png" alt="Zoomed preview of spectrum at {0}">'.format(comp_name))
    f.write('\n</a>\n</div>')
    return

    
def output_non_zoom_img(f, comp_name, rating):
    f.write('\n<div class="col-lg-3 col-md-4 col-6 px-2">')
    f.write('\n<a href="spectra/{0}_spec.png" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(comp_name, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="spectra/{0}_spec.png" alt="Preview of spectrum at {0}">'.format(comp_name))
    f.write('\n</a>\n</div>')
    return

def output_footer(f):
    f.write('\n\n</div>\n</div>\n</body>\n</html>')
    return

def output_j19_img(f, gaskap_name, j19_name, rating, sep=None):
    name_text = gaskap_name
    if sep:
        name_text += ' at {:.1f} arcsec'.format(sep)
    f.write('\n<div class="col-4">')
    f.write('\n<a href="spectra/{0}_spec.png" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(gaskap_name, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="spectra/{0}_spec_zoom.png" alt="Zoomed preview of spectrum at {0}">'.format(gaskap_name))
    f.write('\n{0}</a>\n</div>'.format(name_text))
    f.write('\n<div class="col-8">')
    j19_filename = '../jameson2019figset2/{}_lr.jpg'.format(j19_name)
    f.write('\n<a href="{0}" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(j19_filename, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="{0}" alt="Zoomed preview of spectrum at {0}">'.format(j19_filename))
    f.write('\n</a>\n</div>')

    return


def output_spectra(sbid, table, title, filename, threshold=None, verbose=False):
    print (title, filename)
    with open(filename, 'w') as f:
        output_header(f, title)

        output_location_plots(f)

        for rating in 'ABCDEF':
            targets = table[table['rating']==rating]
            if threshold:
                targets = targets[(1-targets['min_opacity'])/targets['sd_cont'] > threshold]
            comp_names = sorted(targets['comp_name'])
            print('Rating {} has {} spectra'.format(rating, len(comp_names)))

            if verbose:
              print (comp_names)
              
            output_block_title(f, rating, rating=='A', len(comp_names))

            for name in comp_names:
                output_img(f, name, rating)
                    
        output_footer(f)

def find_mw_vel_comp_names(comp_names, parent_folder, threshold, max_velocity):
    mw_comp_names = []
    for name in comp_names:
        filename = '{0}/spectra/{1}_spec.vot'.format(parent_folder, name)
        spectrum_votable = votable.parse(filename, pedantic=False)
        spectrum_tab = spectrum_votable.get_first_table().to_table()
        non_smc_vel_range = spectrum_tab['velocity'] < max_velocity
        non_smc_spectrum = spectrum_tab[non_smc_vel_range]
        non_smc_min_idx = np.argmin(non_smc_spectrum['opacity'])
        signal = (1-non_smc_spectrum['opacity'][non_smc_min_idx])
        sigma = signal/non_smc_spectrum['sigma_opacity'][non_smc_min_idx]
        velocity = non_smc_spectrum['velocity'][non_smc_min_idx]*(u.m/u.s).to(u.km/u.s)
        if sigma > threshold:
            mw_comp_names.append(name)
    return mw_comp_names




def output_mw_spectra(sbid, table, parent_folder, title, filename, threshold=None, max_velocity=70*u.km/u.s):
    print (title, filename)
    with open(filename, 'w') as f:
        output_header(f, title)

        output_location_plots(f)

        for rating in 'ABCDEF':
            targets = table[table['rating']==rating]
            if threshold:
                targets = targets[(1-targets['min_opacity'])/targets['sd_cont'] > threshold]
            comp_names = sorted(targets['comp_name'])
            mw_comp_names = find_mw_vel_comp_names(comp_names, parent_folder, threshold, max_velocity)
            print('Rating {} has {} spectra'.format(rating, len(mw_comp_names)))

            output_block_title(f, rating, rating=='A', len(mw_comp_names))

            for name in mw_comp_names:
                output_non_zoom_img(f, name, rating)
                    
        output_footer(f)


def find_j19_matches(gaskap_table, no_match_cat=None):
    print ('\nCross-matching with Jameson et al 2019', no_match_cat)
    j19_table = ascii.read('jameson2019.csv', format='csv')
    col_index = Column(name='index', data=1+np.arange(len(j19_table)))
    j19_table.add_column(col_index)

    coo_j19 = SkyCoord(j19_table['ra']*u.deg, j19_table['dec']*u.deg)
    coo_gaskap = SkyCoord(gaskap_table['ra']*u.deg, gaskap_table['dec']*u.deg)

    idx_j19, d2d_j19, d3d_j19 = coo_gaskap.match_to_catalog_sky(coo_j19)
    matched = d2d_j19 <= 18.5*u.arcsec # This cutoff allows for the widest separation without adding duplicates

    matched_j19_idx = idx_j19[matched]
    un_matched_j19_idx = [i for i in np.arange(len(j19_table)) if i not in matched_j19_idx]
    j19_unmatched = j19_table[un_matched_j19_idx]
    print ("Found {} sources in Jameson et al 2019 not in GASKAP data.".format(len(j19_unmatched)))
    coo_j19_unm = SkyCoord(j19_unmatched['ra']*u.deg, j19_unmatched['dec']*u.deg)
    idx_gaskap, d2d_gaskap, d3d_gaskap = coo_j19_unm.match_to_catalog_sky(coo_gaskap)
    close_gaskap_comp_names = gaskap_table[idx_gaskap]['comp_name']
    col_closest = Column(name='closest_gaskap', data=close_gaskap_comp_names)
    col_gaskap_ra = Column(name='gaskap_ra', data=gaskap_table[idx_gaskap]['ra'])
    col_gaskap_dec = Column(name='gaskap_dec', data=gaskap_table[idx_gaskap]['dec'])
    col_sep = Column(name='gaskap_sep', data=d2d_gaskap.to(u.arcsec))
    j19_unmatched.add_columns([col_closest, col_gaskap_ra, col_gaskap_dec, col_sep])
    if no_match_cat:
        print (j19_unmatched)
        j19_unm_vo_table = from_table(j19_unmatched)
        writeto(j19_unm_vo_table, no_match_cat)

    return j19_table, idx_j19, d2d_j19, matched, j19_unmatched


def output_j19_comparison(sbid, gaskap_table, j19_table, idx_j19, d2d_j19, j19_match, j19_unmatched, title, filename, match_cat=None): 
    print (title, filename)

    gaskap_targets = gaskap_table[j19_match]
    j19_targets = j19_table[idx_j19]
    j19_targets = j19_targets[j19_match]
    sort_order = gaskap_targets.argsort(['comp_name'])
    #comp_names = sorted(targets['comp_name'])
    gaskap_tgt_ordered = gaskap_targets[sort_order]
    j19_tgt_ordered = j19_targets[sort_order]

    with open(filename, 'w') as f:
        output_header(f, title)

        for rating in 'ABCDEF':
            mask = gaskap_tgt_ordered['rating']==rating
            subset = gaskap_tgt_ordered[mask]
            j19_subset = j19_tgt_ordered[mask]
            print('Rating {} has {} spectra'.format(rating, len(subset)))

            output_block_title(f, rating, rating=='A', len(subset))

            for idx, gaskap_src in enumerate(subset):
                gaskap_name  = gaskap_src['comp_name']
                j19_name = j19_subset[idx]['Source']
                output_j19_img(f, gaskap_name, j19_name, rating)

        # Add a section for missed spectra
        output_block_title(f, None, False, len(j19_unmatched))
        for row in j19_unmatched:
            gaskap_name  = row['closest_gaskap']
            j19_name = row['Source']
            output_j19_img(f, gaskap_name, j19_name, rating, sep=row['gaskap_sep'])

        output_footer(f)

    if match_cat:
        augmented_table = Table(gaskap_tgt_ordered)
        close_j19_comp_names = j19_tgt_ordered['Source']
        col_closest = Column(name='closest_j19', data=close_j19_comp_names)
        col_gaskap_ra = Column(name='j19_ra', data=j19_tgt_ordered['ra']*u.deg)
        col_gaskap_dec = Column(name='j19_dec', data=j19_tgt_ordered['dec']*u.deg)
        sep_vals = d2d_j19[j19_match]
        sep_vals_sorted = sep_vals[sort_order]
        col_sep = Column(name='j19_sep', data=sep_vals_sorted.to(u.arcsec))
        augmented_table.add_columns([col_closest, col_gaskap_ra, col_gaskap_dec, col_sep])
        print (augmented_table)
        j19_match_vo_table = from_table(augmented_table)
        writeto(j19_match_vo_table, match_cat)


def main():
    args = parseargs()

    start = time.time()
    print("#### Creating Started ASKAP spectra extraction for sbid {} at {} ####".format(args.sbid,
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))))

    parent_folder = 'sb{}/'.format(args.sbid)
    if args.parent:
        parent_folder = args.parent
    if not os.path.exists(parent_folder):
        print("Error: Folder {} does not exist.".format(parent_folder))
        return 1

    spectra_votable = votable.parse('{}/askap_spectra.vot'.format(parent_folder), pedantic=False)
    spectra_table = spectra_votable.get_first_table().to_table()

    output_spectra(args.sbid, spectra_table, 'Absorption spectra for SBID {}'.format(
        args.sbid), '{}/index.html'.format(parent_folder))
    output_spectra(args.sbid, spectra_table, 'Absorption spectra for SBID {} with {}σ candidate detections'.format(
        args.sbid, args.good), '{}/detections.html'.format(parent_folder), threshold=args.good)
    output_spectra(args.sbid, spectra_table, 'Absorption spectra for SBID {} with {}σ candidate detections'.format(
        args.sbid, args.best), '{}/best.html'.format(parent_folder), threshold=args.best, verbose=True)

    output_mw_spectra(args.sbid, spectra_table, parent_folder, 'Absorption spectra for SBID {} with {}σ candidate Milky Way detections'.format(
        args.sbid, args.good), '{}/mw_detections.html'.format(parent_folder), threshold=args.good)

    if args.sbid in (8906, 10941, 10944):
        j19_table, idx_j19, d2d_j19, j19_match, j19_unmatched = find_j19_matches(spectra_table, no_match_cat='{}/j19_not_matched.vot'.format(parent_folder))
        output_j19_comparison(args.sbid, spectra_table, j19_table, idx_j19, d2d_j19, j19_match, j19_unmatched,
            'Absorption spectra for SBID {} also in Jameson 19'.format(args.sbid), '{}/j19.html'.format(parent_folder), match_cat='{}/askap_spectra_in_j19.vot'.format(parent_folder))

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Extracted %d spectra in %.02f s' %
          (len(spectra_table), end - start))

    return 0



if __name__ == '__main__':
    exit(main())
