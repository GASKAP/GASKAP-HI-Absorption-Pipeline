# GenerateSpectraPages.py

import argparse
import os
import time

from astropy.coordinates import SkyCoord
from astropy.io import ascii, votable
from astropy.io.votable import parse_single_table, from_table, writeto
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
                        type=int, required=False)
    parser.add_argument("-f", "--field", help="The id of the GASKAP field",
                        type=str, required=False)
    parser.add_argument("-g", "--good", help="The sigma threshold for spectra to be included in the detections.html page",
                        type=float, default=3.0)
    parser.add_argument("-b", "--best", help="The sigma threshold for spectra to be included in the best.html page",
                        type=float, default=5.0)
    parser.add_argument("-p", "--parent", help="The parent folder for the processing, will default to sbnnn/ where nnn is the sbid.",
                        required=False)
    parser.add_argument("-c", "--catalog", help="The specytra catalgo to be used, requred for fields.",
                        required=False)
    args = parser.parse_args()

    if args.sbid and args.field:
        print("Error: Either an sbid or a field must be provided, but not both.")
        parser.print_usage()
        return 1
    if not args.sbid and not args.field:
        print("Error: Either an sbid or a field must be provided.")
        parser.print_usage()
        return 1
    if args.field and not args.catalog:
        print("Error: The catalog must be specified when a field is provided.")
        parser.print_usage()
        return 1

    return args

def output_header(f, title):
    f.write('<!doctype html>\n<html lang="en">\n<head>\n<title>{}</title>'.format(title))
    with open('style.html') as style:
        f.write(style.read())
    f.write('\n</head>\n<body>')
    f.write('\n<div class="container-fluid">')
    f.write('\n<h1 align="middle">{}</h1>'.format(title))
    return

def output_location_plots(f, source_map='figures/source_loc.png'):
    f.write('\n<div class="row px-3" id="maps">')
    f.write('\n<div class="col-md-auto"><h2 class="d-inline font-weight-light text-center text-lg-left mt-4 mb-0">Location</hs></div>')
    f.write('\n<div class="col-md-auto">')
    f.write('\nField Location')
    f.write('\n<a href="figures/field_loc.png" class="d-block mb-4 h-100"  data-lightbox="maps">')
    f.write('\n<img class="img-fluid img-thumbnail" style="height: 180px" src="figures/field_loc.png" alt="Map of the location of the field.">')
    f.write('\n</a>\n</div>')
    f.write('\n<div class="col-md-auto">')
    has_mw_loc_plot = os.path.exists(os.path.dirname(f.name)+'/figures/source_loc_mw.png')
    f.write('\n{}Absorption Locations'.format('Magellanic ' if has_mw_loc_plot else ''))
    f.write('\n<a href="{}" class="d-block mb-4 h-100"  data-lightbox="maps">'.format(source_map))
    f.write('\n<img class="img-fluid img-thumbnail" style="height: 180px" src="{}" alt="Map of the location of the sources.">'.format(source_map))
    f.write('\n</a>\n</div>')
    if has_mw_loc_plot:
        f.write('\n<div class="col-md-auto">')
        f.write('\nMilky Way Absorption Locations')
        f.write('\n<a href="{}" class="d-block mb-4 h-100"  data-lightbox="maps">'.format('figures/source_loc_mw.png'))
        f.write('\n<img class="img-fluid img-thumbnail" style="height: 180px" src="{}" alt="Map of the location of the Milky Way sources.">'.format('figures/source_loc_mw.png'))
        f.write('\n</a>\n</div>')
    print(os.path.dirname(f.name)+'/figures/long_vel.png')
    if os.path.exists(os.path.dirname(f.name)+'/figures/long_vel.png'):
        f.write('\n<div class="col-md-auto">')
        f.write('\n<a href="figures/long_vel.png" class="d-block mb-4 h-100"  data-lightbox="maps">')
        f.write('\n<img class="img-fluid img-thumbnail" style="height: 180px" src="figures/long_vel.png" alt="Longitude-velocity plot of the spectra.">')
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

    
def output_img(f, comp_name, rating, id, comment, combined=False):
    zoom_file_pattern = 'figures/{0}_combined.png' if combined else 'figures/{0}_spec_zoom.png'
    zoom_filename = zoom_file_pattern.format(comp_name)
    file_pattern = 'figures/{0}_combined.png' if combined else 'figures/{0}_spec.png'
    filename = file_pattern.format(comp_name)
    f.write('\n<div class="col-lg-3 col-md-4 col-6 px-2">')
    f.write('<figure class="figure d-block">')
    f.write('\n<a href="{0}" class="mb-4"  data-lightbox="rating{1}">'.format(filename, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="{0}" alt="Zoomed preview of spectrum at {1}">'.format(zoom_filename, comp_name))
    f.write('\n</a>')
    f.write('<figcaption class="figure-caption text-right">Source #{} {} {}</figcaption>'.format(id, comp_name, comment))
    f.write('\n</figure></div>')
    return

    
def output_non_zoom_img(f, comp_name, rating, id):
    file_pattern = 'figures/{0}_spec.png'
    filename = file_pattern.format(comp_name)
    f.write('\n<div class="col-lg-3 col-md-4 col-6 px-2">')
    f.write('<figure class="figure d-block">')
    f.write('\n<a href="{0}" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(filename, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="{0}" alt="Preview of spectrum at {0}">'.format(filename))
    f.write('\n</a>')
    f.write('<figcaption class="figure-caption text-right">Source #{} {}</figcaption>'.format(id, comp_name))
    f.write('\n</figure></div>')
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


def output_spectra(table, title, filename, threshold=None, has_other_abs=False, has_mw_abs=False, 
        verbose=False, source_map=None, max_noise=None):
    print (title, filename)
    with open(filename, 'w') as f:
        output_header(f, title)

        if source_map:
            output_location_plots(f, source_map=source_map)
        else:    
            output_location_plots(f)

        for rating in 'ABCDEF':
            targets = table[table['rating']==rating]
            if max_noise:
                targets = targets[targets['sd_cont'] < max_noise]
            if has_other_abs:
                targets = targets[targets['has_other_abs'] == 1]
            elif has_mw_abs:
                targets = targets[targets['has_mw_abs'] == 1]
            elif threshold:
                targets = targets[(1-targets['min_opacity'])/targets['sd_cont'] > threshold]
            sort_order = targets.argsort(['comp_name'])
            sorted_targets = targets[sort_order]
            comp_names = sorted_targets['comp_name']
            ids = sorted_targets['id']
            maj_axes = sorted_targets['semi_maj_axis']*2
            min_axes = sorted_targets['semi_min_axis']*2
            fluxes_int = sorted_targets['flux_int']

            print('Rating {} has {} spectra'.format(rating, len(comp_names)))

            if verbose:
              print (comp_names)
              
            output_block_title(f, rating, rating=='A', len(comp_names))

            for idx, name in enumerate(comp_names):
                comment = '{:.0f}x{:.0f}" {:.0f} mJy'.format(maj_axes[idx], min_axes[idx], fluxes_int[idx])
                output_img(f, name, rating, ids[idx], comment, combined=True)
                    
        output_footer(f)


def output_listed_spectra(table, title, filename, comp_names_list, verbose=False, source_map=None, zoomed=True):
    print (title, filename)
    with open(filename, 'w') as f:
        output_header(f, title)

        if source_map:
            output_location_plots(f, source_map=source_map)
        else:    
            output_location_plots(f)

        for rating in 'ABCDEF':
            targets = table[table['rating']==rating]
            targets = targets[np.in1d(targets['comp_name'], comp_names_list)]
            sort_order = targets.argsort(['comp_name'])
            sorted_targets = targets[sort_order]
            comp_names = sorted_targets['comp_name']
            ids = sorted_targets['id']
            maj_axes = sorted_targets['semi_maj_axis']*2
            min_axes = sorted_targets['semi_min_axis']*2
            fluxes_int = sorted_targets['flux_int']
            print('Rating {} has {} spectra'.format(rating, len(comp_names)))

            if verbose:
              print (comp_names)
              
            output_block_title(f, rating, rating=='A', len(comp_names))

            for idx, name in enumerate(comp_names):
                comment = '{:.0f}x{:.0f}" {:.0f} mJy'.format(maj_axes[idx], min_axes[idx], fluxes_int[idx])
                if zoomed:
                    output_img(f, name, rating, ids[idx], comment, combined=True)
                else:
                    output_non_zoom_img(f, name, rating, ids[idx])
                    
        output_footer(f)


def output_diff_sigma_spectra(table, title, filename, verbose=False, source_map=None, zoomed=True):
    print (title, filename)
    with open(filename, 'w') as f:
        output_header(f, title)

        if source_map:
            output_location_plots(f, source_map=source_map)
        else:    
            output_location_plots(f)

        sigma_name_map = {2.8: ['J005518-714450', 'J010401-720206', 'J005116-734000', 'J010431-720726', 'J011157-734129', 'J010532-721331', 'J002620-743741'],
        2.7:['J011332-740758', 'J003037-742903', 'J013218-715348', 'J005448-725353', 'J010556-714607', 'J012924-733153', 'J003208-735038', 'J012037-703843', 'J004306-732828'],
        2.6:['J011134-711414', 'J005715-704046', 'J003936-742018', 'J002411-735717', 'J012306-695600', 'J005014-730326', 'J002222-742825', 'J010932-713453'],
        2.5:['J014924-730231', 'J012945-701803', 'J005141-725545', 'J002826-703501', 'J002034-705526'],
        3: ['J010532-721331', 'J005448-725353', 'J010556-714607', 'J005715-704046']}

        for k,v in sigma_name_map.items():
            #v = sigma_name_map[k]
            print (k, v)
            comp_names_list = v
            if k < 2.8:
                f.write('\n</div>')
            if k < 3:
                f.write('\n<h2>Spectra included at 2+ channels of {} sigma cutoff</h2>'.format(k))
            else:
                f.write('\n</div>\n<h2>Spectra included at 3+ channels of 2.5 sigma cutoff</h2>')

            first = True
            # TODO: Switch to use source lists
            for rating in 'ABCDEF':
                targets = table[table['rating']==rating]
                targets = targets[np.in1d(targets['comp_name'], comp_names_list)]
                sort_order = targets.argsort(['comp_name'])
                sorted_targets = targets[sort_order]
                comp_names = sorted_targets['comp_name']
                ids = sorted_targets['id']
                maj_axes = sorted_targets['semi_maj_axis']*2
                min_axes = sorted_targets['semi_min_axis']*2
                fluxes_int = sorted_targets['flux_int']
                print('Rating {} has {} spectra'.format(rating, len(comp_names)))
                if len(comp_names) == 0:
                    continue

                if verbose:
                    print (comp_names)
                
                output_block_title(f, rating, first, len(comp_names))
                first = False

                for idx, name in enumerate(comp_names):
                    comment = '{:.0f}x{:.0f}" {:.0f} mJy'.format(maj_axes[idx], min_axes[idx], fluxes_int[idx])
                    if zoomed:
                        output_img(f, name, rating, ids[idx], comment, combined=True)
                    else:
                        output_non_zoom_img(f, name, rating, ids[idx])
                    
        output_footer(f)


def find_j19_matches(gaskap_table, no_match_cat=None):
    print ('\nCross-matching with Jameson et al 2019', no_match_cat)
    j19_table = ascii.read('jameson2019.csv', format='csv')
    col_index = Column(name='index', data=1+np.arange(len(j19_table)))
    j19_table.add_column(col_index)

    coo_j19 = SkyCoord(j19_table['ra']*u.deg, j19_table['dec']*u.deg)
    coo_gaskap = SkyCoord(gaskap_table['ra'], gaskap_table['dec'])

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


def output_j19_comparison(gaskap_table, j19_table, idx_j19, d2d_j19, j19_match, j19_unmatched, title, filename, match_cat=None): 
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
        #print (augmented_table)
        j19_match_vo_table = from_table(augmented_table)
        writeto(j19_match_vo_table, match_cat)


def main():
    args = parseargs()

    start = time.time()
    target = 'sbid {}'.format(args.sbid) if args.sbid else 'field {}'.format(args.field)
    print("#### Started generating spectra pages for {} at {} ####".format(target,
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))))

    parent_folder = 'sb{}/'.format(args.sbid) if args.sbid else args.field
    if args.parent:
        parent_folder = args.parent
    if not os.path.exists(parent_folder):
        print("Error: Folder {} does not exist.".format(parent_folder))
        return 1

    if args.sbid:
        title_prefix = 'Absorption spectra for SBID {}'.format(args.sbid)
    else:
        title_prefix = 'Absorption spectra for field {}'.format(args.field)

    catalog = args.catalog if args.catalog else '{}/gaskap_sb{}_abs_spectra.vot'.format(parent_folder, args.sbid)
    spectra_votable = votable.parse(catalog, pedantic=False)
    spectra_table = spectra_votable.get_first_table().to_table()
    if args.field:
        spectra_table = spectra_table[spectra_table['field']==args.field]

    output_spectra(spectra_table, title_prefix, '{}/all.html'.format(parent_folder))
    output_spectra(spectra_table, '{} with non MW absorption features'.format(title_prefix, 
        args.good), '{}/detections.html'.format(parent_folder), has_other_abs=True)
    output_spectra(spectra_table, '{} with {}σ candidate detections'.format(title_prefix, 
        args.best), '{}/best.html'.format(parent_folder), threshold=args.best)

    output_spectra(spectra_table, '{} with MW absorption features'.format(title_prefix,
        args.best), '{}/mw_detections.html'.format(parent_folder), has_mw_abs=True)

    max_noise=0.03
    output_spectra(spectra_table, '{} with less than {} noise level'.format(title_prefix, 
        max_noise), '{}/quiet.html'.format(parent_folder), max_noise=max_noise)

    if args.sbid == 10944:
        output_diff_sigma_spectra(spectra_table, 'Comparison of sigma cutoffs', '{}/sigmacomp.html'.format(parent_folder))

        missed_sources = ['J005448-725353', 'J010532-721331', 'J005014-730326', 'J012924-733153', 'J005217-730157', 'J010556-714607', 'J005141-725545', 'J004306-732828', 'J010401-720206', 
            'J010359-720144', 'J010404-720145', 'J013032-731741', 'J003524-732223', 'J010919-725600', 'J013218-715348', 'J004718-723947', 'J010431-720726', 'J005116-734000', 'J003037-742903', 
            'J003037-742901', 'J012733-713639', 'J010932-713453', 'J003936-742018', 'J004808-741206', 'J002411-735717', 'J002143-741500']
        output_listed_spectra(spectra_table, 'Absorption spectra for SBID {} excluded by changed noise'.format(
            args.sbid), '{}/excluded.html'.format(parent_folder), missed_sources)

        wide_added = ['J012639-731502', 'J012639-731502', 'J005644-725200', 'J011408-732006', 'J005217-730157']
        output_listed_spectra(spectra_table, 'Absorption spectra for SBID {} added by using 3 channels with 2.3 sigma match'.format(
            args.sbid), '{}/wide.html'.format(parent_folder), wide_added)

        bad_noise = ['J003749-735128',
            'J010932-713453',
            'J013134-700042',
            'J013742-733050',
            'J014105-722748']
        output_listed_spectra(spectra_table, 'Absorption spectra for SBID {} with poor noise estimates'.format(
            args.sbid), '{}/bad_noise.html'.format(parent_folder), bad_noise)

    if args.sbid in (8906, 10941, 10944):
        j19_table, idx_j19, d2d_j19, j19_match, j19_unmatched = find_j19_matches(spectra_table, no_match_cat='{}/j19_not_matched.vot'.format(parent_folder))
        output_j19_comparison(spectra_table, j19_table, idx_j19, d2d_j19, j19_match, j19_unmatched,
            'Absorption spectra for SBID {} also in Jameson 19'.format(args.sbid), '{}/j19.html'.format(parent_folder), match_cat='{}/askap_spectra_in_j19.vot'.format(parent_folder))
        non_j19_table = gaskap_targets = spectra_table[~j19_match]
        print (len(non_j19_table))
        output_spectra(non_j19_table, 'Absorption spectra for SBID {} not in J19 with absorption features'.format(
            args.sbid), '{}/non_j19_detections.html'.format(parent_folder), has_other_abs=True, source_map='figures/source_loc_nonj19.png')
        output_spectra(non_j19_table, 'Absorption spectra for SBID {} not in J19 with {}σ candidate detections'.format(
            args.sbid, args.best), '{}/non_j19_best.html'.format(parent_folder), threshold=args.best, source_map='figures/source_loc_nonj19.png')

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Reported %d spectra in %.02f s' %
          (len(spectra_table), end - start))

    return 0



if __name__ == '__main__':
    exit(main())
