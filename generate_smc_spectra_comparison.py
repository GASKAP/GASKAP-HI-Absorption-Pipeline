# GenerateSpectraPages.py

import argparse
import os
import time

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii, votable
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

def output_block_title(f, rating, first, count, total=None):
    if not first:
        f.write('\n\n</div><br/>\n')
    spec = 'spectrum' if count == 1 else 'spectra'
    f.write('\n<div>')
    f.write('\n<div class="col-9 d-inline"><h2 class="d-inline font-weight-light text-center text-lg-left mt-4 mb-0">')
    if total:
        f.write('{} (of {}) Rating {} {}</h2></div>'.format(count, total, rating, spec))
    else:
        f.write('{} Rating {} {}</h2></div>'.format(count, rating, spec))
    f.write('\n<div class="col-3 pull-right d-inline"><a class="btn btn-primary" data-toggle="collapse" href="#spectra{0}" role="button" aria-expanded="false" aria-controls="spectra{0}" style="font-size: x-small;">Toggle</a></div>'.format(rating))
    f.write('\n</div>')
    f.write('\n<div class="row text-center text-lg-left collapse show" id="spectra{}">'.format(rating))

    
def output_img(f, comp_name, rating):
    f.write('\n<div class="col-lg-3 col-md-4 col-6 px-2">')
    f.write('\n<a href="spectra/{0}_spec.png" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(comp_name, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="spectra/{0}_spec_zoom.png" alt="Zoomed preview of spectrum at {0}">'.format(comp_name))
    f.write('\n</a>\n</div>')
    return

def output_footer(f):
    f.write('\n\n</div>\n</div>\n</body>\n</html>')
    return

def output_comp_img(f, comp_name, rating, othersbid, other_name):
    f.write('\n<div class="col-6">')
    f.write('\n<a href="spectra/{0}_spec.png" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(comp_name, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="spectra/{0}_spec_zoom.png" alt="Zoomed preview of spectrum at {0}">'.format(comp_name))
    f.write('\n</a>\n</div>')
    f.write('\n<div class="col-6">')
    f.write('\n<a href="../sb{2}/spectra/{0}_spec.png" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(other_name, rating, othersbid))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="../sb{1}/spectra/{0}_spec_zoom.png" alt="Zoomed preview of spectrum at {0}">'.format(other_name, othersbid))
    f.write('\n</a>\n</div>')

    return

def output_j19_img(f, comp_name, j19_comp_name, rating):
    f.write('\n<div class="col-6">')
    f.write('\n<a href="spectra/{0}_spec.png" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(comp_name, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="spectra/{0}_spec_zoom.png" alt="Zoomed preview of spectrum at {0}">'.format(comp_name))
    f.write('\n</a>\n</div>')
    f.write('\n<div class="col-6">')
    f.write('\n<a href="../jameson2019/spectra/{0}_zoom.png" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(j19_comp_name, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="../jameson2019/spectra/{0}_zoom.png" alt="Zoomed preview of spectrum at {0}">'.format(j19_comp_name))
    f.write('\n</a>\n</div>')

    return

def output_d20_img(f, comp_name, d20_comp_name, rating):
    f.write('\n<div class="col-6">')
    f.write('\n<a href="spectra/{0}_spec.png" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(comp_name, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="spectra/{0}_spec_zoom.png" alt="Zoomed preview of spectrum at {0}">'.format(comp_name))
    f.write('\n</a>\n</div>')
    f.write('\n<div class="col-6">')
    f.write('\n<a href="../dempsey2020/spectra/{0}_plot.png" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(d20_comp_name, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="../dempsey2020/spectra/{0}_plot_zoom.png" alt="Zoomed preview of spectrum at {0}">'.format(d20_comp_name))
    f.write('\n</a>\n</div>')

    return


def output_spectra(sbid, table, title, filename, threshold=None):
    print (title, filename)
    with open(filename, 'w') as f:
        output_header(f, title)

        for rating in 'ABCDEF':
            targets = table[table['rating']==rating]
            if threshold:
                targets = targets[(1-targets['min_opacity'])/targets['sd_cont'] > threshold]
            comp_names = sorted(targets['comp_name'])
            print('Rating {} has {} spectra'.format(rating, len(comp_names)))

            output_block_title(f, rating, rating=='A', len(comp_names))

            for name in comp_names:
                output_img(f, name, rating)
                    
        output_footer(f)


def crossmatch_8906(spectra, sep_limit=4*u.arcsec):
    vot_8906 = votable.parse_single_table('sb8906/askap_spectra.vot', pedantic=False)
    table_8906 = vot_8906.to_table()
    sky_8906 = SkyCoord(ra=table_8906['ra']*u.deg, dec=table_8906['dec']*u.deg)
    sky_tgt = SkyCoord(ra=spectra['ra']*u.deg, dec=spectra['dec']*u.deg)
    
    idx_8906, d2d_8906, d3d_8906 = sky_tgt.match_to_catalog_sky(sky_8906)
    match = d2d_8906 < sep_limit
    not_match = np.logical_not(match)
    
    names_8906 = np.asarray(table_8906['comp_name'])
    names_match_8906 = names_8906[idx_8906]
    names_match_8906[not_match] = ''
    return names_match_8906


def output_smc_comparison(sbid, othersbid, table, title, filename, threshold=None):
    print ('{} ({})'.format(title, filename))
    match_map = None
    if othersbid == 8906:
        crossmatch = crossmatch_8906(table)
        match_map = {k: v for k, v in zip(table['comp_name'], crossmatch)}


    with open(filename, 'w') as f:
        output_header(f, title)

        for rating in 'ABCDEF':
            targets = table[table['rating']==rating]
            if threshold:
                targets = targets[(1-targets['min_opacity'])/targets['sd_cont'] > threshold]
            comp_names = sorted(targets['comp_name'])

            if match_map:
                match_comp_names = []
                for name in comp_names:
                    sb8906_name = match_map[name]
                    if len(sb8906_name) > 0:
                        match_comp_names.append(name)
                print('Rating {} has {} matched spectra (from {})'.format(rating, len(match_comp_names), len(comp_names)))
                output_block_title(f, rating, rating=='A',  len(match_comp_names), total=len(comp_names))
            else:
                print('Rating {} has {} spectra'.format(rating, len(comp_names)))
                output_block_title(f, rating, rating=='A', len(comp_names))

            f.write('\n<div class="col-6">SB {} (Combined)</div>'.format(sbid))
            f.write('\n<div class="col-6">SB {}</div>'.format(othersbid))


            for name in comp_names:
                other_name = name
                if match_map:
                    other_name = match_map[name]
                if other_name:
                    output_comp_img(f, name, rating, othersbid, other_name)
                    
        output_footer(f)

def crossmatch_j19(spectra, sep_limit=4*u.arcsec):
    table_j19 = ascii.read("jameson2019/table1.csv")
    sky_j19 = SkyCoord(ra=table_j19['RA']*u.deg, dec=table_j19['Dec']*u.deg)
    sky_tgt = SkyCoord(ra=spectra['ra']*u.deg, dec=spectra['dec']*u.deg)
    
    idx_j19, d2d_j19, d3d_j19 = sky_tgt.match_to_catalog_sky(sky_j19)
    match = d2d_j19 < sep_limit
    not_match = np.logical_not(match)
    
    # 0029-7228_src3-0_plot.png
    names_j19 = np.asarray([row['Source'] for row in table_j19])
    names_match_j19 = names_j19[idx_j19]
    names_match_j19[not_match] = ''
    return names_match_j19


def output_j19_comparison(sbid, table, title, filename): 
    print ('{} ({})'.format(title, filename))
    crossmatch = crossmatch_j19(table)
    match_map = {k: v for k, v in zip(table['comp_name'], crossmatch)}

    with open(filename, 'w') as f:
        output_header(f, title)

        for rating in 'ABCDEF':
            targets = table[table['rating']==rating]
            #mask = np.isin(targets['comp_name'], crossmatch)
            #targets = targets[mask]
            comp_names = sorted(targets['comp_name'])
            match_comp_names = []
            for name in comp_names:
                j19_name = match_map[name]
                if len(j19_name) > 0 and os.path.exists('jameson2019/spectra/{0}_zoom.png'.format(j19_name)):
                    match_comp_names.append(name)

            print('Rating {} has {} matched spectra (from {})'.format(rating, len(match_comp_names), len(comp_names)))

            output_block_title(f, rating, rating=='A', len(match_comp_names), total=len(comp_names))

            f.write('\n<div class="col-6">SB {} (Combined) (1 km/s)</div>'.format(sbid))
            f.write('\n<div class="col-6">Jameson 2019 (0.8 km/s)</div>')

            for name in match_comp_names:
                j19_name = match_map[name]
                output_j19_img(f, name, j19_name, rating)

        output_footer(f)

def crossmatch_d20(spectra, sep_limit=4*u.arcsec):
    vot_d20 = votable.parse_single_table('dempsey2020/d20_spectra.vot', pedantic=False)
    table_d20 = vot_d20.to_table()
    sky_d20 = SkyCoord(ra=table_d20['RA']*u.deg, dec=table_d20['Dec']*u.deg)
    sky_tgt = SkyCoord(ra=spectra['ra']*u.deg, dec=spectra['dec']*u.deg)
    
    idx_d20, d2d_d20, d3d_d20 = sky_tgt.match_to_catalog_sky(sky_d20)
    match = d2d_d20 < sep_limit
    not_match = np.logical_not(match)
    
    # 0029-7228_src3-0_plot.png
    names_d20 = np.asarray(['{}_src{}'.format(row['Field'],row['Source']) for row in table_d20])
    names_match_d20 = names_d20[idx_d20]
    names_match_d20[not_match] = ''
    return names_match_d20


def output_d20_comparison(sbid, table, title, filename): 
    print ('{} ({})'.format(title, filename))
    crossmatch = crossmatch_d20(table)
    match_map = {k: v for k, v in zip(table['comp_name'], crossmatch)}

    with open(filename, 'w') as f:
        output_header(f, title)

        for rating in 'ABCDEF':
            targets = table[table['rating']==rating]
            comp_names = sorted(targets['comp_name'])
            match_comp_names = []
            for name in comp_names:
                d20_name = match_map[name]
                if len(d20_name) > 0 and os.path.exists('dempsey2020/spectra/{0}_plot.png'.format(d20_name)):
                    match_comp_names.append(name)

            print('Rating {} has {} matched spectra (from {})'.format(rating, len(match_comp_names), len(comp_names)))

            output_block_title(f, rating, rating=='A', len(match_comp_names), total=len(comp_names))

            for name in match_comp_names:
                d20_name = match_map[name]
                output_d20_img(f, name, d20_name, rating)

        output_footer(f)

    



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

    output_smc_comparison(args.sbid, 10941, spectra_table, 'Absorption spectra for SBID {} compared with {}'.format(
        args.sbid, 10941), '{}smc.html'.format(parent_folder))

    output_smc_comparison(args.sbid, 8906, spectra_table, 'Absorption spectra for SBID {} compared with {}'.format(
        args.sbid, 8906), '{}smc-8906.html'.format(parent_folder))

    output_j19_comparison(args.sbid, spectra_table, 'Absorption spectra for SBID {} also in Jameson 19'.format(
        args.sbid), '{}j19.html'.format(parent_folder))

    output_d20_comparison(args.sbid, spectra_table, 'Absorption spectra for SBID {} also in Dempsey 20'.format(
        args.sbid), '{}d20.html'.format(parent_folder))

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Extracted %d spectra in %.02f s' %
          (len(spectra_table), end - start))

    return 0



if __name__ == '__main__':
    exit(main())
