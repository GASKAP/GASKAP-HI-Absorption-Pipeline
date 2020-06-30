# GenerateSpectraPages.py

import argparse
import os
import time

from astropy.io import votable
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

def output_block_title(f, rating, first, count):
    if not first:
        f.write('\n\n</div><br/>\n')
    spec = 'spectrum' if count == 1 else 'spectra'
    f.write('\n<div>')
    f.write('\n<div class="col-9 d-inline"><h2 class="d-inline font-weight-light text-center text-lg-left mt-4 mb-0">{} Rating {} {}</h2></div>'.format(count, rating, spec))
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

def output_j19_img(f, comp_name, rating):
    f.write('\n<div class="col-4">')
    f.write('\n<a href="spectra/{0}_spec.png" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(comp_name, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="spectra/{0}_spec_zoom.png" alt="Zoomed preview of spectrum at {0}">'.format(comp_name))
    f.write('\n</a>\n</div>')
    f.write('\n<div class="col-8">')
    f.write('\n<a href="jameson19plots/{0}.jpg" class="d-block mb-4 h-100"  data-lightbox="rating{1}">'.format(comp_name, rating))
    f.write('\n<img class="img-fluid img-thumbnail" ')
    f.write('src="jameson19plots/{0}.jpg" alt="Zoomed preview of spectrum at {0}">'.format(comp_name))
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


def output_j19_comparison(sbid, table, title, filename): 
    print (title, filename)
    crossmatch = ['J010029-713825', 'J003809-735023', 'J012639-731501', 'J012350-735041', 'J011132-730209', 
                'J005611-710707', 'J011005-722648', 'J012930-733310', 'J013243-734413', 'J011919-710522', 
                'J011917-710537', 'J012924-733152', 'J012629-732714', 'J013032-731740', 'J005820-713039',
                'J010919-725600', 'J004808-741205', 'J004956-723553', 'J005337-723143', 'J005557-722604', 
                'J010931-713454', 'J011056-731403', 'J011432-732142', 'J013147-734941', 'J003947-713734',
                'J005238-731244', 'J003824-742211', 'J003939-714141', 'J011629-731438', 'J003754-725156',
                'J003801-725210', 'J004047-714600', 'J005636-740315', 'J005652-712300', 'J005732-741243',
                'J011049-731425', 'J011049-731428']

    mask = np.isin(table['comp_name'], crossmatch)
    targets = table[mask]
    comp_names = sorted(targets['comp_name'])

    with open(filename, 'w') as f:
        output_header(f, title)

        for rating in 'ABCDEF':
            targets = table[table['rating']==rating]
            mask = np.isin(targets['comp_name'], crossmatch)
            targets = targets[mask]
            comp_names = sorted(targets['comp_name'])
            print('Rating {} has {} spectra'.format(rating, len(comp_names)))

            output_block_title(f, rating, rating=='A', len(comp_names))

            for name in comp_names:
                output_j19_img(f, name, rating)

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

    output_spectra(args.sbid, spectra_table, 'Absorption spectra for SBID {}'.format(
        args.sbid), '{}/index.html'.format(parent_folder))
    output_spectra(args.sbid, spectra_table, 'Absorption spectra for SBID {} with {}σ candidate detections'.format(
        args.sbid, args.good), '{}/detections.html'.format(parent_folder), threshold=args.good)
    output_spectra(args.sbid, spectra_table, 'Absorption spectra for SBID {} with {}σ candidate detections'.format(
        args.sbid, args.best), '{}/best.html'.format(parent_folder), threshold=args.best)

    if args.sbid == 8906:
        output_j19_comparison(args.sbid, spectra_table, 'Absorption spectra for SBID {} also in Jameson 19'.format(
            args.sbid), '{}/j19.html'.format(parent_folder))

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Extracted %d spectra in %.02f s' %
          (len(spectra_table), end - start))

    return 0



if __name__ == '__main__':
    exit(main())
