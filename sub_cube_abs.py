# CASA Python script to extract a subcube around a source. 
# The subcube will exclude short baselines, thus excluding large scale emission and optimised for 
# absorption studies.
# e.g. casa --nologger --log2term -c sub_cube_abs.py

# Author James Dempsey
# Date 27 Jan 2020

from __future__ import print_function

import csv 
import os
import shutil

def cleanup_prev(comp_name):
    image_name='sb8906/{}_sl'.format(comp_name)
    for suffix in ['image', 'model', 'pb', 'psf', 'residual', 'sumwt']:
        filename = '{}.{}'.format(image_name, suffix)
        if os.path.exists(filename):
            print ('Deleting', filename)
            shutil.rmtree(filename)
    for suffix in ['fits']:
        filename = '{}.{}'.format(image_name, suffix)
        if os.path.exists(filename):
            print ('Deleting', filename)
            os.remove(filename)


def get_target_params(sbid, sample_id):
    with open('targets_{}.csv'.format(sbid), 'r') as csvfile:
        tgt_reader = csv.reader(csvfile, delimiter=',',
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for row in tgt_reader:
            if (tgt_reader.line_num == 1):
                # skip header
                continue
            if int(row[0]) == sample_id:
                #print (row)
                comp_name = row[1]
                ra = float(row[2])
                dec = float(row[3])
                beams = row[4:]
                return comp_name, ra, dec, beams
    raise Exception('Unknown sample id={}'.format(sample_id))


def get_ms_pattern(sbid):
    with open('data_loc.csv', 'r') as csvfile:
        tgt_reader = csv.reader(csvfile, delimiter=',',
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for row in tgt_reader:
            if (tgt_reader.line_num == 1):
                # skip header
                continue
            if int(row[0]) == sbid:
                #print (row)
                pattern = row[1]
                return pattern
    raise Exception('Unknown sbid={}'.format(sbid))


def main():
    sample_id = int(os.environ['SAMPLE_ID'])
    sbid = int(os.environ['SBID'])

    comp_name, ra, dec, beams = get_target_params(sbid, sample_id)
    print ("Starting extract of subcube for sbid {} sample {} comp {}".format(sbid, sample_id, comp_name))
    print (comp_name, ra, dec, beams)
    vis = []
    pattern = get_ms_pattern(sbid)
    for beam in beams:
        num = beam[0:2]
        interleave = beam[2]
        vis_name = pattern.format(interleave, num)
        vis.append(vis_name)

    cleanup_prev(comp_name)

    print ('Processing input visibilities of ' + str(vis))
    klambda = '1.6'
    image_name='sb{}}/{}_sl'.format(sbid, comp_name)
    uvdist='>{}Klambda'.format(klambda)
    fits_name=image_name + '.fits'
    phasecenter='J2000 {}deg {}deg'.format(ra, dec)
    tclean (vis=vis,specmode='cube',imagename=image_name,reffreq='1.42040571183GHz',restfreq='1.42040571183GHz',
      phasecenter=phasecenter,imsize=50,uvrange=uvdist,
      gridder='standard', width='1km/s',
      vptable='../ASKAP_AIRY_BP.tab')
    if not os.path.exists(image_name+'.image'):
        print ('ERROR: tclean did not produce an image')
        return 1

    exportfits(imagename=image_name+'.image', fitsimage=fits_name,velocity=True)
    if not os.path.exists(fits_name):
        print ('ERROR: exportfits did not produce a fits file')
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
