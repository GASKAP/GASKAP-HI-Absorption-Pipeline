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

def cleanup_prev(sbid,comp_name, image_name, fits_name):
    #image_name='sb{}/{}_sl'.format(sbid, comp_name)
    for suffix in ['image', 'model', 'pb', 'psf', 'residual', 'sumwt', 'mask', 'image.pbcor']:
        filename = '{}.{}'.format(image_name, suffix)
        if os.path.exists(filename):
            print ('Deleting', filename)
            shutil.rmtree(filename)
    if os.path.exists(fits_name):
        print ('Deleting', fits_name)
        os.remove(fits_name)


def get_target_params(sbid, sample_id):
    with open('sb{0}/targets_{0}.csv'.format(sbid), 'r') as csvfile:
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
                velocity = row[2]
                return pattern, velocity
    raise Exception('Unknown sbid={}'.format(sbid))


def main():
    sample_id = int(os.environ['SAMPLE_ID'])
    sbid = int(os.environ['SBID'])
    chan_width = os.environ['CHAN_WIDTH'] if 'CHAN_WIDTH' in os.environ else '' # e.g. '1km/s'
    print ("Processing sbid {} sample {}".format(sbid, sample_id))

    comp_name, ra, dec, beams = get_target_params(sbid, sample_id)
    print ("Starting extract of subcube for sbid {} sample {} comp {} width {}".format(sbid, sample_id, comp_name, chan_width))
    print (comp_name, ra, dec, beams)
    vis = []
    pattern, velocity = get_ms_pattern(sbid)
    for beam in beams:
        num = beam[0:2]
        interleave = beam[2]
        vis_name = pattern.format(interleave, num)
        vis.append(vis_name)

    image_name='sb{}/work/{}_sl'.format(sbid, comp_name)
    fits_name='sb{}/cutouts/{}_sl.fits'.format(sbid, comp_name)
    cleanup_prev(sbid,comp_name, image_name, fits_name)

    print ('Processing input visibilities of ' + str(vis))
    klambda = '1.5'
    uvdist='>{}Klambda'.format(klambda)
    start_vel=velocity+'m/s'
    phasecenter='J2000 {}deg {}deg'.format(ra, dec)
    tclean (vis=vis,specmode='cube',imagename=image_name,reffreq='1.42040571183GHz',restfreq='1.42040571183GHz',
      phasecenter=phasecenter,imsize=50,uvrange=uvdist,weighting='natural', start=start_vel,
      gridder='standard', pbcor=True, width=chan_width, niter=1000, cell='1arcsec', 
      vptable='../ASKAP_AIRY_BP.tab')
    if not os.path.exists(image_name+'.image'):
        print ('ERROR: tclean did not produce an image')
        raise Exception('tclean failed')

    exportfits(imagename=image_name+'.image', fitsimage=fits_name,velocity=True)
    if not os.path.exists(fits_name):
        print ('ERROR: exportfits did not produce a fits file')
        raise Exception('exportfits failed')

    return


if __name__ == '__main__':
    main()
    exit()
