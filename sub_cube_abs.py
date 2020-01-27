# CASA Python script to extract a subcube around a source. 
# The subcube will exclude short baselines, thus excluding large scale emission and optimised for 
# absorption studies.
# e.g. casa --nologger --log2term -c sub_cube_abs.py

# Author James Dempsey
# Date 27 Jan 2020

import csv 
import os

def get_target_params(sample_id):
    with open('targets.csv', 'r') as csvfile:
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

def main():
    sample_id = int(os.environ['SAMPLE_ID'])
    comp_name, ra, dec, beams = get_target_params(sample_id)
    print ("Starting extract of subcube for sample {} comp {}".format(sample_id, comp_name))
    print (comp_name, ra, dec, beams)
    vis = []
    for beam in beams:
        num = beam[0:2]
        interleave = beam[2]
        vis_name='scienceData_SB8906_SMC1-0_M344-11{}.beam{}_SL.ms'.format(interleave, num)
        vis.append(vis_name)

    print 'Processing input visibilities of ' + str(vis)
    klambda = '1.6'
    image_name='sb8906/{}_sl'.format(comp_name)
    uvdist='>{}Klambda'.format(klambda)
    fits_name=image_name + '.fits'
    phasecenter='J2000 {}deg {}deg'.format(ra, dec)
    tclean (vis=vis,specmode='cube',imagename=image_name,reffreq='1.42040571183GHz',restfreq='1.42040571183GHz',
      phasecenter=phasecenter,imsize=50,uvrange=uvdist,
      gridder='standard', width='1km/s',
      vptable='ASKAP_AIRY_BP.tab')
    exportfits(imagename=image_name+'.image', fitsimage=fits_name,velocity=True)

    return


if __name__ == '__main__':
    main()
    exit()