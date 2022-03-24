# Test python script for use within CASA
# e.g. casa --nologger --log2term -c test_script.py

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
    return

if __name__ == '__main__':
    main()
    exit()
