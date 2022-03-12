import csv
import glob
import os

def get_ms_pattern(sbid):
    with open('data_loc.csv', 'r') as csvfile:
        tgt_reader = csv.reader(csvfile, delimiter=',',
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for row in tgt_reader:
            if (tgt_reader.line_num == 1):
                # skip header
                continue
            if str(row[0]) == sbid:
                #print (row)
                pattern = row[1]
                return pattern
    raise Exception('Unknown sbid={}'.format(sbid))


sbid = str(os.environ['SBID'])
ms_template = get_ms_pattern(sbid)
target_path = ms_template.replace('{0}', '*').replace('{1}', '*')
print("working on ", target_path)
file_list = glob.glob(target_path)
file_list.sort()
#print("files for:", target_path, file_list)
beam_listing = 'sb{0}/beam_listing_SB{0}.csv'.format(sbid)

with open(beam_listing, 'w') as outfile:
    for beam_file_name in file_list:
		print (beam_file_name)
                res = vishead( vis=beam_file_name, mode='list',listitems=['field', 'ptcs'])
                for i in range(0,3):
                        coords = res['ptcs'][0]['r'+str(i+1)]
                        outfile.write('"{}","{}",{},{}\n'.format(beam_file_name,res['field'][0][i], coords[0][0][0], coords[1][0][0]))
