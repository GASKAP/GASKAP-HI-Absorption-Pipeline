"""
12/4/19
This script (slightly edited version of Emil Lenc's) reads in a measurement set (MS) file output directly from the ASKAPsoft pipeline
to apply the following alterations to various sub-tables:
1 - save original phase center as column as a binary .field file (in case we need to use it again)
2 - uses BEAM_OFFSET values for associated beam to apply necessary shift to PHASE_DIR, DELAY_DIR, and REFERENCE_DIR values in FIELD table
3 - deletes rows from POINTING table 

Note that this must be run in the ACES virtual environment in an interactive CASA session (as it calls casacore tasks). When the altered tables need to be applied, 
run an associated script, changeMetaData.py with the flag -casa, in a CASA session to load the altered subtables. If the original subtables need to be used, run the
same script with the -askapsoft flag.  

Usage (in CASA interactive session): 
ipython fixPhaseCenters.py /path/to/input/measurement/set
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "Nickolas.Pingel@anu.edu.au"
__status__ = "Production"
"""

##imports

#!/usr/bin/env python

import os
import sys
import numpy as np
from casacore.tables import *
from askap.footprint import Skypos

if len(sys.argv) != 2:
    sys.exit("Usage %s <ms_file>" %(sys.argv[0]))

## load MS
ms = sys.argv[1]

# Check that the observation wasn't in pol_fixed mode
ta = table("%s/ANTENNA" %(ms), readonly=True, ack=False)
ant_mount = ta.getcol("MOUNT", 0, 1)
if ant_mount[0] != 'equatorial':
    sys.exit("%s doesn't support pol_fixed mode" %(sys.argv[0]))
ta.close()

# Work out which beam is in this MS
t = table(ms, readonly=True, ack=False)
vis_feed = t.getcol('FEED1', 0, 1)
beam = vis_feed[0]

# Dump the FIELD table into a field file if it doesn't already exist.
# This ensures that there is a backup of the table in case it needs to be restored and
# also to ensure that if fix_dir.py is run again that it doesn't apply the wrong
# beam position (since the FIELD table is over-written the first time this is run).
field_file = ms.replace(".ms", ".field")
if os.path.exists(field_file) == False:
    tp = table("%s/FIELD" %(ms), readonly=True, ack=False)
    p_phase = tp.getcol('PHASE_DIR')
    p_phase.dump(ms.replace(".ms", ".field"))
    tp.close()
    print("Dumped field table for beam %d" %(beam))

## open POINTING table to delete row
tpoint = table("%s/POINTING" % (ms), readonly = False, ack=False)
a = tpoint.rownumbers()
tpoint.removerows(a)
print('Deleted POINTING table from: %s' % (ms))

# Load the FIELD table from the field file.
ms_phase = np.load(field_file, allow_pickle=True)
# Work out how many fields are in the MS.
n_fields = ms_phase.shape[0]

# Open up the MS FIELD table so it can be updated.
tp = table("%s/FIELD" %(ms), readonly=False, ack=False)
# Open up the MS FEED table so we can work out what the offset is for the beam.
tf = table("%s/FEED" %(ms), readonly=True, ack=False)

# The offsets are assumed to be the same for all antennas so get a list of all
# the offsets for one antenna and for the current beam. This should return offsets
# required for each field.
t1 = taql("select from $tf where ANTENNA_ID==0 and FEED_ID==$beam")
print("%d fields found" %(n_fields))
n_fields2 = t1.getcol("BEAM_OFFSET").shape[0]
# Update the beam position for each field
for field in range(n_fields):
    # Obtain the offset for the current field.
    offset = t1.getcol("BEAM_OFFSET")[field]
    # Get the pointing direction for the field
    p_phase = ms_phase[field]
    
    # Shift the pointing centre by the beam offset
    phase = Skypos(p_phase[0][0], p_phase[0][1], 9, 9)
    new_pos = phase.shift(-offset[0][0], offset[0][1])
    new_pos.rn = 15
    new_pos.dn = 15
    new_pos_str = "%s" %(new_pos)
    print("Setting position of beam %d, field %d to %s" %(beam, field, new_pos_str))
    # Update the FIELD table with the beam position
    new_ra = new_pos.ra
    if new_ra > np.pi:
        new_ra -= 2.0 * np.pi
    ms_phase[field][0][0] = new_ra
    ms_phase[field][0][1] = new_pos.dec

# Write the updated beam positions in to the MS.
tp.putcol("DELAY_DIR", ms_phase)
tp.putcol("PHASE_DIR", ms_phase)
tp.putcol("REFERENCE_DIR", ms_phase)
tp.close()
tpoint.close()
t.close()


