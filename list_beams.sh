#!/usr/bin/bash
# Produce a list of the beam measurement sets for an sbid, along with their beam centres
# Measurement sets are located based on the location listed for the sbid in data_loc.csv in the working  directory

if [ $# -lt 1 ]; then
  echo "Usage: ${0} sbid"
  exit 1
fi

SBID=${1}
export SBID
SCRIPTDIR=`dirname "$0"`

##### Execute Program #####
casa -c ${SCRIPTDIR}/list_beams.py
