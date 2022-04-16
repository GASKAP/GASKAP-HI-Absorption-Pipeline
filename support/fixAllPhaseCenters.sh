#!/bin/bash
set -e
echo "Updating phase centres of $# measurement sets"
date

# Loop over a set of measurement sets and correct their pointings for use with CASA
SCRIPTDIR=`dirname "$0"`
for file in "$@"; do
    echo $file
    python3 ${SCRIPTDIR}/fixPhaseCenters.py $file
    echo "---------"
done
date
