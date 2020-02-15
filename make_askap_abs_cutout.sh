#!/usr/bin/bash

# Creates a cutout from an ASKAP measurement set. 
# The cutout will exclude large scale emission and thus be tuned for absorption usage.

SAMPLE_ID="$1"
status_folder="$2"

# Run the cutout script in CASA
export SAMPLE_ID
#casa -c sub_cube_abs.py
casa --nologger --log2term -c sub_cube_abs.py

# Mark the job as complete
retval=$?
echo "Return value=${retval}"
if [ $retval -ne 0 ]; then
    echo "Job failed"
    echo `date` > ${status_folder}/${SAMPLE_ID}.FAILED
    exit 1
fi

echo `date` > ${status_folder}/${SAMPLE_ID}.COMPLETED
exit 0
