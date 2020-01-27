#!/bin/bash

# Creates a cutout from an ASKAP measurement set. 
# The cutout will exclude large scale emission and thus be tuned for absorption usage.

sample_id="$1"
status_folder="$2"
# Run the cutout script in CASA

# Mark the job as complete
echo `date` > ${status_folder}/${sample_id}.COMPLETED