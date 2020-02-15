# askap-bulk-cutouts
Produce bulk cutouts from ASKAP measurement sets

Notebooks and scripts to produce bulk cutouts from an ASKAP dataset using CASA. The scheduler ensures that each measurement set is only used by on e CASA job at a time.

Scripts are:
* *ProduceSourceCutoutScripts.ipynb* - Builds the list of targets and the measurement sets they will use based on the Selavy catalogue produced for the measurement set.
* *askap_cutout_daemon.py* - Daemon process to manage the production of cutouts from an ASKAP scheduing block. It will start and monitor a series of cutout jobs ensuring that measurement sets are only used by one job at a time.
* *make_askap_abs_cutout.sh* - Bash script to run CASA
* *start_job.sh* - PBS control Bash script - run by PBS to start each cutout job.
* *sub_cube_abs.py* - CASA control script to produce a cutout 
