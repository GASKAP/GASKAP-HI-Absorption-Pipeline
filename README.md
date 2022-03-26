# GASKAP HI Absorption Pipeline
Pipeline to produce HI absorption spectra from calibrated ASKAP visibilities for the GASKAP project.

## Scripts

The scripts are listed in the order in which they are normally used.

1. `prep_gaskap_abs.py` - Script to prepare a GASKAP absorption cutout run for an ASKAP scheduling block.
This script will build up metadata files describing the observation and data for use by the cutout processor. 
    1. `list_beams.sh` - Control script to invoke CASA to list the beams in an observation and their pointing centres
    1. `list_beams.py` - CASA script to list the beams in an observation and their pointing centres to a csv file
1. `askap_cutout_daemon.py` - Daemon process to manage the production of cutouts for each source from an ASKAP scheduling block.
The daemon dynamically schedules jobs to produce cutouts in order to minimise I/O contention on measurement sets and 
to balance the resource usage of the process. 
It will monitor the progress of jobs, record which ones have succeeded and failed, and track which beams are currently being used for processing. When the number of active jobs allows, the daemon will select a job that doesn't use any beams already in use and will start the job, recording the beams that are now being used. A minimum number of active jobs parameter specifies when the multiple jobs will be allowed to use the same beams in order to avoid a long processing tail all waiting for central beams.
    1. `start_job.sh` - PBS control script to start a single cutout job
    1. `make_askap_abs_cutout.sh` - Control script to produce a cutout for a single source
    1. `sub_cube_abs.py` - CASA python script to extract a cutout cube/subcube around a single source
1. `extract_emission.py` - Python script that uses astropy to produce emission spectra for a set of sources
1. `extract_pb_emission.py` - Python script that uses astropy to produce primary beam emission spectra for a set of sources
1. `extract_spectra.py` - Produce absorption spectra from a set of subcubes
1. `generate_spectra_pages.py` - Create a set of web pages suitable for browsing the spectra for an observation.

## Typical process

### Inputs

* Continuum source catalogue for the observation, normally the Selavy continuum component catalogue produced by the ASKAPSoft pipeline.
* Spectral line measurements sets - the measurement sets produced by the ASKAPSoft pipeline. These will need to be adjusted to
    1. Update the phase centres of each beam measurement set to reflect the beam pointing (a difference in processing requirements between ASKAPSoft and CASA) . e.g. using `GASKAP_Imaging/pbsScripts/alterAllMS.sh` 
    2. Concatenate beam measurement sets from different observations with the same pointings. e.g. using `GASKAP_Imaging/pbsScripts/concat_All.pbs`

### Processing

A typical processing run will use the following steps

1. Run (as a job) `prep_gaskap_abs.py` supplying the sbid, the continuum catalogue and the folder containing the measurement sets (will take ~1 min).
1. Start the cutout daemon `askap_cutout_daemon.py` this may run for a 2-7 days depending on demand, number of sources, compute power and max allowed processes.
1. Produce emission spectra using `extract_emission.py` and input an emission cube and Selavy catalogue
1. Produce primary beam emission spectra using `extract_pb_emission.py` and input GASS data and the Selavy catalogue
1. Produce absorption spectra and catalogues using `extract_spectra.py` with the cutouts, emission spectra and primary beam emission spectra produced above and the input Selavy catalogue
1. Create preview web pages using `generate_spectra_pages.py`

## Dependencies

* astropy
* CASA
* https://github.com/nipingel/GASKAP_Imaging

## Reference Data

* GASS ZEA cubes