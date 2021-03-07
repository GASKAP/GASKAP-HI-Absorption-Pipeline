#!/usr/bin/bash
##### Select resources #####
#PBS -N ASKAP_Abs_Prep
#PBS -l select=1:ncpus=1
#PBS -l place=vscatter
#PBS -l walltime=1:00:00
##### Queue #####
#PBS -q smallmem
##### Mail Options #####
# Send an email at job start, end and if aborted
#PBS -m a
##### Change to your working directory #####
cd /avatar/jdempsey/abs_cutouts

date

# If not run in PBS, get the sbid as the second param
if [ -z "${SBID}" ]; then SBID=${1}; fi
export SBID

#  Show list of CPUs you ran on, if you're running under PBS
if [ -n "$PBS_NODEFILE" ]; then cat $PBS_NODEFILE; fi

##### Execute Program #####
casa -c list_beams.py
