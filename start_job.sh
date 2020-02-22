#!/usr/bin/bash
##### Select resources #####
#PBS -N ASKAP_Abs_Cutout
#PBS -l select=1:ncpus=3:mem=40g
#PBS -l place=vscatter
#PBS -l walltime=10:00:00
##### Queue #####
#PBS -q smallmem
##### Mail Options #####
# Send an email at job start, end and if aborted
#PBS -m a
##### Change to your working directory #####
cd /avatar/jdempsey/smc/abs_cutouts

date

# If not run in PBS, get the array id as the first param
if [ -z "${COMP_INDEX}" ]; then COMP_INDEX=${1}; fi

# If not run in PBS, get the sbid as the second param
if [ -z "${SBID}" ]; then SBID=${2}; fi

#  Show list of CPUs you ran on, if you're running under PBS
if [ -n "$PBS_NODEFILE" ]; then cat $PBS_NODEFILE; fi

##### Execute Program #####
bash ./make_askap_abs_cutout.sh ${COMP_INDEX} ${SBID} 'status/${SBID}'
