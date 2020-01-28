##### Select resources #####
#PBS -N ASKAP_Abs_Cutout
#PBS -l ncpus=2
#PBS -l walltime=06:00:00
##### Queue #####
#PBS -q smallmem
##### Mail Options #####
# Send an email at job start, end and if aborted
#PBS -m abe
##### Change to your working directory #####
cd /avatar/jdempsey/smc/abs_cutouts

# If not run in PBS, get the aray id as the first param
if [ -z "${PBS_ARRAY_INDEX}" ]; then PBS_ARRAY_INDEX=${1}; fi

#  Show list of CPUs you ran on, if you're running under PBS
if [ -n "$PBS_NODEFILE" ]; then cat $PBS_NODEFILE; fi

##### Execute Program #####
which bash
bash ./make_askap_abs_cutout.sh ${PBS_ARRAY_INDEX} 'status/8906'
