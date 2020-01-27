##### Select resources #####
#PBS -N ASKAP_Abs_Cutout
#PBS -l ncpus=2
#PBS -l walltime=02:00:00
##### Queue #####
#PBS -q smallmem
##### Mail Options #####
# Send an email at job start, end and if aborted
#PBS -m abe
##### Change to your working directory #####
#cd /avatar/jdempsey/smc/abs/8906

# If not run in PBS, get the aray id as the first param
if [ -z "${PBS_ARRAYID}" ]; then PBS_ARRAYID=${1}; fi

#  Show list of CPUs you ran on, if you're running under PBS
if [ -n "$PBS_NODEFILE" ]; then cat $PBS_NODEFILE; fi

##### Execute Program #####
bash ./make_askap_abs_cutout.sh ${PBS_ARRAYID} 'status/8906'