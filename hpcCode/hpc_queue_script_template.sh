#!/bin/sh
#PBS -N BAM_template
#PBS -q hpc
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -o HpcOutput/${PBS_JOBNAME}_${PBS_JOBID}.out
#PBS -e HpcOutput/${PBS_JOBNAME}_${PBS_JOBID}.err
#PBS -m ae
#PBS -M s113245@student.dtu.dk
#PBS -t 1-16
# set above interval of the indices for output files and number of jobs
# -l mem=24GB USE THIS IF RUN ON DRUGS DATASET

cd $PBS_O_WORKDIR

# nFile:	1 = drugs and sideeffects
# 		2 = davis
#		3 = Animal with Attributes
#		4 = synthetic


matlab -nodisplay -r "nFile=2; T=100; d=6; k=3; method='kflip'; Percent=0.1; missing_data=1; save_f=1; iter=${PBS_ARRAYID}; Run_script_template"
