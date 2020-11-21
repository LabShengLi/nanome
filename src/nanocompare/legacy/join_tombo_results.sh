#!/bin/bash
#PBS -q batch # batch - up to 72h, long - up to 504h. 
#PBS -l nodes=1:ppn=1   # if I want 32 cores, using machines that has 16 cores, I should specify nodes=2:ppn=16 and be sure 
##PBS -l mem=1000mb # this will ask for 1000mb of RAM; 300gb
#PBS -l walltime=72:00:00
#PBS -M wojciech.rosikiewicz@jax.org
##PBS -m abe # be emailed at the beginning, end and error
#PBS -o /home/rosikw/qsubReports
#PBS -e /home/rosikw/qsubReports
##PBS -t 1-24 # array 
cd $PBS_O_WORKDIR # change the directory for the one, from which the job was submitted

tomboNA19240BaseDir=/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/automated_Tombo_runs/NA19240_perChr

for f in $(ls NA19240_perChr.chr*.bed)
do
  cat $f >> NA19240_allChr.bed
  echo "$f finished"
done

