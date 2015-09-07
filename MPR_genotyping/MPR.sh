#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=100:00:00
#PBS -d ./

export R_LIBS="/bigdata/stajichlab/cjinfeng/software/R_package_MPR/"
#cd $PBS_O_WORKDIR

cat MPR_run.R | R --slave 

echo "Done"
