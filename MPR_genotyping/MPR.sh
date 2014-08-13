#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=100:00:00

export R_LIBS="/rhome/cjinfeng/software/tools/R_package_MPR/"
cd $PBS_O_WORKDIR

cat MPR_run.R | /rhome/cjinfeng/software/tools/R-2.15.3/bin/R --slave 

echo "Done"
