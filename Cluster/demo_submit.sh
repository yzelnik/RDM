#!/bin/bash

# creates arrays ($(seq BeginningValue Step EndValue)) , includes beginning value AND includes end value
par_vec=($(seq 1.0 3.0 7.0))

# length of the created arrays
length_par=${#par_vec[@]}

for (( i=0; i<${length_par}; i++ ));
do
        JOBNAME="scan_par_${par_vec[i]}"
	JOB=`qsub -S /bin/sh -N $JOBNAME -cwd -o output.txt -e error.txt -q intel_all.q << EOJ
	/storage/matlab/bin/matlab << M_PROG
	demo(${par_vec[i]});
	M_PROG
	EOJ
	`
	echo "JobID = ${JOB} for parameters ${par_vec[i]} ${name} submitted on `date`" 
done 

exit 
