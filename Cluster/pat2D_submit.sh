#!/bin/bash

# creates arrays ($(seq BeginningValue Step EndValue)) , includes beginning value AND includes end value
par_vec=($(seq 0.9 0.1 2.1))

# length of the created arrays
length_par=${#par_vec[@]}

for (( i=0; i<${length_par}; i++ ));
do
        JOBNAME="scan_par_${par_vec[i]}"
	JOB=`qsub -S /bin/sh -N $JOBNAME -cwd -o output.txt -e error.txt -q ehud.q << EOJ
	/storage/matlab/bin/matlab << M_PROG
	pat2D(${par_vec[i]});
	M_PROG
	EOJ
	`
	echo "JobID = ${JOB} for parameters ${par_vec[i]} ${name} submitted on `date`" 
done 

exit 
