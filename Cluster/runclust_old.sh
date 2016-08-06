#!/bin/bash

######### setup proper path locations #############

#matloc="storage/matlab/bin/matlab"
#matloc="/usr/local/bin/matlab"  # For Linux
#matloc="/Applications/MATLAB_R2014b.app/bin/matlab"

rdmloc="~/Dropbox/Tools/RDM"
pathcom="tmppath=genpath('$rdmloc'); addpath(tmppath); clear tmppath;"

kernelversion="$(uname)"
servername="$(hostname)"
echo $servername
if [ "$kernelversion" == "Darwin" ]; then
#echo "we're using mac?"
	matloc="/Applications/MATLAB_R2014b.app/bin/matlab" # For Mac
else
	if [[ $servername =~ "sge" ]]; then
#echo "we're on the cluster?"
		matloc="storage/matlab/bin/matlab"  # For cluster computer
	else
#echo "we're using my office computer?"
		matloc="/usr/local/bin/matlab"  # For Linux
	fi
fi

echo $matloc

# Make sure there is sufficient input
if [ "$#" -lt 5 ]; then
    if [ "$#" -lt 1 ]; then
        printf "general run format is: \n$0 base_run_name totruns from_part to_part queue comm1 comm2 ... \n  totruns is the total division of the runs to make (totruns==0 means run locally)\n  where the parts in this division to run are [from_part .. to_part]\n  queue should give the name of the queue to run (if totruns==0, queue is ignored)\n  further arguments (comm#) are either a matlab command or a .m/.mat file to run/load\n  if no comm# are supplied, a trial run is made (not running anything, just printing)\n"
	exit
    else
	# File collection mode
	printf "attempting to load files of: $1. not implemented yet.\n"
	exit
    fi
fi

######### If we reached here, enough input was given for a "normal" or "trial" run #############

# Load parameters and get rid of them (from the $@ vector)
projname=$1
shift
totruns=$1
shift
from_part=$1
shift
to_part=$1
shift 

if [ $totruns -eq 0 ]; then	# if totruns=0, run locally (not on cluster, this ignores queue)
    queue=0
    totruns=$to_part
else
    queue=$1".q"
fi
shift

if [ "$#" -lt 1 ]; then
    echo "Working in trial run mode (no actual runs are made)"
    input=""
else
    for tmparg in "$@" # Loop over arguments
    do
        input+=",'$tmparg'"	# Concat arguments
    done
    echo "Setting up base matlab file (arguments: $input)"
fi



#echo $totruns
#echo $from_part
#echo $to_part
#echo $projname

######### Prepare base file #############

# Run matlab, setup path, and create a base file to work with
$matloc -nosplash << M_PROG
$pathcom
PrepForCluster('$projname'$input);
M_PROG

######### Starting running things!  #############

if test -e "$projname"_tmp.mat""; then  # Make sure matlab created a file just above

for (( i=$from_part; i<=$to_part; i++ ));
do
	JOBNAME="rp_$i$projname"
	
	# load base file, setup which part to run, then use the runar function, and finally save data
	matcom=" $pathcom load('$projname"_tmp.mat"'); Es.RunsChoice=[$i $totruns];  runpar(Vs,Ps,Es);"
#Es.Tmax=20*$i; tmpout=run2ss(Vs,Ps,Es); save('"$projname"_part$i"of$totruns.mat"'); "
	if [ -n "$input" ]; then
	    # Actual run here
	    if [ $queue -eq $queue 2> /dev/null ]; then	# if queue==0, run locally
    		$matloc -nosplash << M_PROG 
	    	$matcom 
M_PROG
	    else			# actual queue, so run on cluster
		JOB=`qsub -S /bin/sh -N $JOBNAME -cwd -o output.txt -e error.txt -q $queue << EOJ
		$matloc -nosplash << M_PROG 
	    	$matcom 
		M_PROG
EOJ
`
	    fi
	else
	    #if [ $queue -eq 0 ]; then
	    if [ $queue -eq $queue 2> /dev/null ]; then
	    #if [ $queue -ne 0 -o $queue -eq 0 2>/dev/null ]; then 
	    	echo "Local run" 
#	    if [ $queue -ne 0 -o $queue -eq 0 2>/dev/null ] # 
#	    if (test $queue -ne 0) then	
	    else   # if queue is not 0, then ouput Job info
	    	printf "JobID = ${JOBNAME} at queue $queue, part $i, would be submitted on `date`\n" 
	    fi
    	    printf "Would run this:   $matcom \n"
	fi
	
done
######### Wrap things up  #############
rm $projname"_tmp.mat"

else # If run failed
echo "PrepForCluster did not create file $projname"_tmp.mat". Quitting!"
fi



exit 
