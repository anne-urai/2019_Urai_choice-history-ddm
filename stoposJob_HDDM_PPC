#!/bin/bash

#PBS -o /home/aeurai/jobs
#PBS -e /home/aeurai/jobs
#PBS -lnodes=1 -lwalltime=5:00:00 # specify in hours, not days

# load necessary modules
module load stopos
source activate python27 # use anaconda

# determine how many parallel jobs we can run on this node
ncores=`sara-get-num-cores`
((ncores -= 1)) # subtract one for system processes
ncores=2
echo "ncores = $ncores"
 
# loop over the cores available
for ((i=1; i<=ncores; i++)) ; do
(

  for ((j=1; j<=1; j++)) ; do
     stopos next -p poolPPC
       if [ "$STOPOS_RC" != "OK" ]; then
        break
     fi
    echo "Running with parameters: $STOPOS_VALUE"

    # see https://userinfo.surfsara.nl/systems/lisa/software/stopos
    a=( $STOPOS_VALUE )
    d=${a[0]}
    v=${a[1]}

    # first, run the model
 	  eval "python /home/aeurai/code/serialDDM/fitHDDM.py -r 2 -d $d -v $v"

    stopos remove -p poolPPC
stopos status -p poolPPC
   done
 ) &
done
wait

echo "Job $PBS_JOBID started/ended at `date`" | mail $USER -s "Job $PBS_JOBID"
