#!/bin/bash

#PBS -o /home/aeurai/jobs
#PBS -e /home/aeurai/jobs
#PBS -lnodes=1 -lwalltime=120:00:00

# load necessary modules
module load stopos
source activate python27 # use anaconda

# determine how many parallel jobs we can run on this node
ncores=`sara-get-num-cores` # 16 in total on LISA normal nodes
((ncores -= 1)) # subtract one for system processes, will have 15
echo "ncores = $ncores"

# loop over the cores available
for ((i=1; i<=ncores; i++)) ; do
(

  for ((j=1; j<=1; j++)) ; do
     stopos next -p pool
       if [ "$STOPOS_RC" != "OK" ]; then
        break
     fi
    echo "Running with parameters: $STOPOS_VALUE"

    # see https://userinfo.surfsara.nl/systems/lisa/software/stopos
    a=( $STOPOS_VALUE )
    d=${a[0]}
    v=${a[1]}
    i=${a[2]}
    s=${a[3]}

    # first, run the model
 	  eval "python /home/aeurai/code/RT_RDK/b1_HDDM_run.py -r 1 -d $d -v $v -i $i -s $s"

    stopos remove -p pool
stopos status -p pool
   done
 ) &
done
wait

# echo "Job $PBS_JOBID finished at `date`" | mail $USER -s "Job $PBS_JOBID"

