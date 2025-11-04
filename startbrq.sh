#!/bin/bash

myname="$1"
XILLVER="$2"

incls=( "5" "30" )
alphas=( "-3" "-6" )
accrates=( "0.01" "0.1" "0.2" "0.3" )
#incls=( "5" )
#alphas=( "-3" )
#accrates=( "0.01" )
if [ "$BINAC2" ]; then
    workdirbase="$WORK"/"$myname"
else
    workdirbase="$HOME"/Gartenzwerg/"$myname"
fi
if [ ! -d "$workdirbase" ]; then
    echo "can't execute, workdir not present"
    exit 1
fi
athfile="$workdirbase/athinput.pp_master_project"
export SPIN=$(grep -E "^\s*a\s*=" "$athfile" | awk -F'=' '{print $2}' | awk '{print $1}')
export XILLVER
for incl in "${incls[@]}"; do
  for alpha in "${alphas[@]}"; do
    for accrate in "${accrates[@]}"; do
      #start the job here
      export WORKDIR="$workdirbase"/br_"$incl"_"$alpha"_"$accrate"
      export INCLINATION="$incl"
      export ALPHA="$alpha"
      export DISKLOC="$workdirbase"/"dshape$accrate.csv"
      mkdir -p "$WORKDIR"
      if [ "$BINAC2" ]; then
        sbatch --partition compute -t 20:00:00 -N 1 --ntasks-per-node=1 --cpus-per-task 2 --mem-per-cpu=8gb -J br_"$myname"_"$incl"_"$alpha"_"$accrate" --output="$WORKDIR"/LOG_BR --error="$WORKDIR"/LOG_BR --export=ALL startjobbr.sh
      else
        ./startjobbr.sh 2>&1 | tee "$WORKDIR"/LOG_BR
      fi
    done
  done
done


