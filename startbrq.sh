#!/bin/bash

myname="$1"
minrange="$2"
maxrange="$3"
XILLVER=$HOME/xillver-a-Ec5.fits

incls=( "5" "50" )
alphas=( "-3" "-6" )
accrates=( "0.05" "0.1" )
if [ "$BINAC2" ]; then
    workdirbase="$WORK"/"$myname"
else
    workdirbase="$HOME"/Gartenzwerg/"$myname"
fi
if [ ! -d "$workdirbase" ]; then
    echo "can't execute, workdir not present"
    exit 1
fi
if [ "$BINAC2" ]; then
    source "$HOME"/miniconda3/etc/profile.d/conda.sh
    conda activate master_project_env
    ../athena/resolveb.sh "$workdirbase" $minrange $maxrange
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
        export OMP_NUM_THREADS=64
        sbatch --partition compute -t 8:00:00 -N 1 --ntasks-per-node=$OMP_NUM_THREADS --cpus-per-task 2 --mem-per-cpu=8gb -J br_"$myname"_"$incl"_"$alpha"_"$accrate" --output="$WORKDIR"/LOG_BR --error="$WORKDIR"/LOG_BR --export=ALL startjobbr.sh
      else
        if [ "$P9000_WORKER_THREADS" ]; then
          echo "Using $P9000_WORKER_THREADS threads"
          export OMP_NUM_THREADS=$P9000_WORKER_THREADS
        else
          export OMP_NUM_THREADS=1
        fi
        ./startjobbr.sh 2>&1 | tee "$WORKDIR"/LOG_BR
      fi
    done
  done
done


