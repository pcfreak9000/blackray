#!/bin/bash
#mkdir -p "$WORKDIR"
if [ "$BINAC2" ]; then
    source "$HOME"/miniconda3/etc/profile.d/conda.sh
    conda activate master_project_env
fi
./run.py "$DISKLOC" $SPIN $INCLINATION $ALPHA "$XILLVER" "$WORKDIR"
