#!/bin/bash

if [ "$BINAC2" ]; then
    module load devel/cmake/3.29
fi

cmake -B build -DCMAKE_BUILD_TYPE=Release
#mkdir -p data
#mkdir -p ironline_data
#mkdir -p output

cd build
make
