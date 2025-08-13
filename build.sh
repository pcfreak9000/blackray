#!/bin/bash

cmake -B build -DCMAKE_BUILD_TYPE=Debug
mkdir -p data
mkdir -p ironline_data
mkdir -p output

cd build
make
