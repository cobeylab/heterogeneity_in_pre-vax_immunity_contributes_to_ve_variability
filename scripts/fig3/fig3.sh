#!/bin/bash

config=$1

echo "RUNNING SIMULATION"
julia run_model.jl --config $1

echo "COLLECTING RESULTS"
julia collect_results.jl --config $1

echo "GENERATING FIGURE"
Rscript fig3.R --config $1