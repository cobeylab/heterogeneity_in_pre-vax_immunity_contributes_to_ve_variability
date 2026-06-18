#!/bin/bash

echo "RUNNING SIMULATION"
julia run_model.jl --config fig3.toml

echo "COLLECTING RESULTS"
julia collect_results.jl --config fig3.toml

echo "GENERATING FIGURE"
Rscript fig3.R --config fig3.toml