#!/bin/bash
set -e

echo "RUNNING SIMULATION"
julia ../fig3/run_model.jl --config config.toml

echo "COLLECTING RESULTS"
julia ../fig3/collect_results.jl --config config.toml --clean true

echo "GENERATING FIGURE"
Rscript sim_validation.R --config config.toml