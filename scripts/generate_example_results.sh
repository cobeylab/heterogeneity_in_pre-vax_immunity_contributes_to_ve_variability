#!/bin/bash
set -e

echo "---------- Figure 1 ----------"
Rscript fig1/fig1.R

echo "---------- Figure 2 ----------"
Rscript fig2/fig2.R

echo "---------- Figure 3 (+ higher/lower true vaccine protection) ----------"
cd fig3
Rscript fig3.R --config smaller_pop_size/fig3.toml
Rscript fig3.R --config smaller_pop_size/fig3_higher_vax_direct_effect.toml
Rscript fig3.R --config smaller_pop_size/fig3_lower_vax_direct_effect.toml
cd ..

echo "---------- Supplemental figures ----------"
echo "---------- comparison_of_final_ve_estimates_as_heterogeneity_varies ----------"
Rscript supplemental_figures/comparison_of_final_ve_estimates_as_heterogeneity_varies.R

echo "---------- final_ve_compared_to_average_vax_protection ----------"
Rscript supplemental_figures/final_ve_compared_to_average_vax_protection.R

echo "---------- final_ve_sensitivity_to_pre-vax_heterogeneity ----------"
Rscript supplemental_figures/final_ve_sensitivity_to_pre-vax_heterogeneity.R

echo "---------- final_ve_sensitivity_to_vax_waning_rate ----------"
Rscript supplemental_figures/final_ve_sensitivity_to_vax_waning_rate.R

echo "---------- final_ve_sensitivity_to_waning_and_pre-vax_means ----------"
Rscript supplemental_figures/final_ve_sensitivity_to_waning_and_pre-vax_means.R

echo "---------- multiyear_sim_ve_comparison ----------"
Rscript supplemental_figures/multiyear_sim_ve_comparison.R --config fig3/smaller_pop_size/fig3.toml --example true

echo "---------- proper_ve_reference_with_waning_vax_protection ----------"
Rscript supplemental_figures/proper_ve_reference_with_waning_vax_protection.R

echo "---------- starting_ve_sensitivity_to_pre-vax_means ----------"
Rscript supplemental_figures/starting_ve_sensitivity_to_pre-vax_means.R

echo "---------- ve_estimate_comparison ----------"
Rscript supplemental_figures/ve_estimate_comparison.R

echo "---------- ve_with_waning_vax_protection_sensitivity_to_pre-vax_means ----------"
Rscript supplemental_figures/ve_with_waning_vax_protection_sensitivity_to_pre-vax_means.R