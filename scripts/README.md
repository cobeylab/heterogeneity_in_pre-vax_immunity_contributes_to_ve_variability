Holds scripts used to generate main-text and supplementary figures. All example code assumes you are in the scripts directory (`cd scripts` from the project's root directory).

---

## Main text figures

### Figure 1: VE vs. differential mean pre-vaccination infection risk

Command: `Rscript fig1/fig1.R `

### Figure 2: VE vs. differential heterogeneity in pre-vaccination infection risk

Command: `Rscript fig2/fig2.R`

### Figure 3: Dynmaic population immunity affects pre-vaccination risk and annual VE

Commands: `cd fig3`, `bash ./fig3.sh`

Note: The single simulation uses one core and about 4.5GB of memory. The simulation and figure generation takes approximately 40 minutes to finish (on AMD Ryzen 5 PRO 6650U).

## Supplemental figures

### `final_ve_compared_to_average_vax_protection.R`

Relates the distance between final VE and average true vaccine protection to whether average protection falls close to the level when maximal effects from susceptible depletion are expected in a no-waning scenario.

Command: `Rscript supplemental_figures/final_ve_compared_to_average_vax_protection.R`

### `final_ve_sensitivity_to_pre-vax_heterogeneity.R`

Heatmap showing distance between final VE and true vaccine protection when heterogeneity in pre-vaccination risk differs between vaccinated and unvaccinated populations.

Command: `Rscript supplemental_figures/final_ve_sensitivity_to_pre-vax_heterogeneity.R`

### `final_ve_sensitivity_to_vax_waning_rate.R`

Heatmap showing distance between final VE and average true vaccine protection across values of the exogeneous infection hazard and vaccine protection waning rate.

Command: `Rscript supplemental_figures/final_ve_sensitivity_to_vax_waning_rate.R`

### `final_ve_sensitivity_to_waning_and_pre-vax_means.R`

Heatmap showing different VE outcomes across scenarios involving starting true vaccine protection, vaccine protection waning rate, and mean pre-vaccination infection risk.

Command: `Rscript supplemental_figures/final_ve_sensitivity_to_waning_and_pre-vax_means.R`

### `lewnard_et_al_2018_fig_1c_recreation.R`

Recreation of Fig. 1C from Lewnard et al. 2018, which shows that susceptible depletion's effect on VE is highest at intermediate levels of true vaccine protection.

Command: `Rscript supplemental_figures/lewnard_et_al_2018_fig_1c_recreation.R`

### `proper_ve_reference_with_waning_vax_protection.R`

Figures showing how VE can be compared to average true vaccine protection when vaccine protection wanes.

Command: `Rscript supplemental_figures/proper_ve_reference_with_waning_vax_protection.R`

### `starting_ve_sensitivity_to_pre-vax_means.R`

Same figure as main text figure 1, but for different levels of true vaccine protection.

Command: `Rscript supplemental_figures/starting_ve_sensitivity_to_pre-vax_means.R`

### `ve_model_comparison.R`

Generates scatterplots showing how similar regression-based and instantaneous-incidence-based VE estimates are to our analytical VE estimate using cumulative attack rates. Sweeps across values of mean pre-vaccination risk, pre-vaccination risk heterogeneity, exogeneous infeciton hazard, true vaccine protection, and vaccine protection waning rate.

Command: `Rscript supplemental_figures/ve_model_comparison.R`

### `ve_with_waning_vax_protection_sensitivity_to_pre-vax_means.R`

Figure showing how differential mean pre-vaccination risk and waning vaccine protection can result in VE crossing zero.

Command: `Rscript supplemental_figures/ve_with_waning_vax_protection_sensitivity_to_pre-vax_means.R`