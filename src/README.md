Holds code used by scripts to generate output.

---

## Code for analytical VE analyses

### `full_VE_model.R`

Includes implementation of mathematical models and code to generate expected counts of cases or controls for use in regression models.

### `data_generation.R`

Includes some helper functions used to generate parameter combinations and perform common calculations.

## Code for individual-based model

### `PreVaxImmunityVsMeasuredVE.jl`

Defines julia module for the simulation module. This helps to load necessary libraries, pre-compile model code, and export model functions for use in the simulation script.

### `abm.jl`

Core codebase for the simulation model.