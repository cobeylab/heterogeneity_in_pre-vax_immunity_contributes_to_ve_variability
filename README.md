# Heterogeneity in pre-vaccination population immunity can contribute to variability in vaccine effectiveness

## Code repository for manuscript analyses and figures

---

## Repository contents

- `src`: Shared code used across figure-generating scripts
- `scripts`: Code used to generate main-text and supplementary figures
- `data`: Location where any data generated from scripts will be stored
- `plots`: Location where scripts will store generated figures

## Dependencies

### Required R libraries

`here`, `tidyverse`, `optparse`, `RcppTOML`, `survival`, `cowplot`, `paletteer`

### Required Julia libraries

`ArgParse`, `CSV`, `DataFramesMeta`, `Distributions`, `DrWatson`, `Random`, `Reexport`, `TOML`

The julia project can be set up by running `julia --project=. --eval "using Pkg; Pkg.instantiate()"` from the project's root directory.