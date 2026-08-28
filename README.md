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

You can run `Rscript scripts/check_R_libraries.R` to check if you have all of the necessary libraries and install any that are missing.

Here is the output of `utils::sessionInfo()` after running the `check_R_libraries.R` script:
```
R version 4.3.3 (2024-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 24.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] paletteer_1.6.0 cowplot_1.1.3   survival_3.8-6  RcppTOML_0.2.3 
 [5] optparse_1.7.5  lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1  
 [9] dplyr_1.1.4     purrr_1.1.0     readr_2.1.5     tidyr_1.3.1    
[13] tibble_3.3.0    ggplot2_3.5.2   tidyverse_2.0.0 here_1.0.1     

loaded via a namespace (and not attached):
 [1] Matrix_1.6-5       gtable_0.3.6       rematch2_2.1.2     compiler_4.3.3    
 [5] Rcpp_1.1.0         tidyselect_1.2.1   splines_4.3.3      scales_1.4.0      
 [9] lattice_0.22-5     R6_2.6.1           generics_0.1.4     rprojroot_2.0.4   
[13] pillar_1.11.0      RColorBrewer_1.1-3 tzdb_0.5.0         rlang_1.1.6       
[17] getopt_1.20.4      stringi_1.8.7      timechange_0.3.0   cli_3.6.5         
[21] withr_3.0.2        magrittr_2.0.3     grid_4.3.3         hms_1.1.3         
[25] lifecycle_1.0.4    vctrs_0.6.5        glue_1.8.0         farver_2.1.2      
[29] tools_4.3.3        pkgconfig_2.0.3
```

### Required Julia libraries

`ArgParse`, `CSV`, `DataFramesMeta`, `Distributions`, `DrWatson`, `Random`, `Reexport`, `TOML`

The julia project can be set up by running `julia --project=. --eval "using Pkg; Pkg.instantiate()"` from the project's root directory.

Most figures can be generated without julia, using only R and the required libraries. See `scripts/README.md` for more details.