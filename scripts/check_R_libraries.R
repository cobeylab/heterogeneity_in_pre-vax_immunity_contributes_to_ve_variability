if (!require(here)) {
    install.packages("here")
    library(here)
}

if (!require(tidyverse)) {
    install.packages("tidyverse")
    library(tidyverse)
}

if (!require(optparse)) {
    install.packages("optparse")
    library(optparse)
}

if (!require(RcppTOML)) {
    install.packages("RcppTOML")
    library(RcppTOML)
}

if (!require(survival)) {
    install.packages("survival")
    library(survival)
}

if (!require(cowplot)) {
    install.packages("cowplot")
    library(cowplot)
}

if (!require(paletteer)) {
    install.packages("paletteer")
    library(paletteer)
}

cat("\n---------- SESSION INFO ----------\n")
print(utils::sessionInfo(), locale = FALSE)