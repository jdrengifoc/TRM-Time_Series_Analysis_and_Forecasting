# Installation or load of libraries
packages <- c("devtools",                      # R-Packages
              "readxl", "writexl", "openxlsx", # Excel files
              "stats", "tseries",              # Timeseries
              "tidyverse", 'latex2exp'         # Escentials
              )
# Install and load packages.              
if (all(packages %in% rownames(installed.packages()))){
  lapply(packages, library, character.only = T)
} else {
  idx <- which((packages %in% rownames(installed.packages()))==F)
  install.packages(packages[idx])
  lapply(packages, library, character.only = T)
}
# Remove unused variables.
rm(base_dir, packages)
