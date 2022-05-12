setwd(paste0(base_dir, 'R'))
# Installation or load of libraries
packages <- c("devtools",                      # R-Packages
              "readxl", "writexl", "openxlsx", # Excel files
              "stats", "tseries",              # Timeseries
              "tidyverse", 'latex2exp', 'xts', # Escentials
              "psych", 'xtable',               # Summary mamado.
              "MASS",                          # Box-Cox
              "pracma",                        # Detrend  
              "forecast",                      # Check residuals  
              "stargazer",                     # Latex tables.  
              "simts",                         # SARMA.  
              "strucchange"                    # Structural change  
              )
# Install and load packages.              
if (!all(packages %in% rownames(installed.packages()))){
  idx <- which((packages %in% rownames(installed.packages()))==F)
  install.packages(packages[idx])
}

lapply(packages, library, character.only = T)
# Remove unused variables.
rm(packages)
