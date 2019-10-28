### load_packages.R

# Pair with LBBC to load required packages for the vignette and figure replication.

### Functions

#ipak function copies from stevenworthington/ipak.R 
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


list_of_packages = c("devtools", "ggplot2", "reshape2", "roxygen2", "MASS", "taxize", "ggpubr", "ineq", "data.table")

ipak(list_of_packages)