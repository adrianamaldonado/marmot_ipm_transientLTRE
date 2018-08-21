# list of all packages
requiredPackageNames <- 
  c("mgcv", "tidyverse", "viridis", "grid", "gridExtra", "ggpubr")

# install (if we don't have it) and load up the package 
sapply(requiredPackageNames, function(name) {
  if(!require(name, quietly=TRUE, character.only=TRUE)) {
    install.packages(name)
    library(name, quietly=TRUE, character.only=TRUE)
  }
})

# clean up
rm(requiredPackageNames)
