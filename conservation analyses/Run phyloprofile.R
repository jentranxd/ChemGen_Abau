if (!requireNamespace("devtools"))
  install.packages("devtools")
devtools::install_github("BIONF/PhyloProfile", INSTALL_opts = c('--no-lock'), build_vignettes = FALSE)


library(dplyr)
library(svglite)
library(PhyloProfile)

runPhyloProfile()

