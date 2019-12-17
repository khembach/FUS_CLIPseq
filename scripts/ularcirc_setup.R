# https://github.com/VCCRI/Ularcirc

library(devtools)
devtools::install_github("VCCRI/Ularcirc", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
library("Ularcirc")

## get gene annotations
all_dbs <- Compatible_Annotation_DBs() # This will return all compatible databases
mmu_dbs <- Compatible_Annotation_DBs(search_term = 'mm10') # returns mm10 compatible databases

# download the databases
BiocManager::install(c(mmu_dbs))

## start shiny app
Ularcirc()
