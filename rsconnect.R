options(repos=BiocManager::repositories())
getOption("repos")
library(rsconnect)

rsconnect::deployApp("RNASeqDGE.Rmd")