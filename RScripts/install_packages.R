#Download the library
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("ChAMP")
BiocManager::install("methylGSA")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
install.packages("stringr")
install.packages("glue")
