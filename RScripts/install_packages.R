if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
BiocManager::install("methylGSA")
BiocManager::install("ChAMP")

install.packages("stringr")
install.packages("glue")
