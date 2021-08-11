# Load packages

library(methylGSA)
library(stringr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(glue)


# Function implementation
enrich <- function(result_dir, input_data, min = 25, max = 250, array_type = "EPIC", padj_threshold = 0.05
                   ){

df <- read.csv(input_data)
df <- df[, c("CpG", "p.value")]
df_pvls <- df$p.value
names(df_pvls) <- df$CpG

dbs = c("Reactome")
results = list()

for (db in dbs){
  res = methylglm(cpg.pval = df_pvls, minsize = min, maxsize = max, GS.type = db, array.type = array_type)
  res = res[res$padj <= padj_threshold, ]
  
  path <- glue(result_dir, "_report_", db, "_.csv")
  write.csv(res, path)
  
  results$db <- res
}

return(results)
}
