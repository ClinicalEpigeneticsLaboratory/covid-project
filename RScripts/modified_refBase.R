#Download the library
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EpiDISH")

# Load
library("EpiDISH")

# Function implementation
modified_refBase <- function(beta, method="RPC"){
  
  data(centDHSbloodDMC.m)
  
  frac.m <- epidish(beta, ref.m = centDHSbloodDMC.m, method = method)
  cellFrac <- frac.m$estF
  
  print(colMeans(cellFrac))
  message(names(which.min(colMeans(cellFrac))), "has smallest cell proportion, all other cell proportions will be corrected by linear regression method.")
  
  lm.o <- lm(t(beta) ~ cellFrac[, -1 * which.min(colMeans(cellFrac))])
  tmp.m <- t(lm.o$res) + rowMeans(beta)
  tmp.m[tmp.m <= 0] <- min(tmp.m[which(tmp.m > 0)])
  tmp.m[tmp.m >= 1] <- max(tmp.m[which(tmp.m < 1)])
  
  frac.m2 <- epidish(tmp.m, ref.m = centDHSbloodDMC.m, method = method)
  cellFrac2 <- frac.m2$estF
  
  return(list(CorrectedBeta = tmp.m, CellFractionBeforeCorrection = cellFrac,
              CellFractionAfterCorrection = cellFrac2))
}

# Load data
mynorm <- data.table::fread("../data/interim/ALL/myNorm.csv", data.table = F)
mynorm <- data.frame(mynorm, row.names = 1)

# Run
cfc <- modified_refBase(mynorm)
