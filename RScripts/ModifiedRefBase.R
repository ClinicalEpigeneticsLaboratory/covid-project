# Function to perform cell-fraction correction based on ChAMP implentation and RPC method from
# EpiDISH package.

library("EpiDISH")
library("arrow")

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

mynorm <- data.table::fread("/home/janbinkowski/Desktop/Projects/covid-project/data/interim/NEW_ALL/myNorm.csv", data.table=F)
rownames(mynorm) <- mynorm[, 1]
mynorm <- mynorm[, -1]

# Split data into chunks
c1 <- modified_refBase(mynorm[, 1:400])
c2 <- modified_refBase(mynorm[, 401:800])
c3 <- modified_refBase(mynorm[, 801:length(colnames(mynorm))])

corrected_mynorm <- cbind(c1$CorrectedBeta, c2$CorrectedBeta, c3$CorrectedBeta)
predicted_cf <- rbind(c1$CellFractionBeforeCorrection, c2$CellFractionBeforeCorrection, c3$CellFractionBeforeCorrection)
cf_after_correction <- rbind(c1$CellFractionAfterCorrection, c2$CellFractionAfterCorrection, c3$CellFractionAfterCorrection)

dim(corrected_mynorm) == dim(mynorm)
all(row.names(corrected_mynorm) == row.names(mynorm))
all(colnames(corrected_mynorm) == colnames(mynorm))

write.csv(corrected_mynorm, "../data/processed/CorrectedMyNorms/mynorm.csv")
write.csv(cf_after_correction, "../data/processed/CF/corrected_CF.csv")
write.csv(predicted_cf, "../data/processed/CF/raw_CF.csv")
