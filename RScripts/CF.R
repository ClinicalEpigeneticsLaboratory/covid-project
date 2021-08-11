source("modified_refBase.R")

###########################################
###########################################

mynorm <- data.table::fread("../data/interim/concated_myNorm_no_CFC/concated_cov2_cov407_hb_no_cfc.csv", data.table = F)
mynorm <- data.frame(mynorm, row.names = 1)

cfc <- modified_refBase(mynorm)

# CF before correction plot
# write.csv(cfc$CellFractionBeforeCorrection, "../data/myNormsAdditional/IBD_450K/CF_IBD_450K_noCFC_Controls.csv")

boxplot(cfc$CellFractionBeforeCorrection)
dev.copy(jpeg,filename="../data/myNormsAdditional/RA_450K/CF_RA_450K_noCFC.csv");
dev.off ();

# CF After correction plot
# write.csv(cfc$CellFractionAfterCorrection, "../data/processed/CF/CF_COV2_COV407_HB_CFC.csv")

boxplot(cfc$CellFractionAfterCorrection)
dev.copy(jpeg,filename="../data/myNormsAdditional/RA_450K/CF_RA_450K_CFC_Controls.jpg");
dev.off ();


# Save corrected myNorm
write.csv(cfc$CorrectedBeta, "../data/processed/CorrectedMyNorms/myNorm_cov2_cov407_hb_CFC.csv")
