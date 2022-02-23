#Launching the library
library("ChAMP")

# Set working directory
# setwd("")

myLoad <- champ.load(directory = ".",
                     method="champ",
                     methValue="B",
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     ProbeCutoff=0,
                     SampleCutoff=0.1,
                     detPcut=0.01,
                     filterBeads=TRUE,
                     beadCutoff=0.05,
                     filterNoCG=TRUE,
                     filterSNPs=TRUE,
                     population=NULL,
                     filterMultiHit=TRUE,
                     filterXY=TRUE,
                     force=TRUE,
                     arraytype="EPIC")

champ.QC(beta = myLoad$beta,
         pheno=myLoad$pd$Sample_Group,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         resultsDir="QC/")

myNorm <- champ.norm(beta=myLoad$beta,
                     rgSet=myLoad$rgSet,
                     mset=myLoad$mset,
                     resultsDir="Norm/",
                     method="BMIQ",
                     plotBMIQ=TRUE,
                     arraytype="EPIC",
                     cores=4)

champ.SVD(beta = myNorm,
          rgSet=NULL,
          pd=myLoad$pd,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="SVD/")

norm <- champ.runCombat(beta = myNorm, pd=myLoad$pd)

# save as csv file
write.table(myNorm, file = "", sep=",")


library("ChAMP")
