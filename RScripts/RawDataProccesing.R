#Launching the library
library("ChAMP")

# Set working directory
setwd("../data/raw/")

myLoad <- champ.load(directory = "CONCATED_Spain_HB/",
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
         resultsDir="../interim/Spain/")

myNorm <- champ.norm(beta=myLoad$beta,
                     rgSet=myLoad$rgSet,
                     mset=myLoad$mset,
                     resultsDir="../interim/Spain/",
                     method="BMIQ",
                     plotBMIQ=TRUE,
                     arraytype="EPIC",
                     cores=2)

champ.SVD(beta = myNorm,
          rgSet=NULL,
          pd=myLoad$pd,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="../SVD/")

# save as csv file
write.table(myNorm, file = "../interim/Spain/Spain_mynorm.csv", sep=",")
