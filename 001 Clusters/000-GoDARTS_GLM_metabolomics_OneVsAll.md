## Load packages

    library(datashieldclient)
    library(dsCDISCclient)
    library(survminer)
    library(survival)
    library(reshape2)
    library(viridis)

    knitr::opts_knit$set(root.dir = "D:/001_Projects/001_RHAPSODY/")

## Login

    source("./000.Final_Scripts/003.Utils/Utils.R")
    source("./000.Final_Scripts/003.Utils/Utils_GB.R")
    source("./000.Final_Scripts/003.Utils/Utils_CoxModels.R")
    datasource <- "opal"
    loginpath <- "logindata_template.txt"
    cohort <- "godarts"
    logIn(datasource, loginpath, cohort)

## Load data

    #Select required attributes 
    myAttributes <- c("METABOLOMICS","CLUSTER")
    myTransform <-c(3,0)
    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,transformVector=myTransform)

# Run models

Using the function `ds.glm` models for the 180 metabolites. metabolites
are `log10` transformed.

    # There are 29 metabolites in GoDARTS
    metabolites <- ds.colnames('ModelData')$godarts
    metabolites <- metabolites[which(metabolites == "LBTESTCD.AADA"):which(metabolites == "LBTESTCD.UDCA")]

    # SIDD comparisons
    GoDARTS.metabolomics.results.SIDD  <- getModelForClusters(metabolites, compCluster=2, refCluster="3,4,5,6")
    GoDARTS.metabolomics.results.SIRD  <- getModelForClusters(metabolites, compCluster=3, refCluster="2,4,5,6")
    GoDARTS.metabolomics.results.MOD  <- getModelForClusters(metabolites, compCluster=4, refCluster="2,3,5,6")
    GoDARTS.metabolomics.results.MARD  <- getModelForClusters(metabolites, compCluster=5, refCluster="2,3,4,6")
    GoDARTS.metabolomics.results.MARDH  <- getModelForClusters(metabolites, compCluster=6, refCluster="2,3,4,5")

    save(GoDARTS.metabolomics.results.SIDD, 
         GoDARTS.metabolomics.results.SIRD, 
         GoDARTS.metabolomics.results.MOD, 
         GoDARTS.metabolomics.results.MARD, 
         GoDARTS.metabolomics.results.MARDH, 
         file="./000.Final_Scripts/002.Clusters.Data/Metabolites GLM_models across clusters_FDB_GoDARTS_OneVsAll.RData")

    datashield.logout(opal)
