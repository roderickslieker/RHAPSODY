## Load packages

    library(datashieldclient)
    library(dsCDISCclient)
    library(survminer)
    library(survival)
    library(reshape2)
    library(viridis)

    setwd("D:/001_Projects/001_RHAPSODY/")

## Login

    source("./000.Final_Scripts/003.Utils/Utils.R")
    source("./000.Final_Scripts/003.Utils/Utils_GB.R")
    source("./000.Final_Scripts/003.Utils/Utils_CoxModels.R")
    datasource <- "opal"
    loginpath <- "logindata_template.txt"
    cohort <- "dcs"
    logIn(datasource, loginpath, cohort)

## Load data

    #Select required attributes 
    myAttributes <- c("METABOLOMICS","CLUSTER")
    myTransform <-c(3,0)
    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,transformVector=myTransform)

# Run models

Using the function `ds.glm` models for the 29 metabolites. Metabolites
are `log10` transformed and `scaled`.

    metabolites <- ds.colnames('ModelData')$dcs
    metabolites <- metabolites[which(metabolites == "LBTESTCD.AADA"):which(metabolites == "LBTESTCD.UDCA")]

    # SIDD comparisons
    DCS.metabolomics.results.SIDD  <- getModelForClusters(metabolites, compCluster=2, refCluster="3,4,5,6")
    DCS.metabolomics.results.SIRD  <- getModelForClusters(metabolites, compCluster=3, refCluster="2,4,5,6")
    DCS.metabolomics.results.MOD  <- getModelForClusters(metabolites, compCluster=4, refCluster="2,3,5,6")
    DCS.metabolomics.results.MARD  <- getModelForClusters(metabolites, compCluster=5, refCluster="2,3,4,6")
    DCS.metabolomics.results.MARDH  <- getModelForClusters(metabolites, compCluster=6, refCluster="2,3,4,5")

    save(DCS.metabolomics.results.SIDD, 
         DCS.metabolomics.results.SIRD, 
         DCS.metabolomics.results.MOD, 
         DCS.metabolomics.results.MARD, 
         DCS.metabolomics.results.MARDH, 
         file="./000.Final_Scripts/002.Clusters.Data/Metabolites GLM_models across clusters_FDB_DCS_OneVsAll.RData")
