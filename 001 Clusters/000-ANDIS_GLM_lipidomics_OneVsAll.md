    library(datashieldclient)
    library(dsCDISCclient)
    library(survminer)
    library(survival)

    source("./000.Final_Scripts/003.Utils/Utils_GB.R")
    source("./000.Final_Scripts/003.Utils/Utils_CoxModels.R")
    loginpath <- "logindata_template.txt"
    cohort <- "andis"
    datasource <- "opal"
    logIn(datasource, loginpath, cohort)

## Load data

    myAttributes <- c("LIPIDOMICS","CLUSTER")
    myTransform <-c(3,0)
    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,transformVector=myTransform)

# Run models

Using the function `ds.glm` models for the 183 lipids. lipids are
`log10` transformed and `scaled`.

    # There are 29 lipids in ANDIS
    lipids <- ds.colnames('ModelData')$andis
    lipids <- lipids[which(lipids == "LBTESTCD.CE.14.0.0"):which(lipids == "LBTESTCD.TAG.58.8.0")]

    # SIDD comparisons
    ANDIS.lipidomics.results.SIDD  <- getModelForClusters(lipids, compCluster=2, refCluster="3,4,5,6")
    ANDIS.lipidomics.results.SIRD  <- getModelForClusters(lipids, compCluster=3, refCluster="2,4,5,6")
    ANDIS.lipidomics.results.MOD  <- getModelForClusters(lipids, compCluster=4, refCluster="2,3,5,6")
    ANDIS.lipidomics.results.MARD  <- getModelForClusters(lipids, compCluster=5, refCluster="2,3,4,6")
    ANDIS.lipidomics.results.MARDH  <- getModelForClusters(lipids, compCluster=6, refCluster="2,3,4,5")

    save(ANDIS.lipidomics.results.SIDD, 
         ANDIS.lipidomics.results.SIRD, 
         ANDIS.lipidomics.results.MOD, 
         ANDIS.lipidomics.results.MARD, 
         ANDIS.lipidomics.results.MARDH, 
         file="./000.Final_Scripts/002.Clusters.Data/Lipids GLM_models across clusters_FDB_ANDIS_OneVsAll.RData")
