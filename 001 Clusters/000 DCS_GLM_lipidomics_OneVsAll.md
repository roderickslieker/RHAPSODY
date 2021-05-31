## Load packages

    library(datashieldclient)
    library(dsCDISCclient)
    library(survminer)
    library(survival)
    library(reshape2)
    library(viridis)

    knitr::opts_knit$set(root.dir = "D:/001_Projects/001_RHAPSODY/")

## Login

Login on the DCS node only for now

    source("./000.Final_Scripts/003.Utils/Utils.R")
    source("./000.Final_Scripts/003.Utils/Utils_GB.R")

    logins <- read.table("logindata_template.txt", header=T)
    logins <- logins[logins$server %in% "dcs",]

    opal <- datashield.login(logins)
    opal

## Load data

    datashield.assign(opal, 'lb','rhapsody.LB')
    datashield.assign(opal, 'vs','rhapsody.VS')
    datashield.assign(opal, 'dm','rhapsody.DM')
    datashield.assign(opal, 'cm','rhapsody.CM')

    #Select required attributes 
    myAttributes <- c("LIPIDOMICS","CLUSTER")
    #Determine which attributes need to be transformed
    #   0 = none
    #   1 = log10(x+1)
    #   2 = scale(x)
    #   3 = scale(log10(x+1))
    myTransform <-c(3,0,0,0)

    #Prepare dataframe
    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,transformVector=myTransform)

# Run models

Using the function `ds.glm` models for the 183 lipids. lipids are
`log10` transformed and `scaled`.

    # There are 29 lipids in DCS
    lipids <- ds.colnames('ModelData')$dcs
    lipids <- lipids[which(lipids == "LBTESTCD.CE.14.0.0"):which(lipids == "LBTESTCD.TAG.58.8.0")]

    # SIDD comparisons
    DCS.lipidomics.results.SIDD  <- getModelForClusters(lipids, compCluster=2, refCluster="3,4,5,6")
    DCS.lipidomics.results.SIRD  <- getModelForClusters(lipids, compCluster=3, refCluster="2,4,5,6")
    DCS.lipidomics.results.MOD  <- getModelForClusters(lipids, compCluster=4, refCluster="2,3,5,6")
    DCS.lipidomics.results.MARD  <- getModelForClusters(lipids, compCluster=5, refCluster="2,3,4,6")
    DCS.lipidomics.results.MARDH  <- getModelForClusters(lipids, compCluster=6, refCluster="2,3,4,5")

    save(DCS.lipidomics.results.SIDD, 
         DCS.lipidomics.results.SIRD, 
         DCS.lipidomics.results.MOD, 
         DCS.lipidomics.results.MARD, 
         DCS.lipidomics.results.MARDH, 
         file="./000.Final_Scripts/002.Clusters.Data/Lipids GLM_models across clusters_FDB_DCS_OneVsAll.RData")
