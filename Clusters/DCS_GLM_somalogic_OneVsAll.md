## Load packages

    library(datashieldclient)
    library(dsCDISCclient)
    library(survminer)
    library(survival)
    library(reshape2)
    library(viridis)

    setwd("/exports/molepi/users/Roderick/Postdoc/RHAPSODY/FDB_Git/Omics_FDB")
    #knitr::opts_knit$set(root.dir = "/Users/roderick/Documents/Projecten en papers/Rhapsody_FDB/")
    #knitr::opts_knit$set(root.dir = "/exports/molepi/users/gbouland/")

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
    myAttributes <- c("SOMALOGIC","CLUSTER","TimeToSpin")
    myTransform <-c(3,0,0)
    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,
                transformVector=myTransform)

# Run models

    peptides <- ds.colnames('ModelData')$dcs
    peptides <- peptides[grep("SL0", peptides)]

    # SIDD comparisons
    DCS.peptidomics.results.SIDD  <- getModelForClusters(peptides, compCluster=2, refCluster="3,4,5,6", formula='ModelData.subset$%s~ModelData.subset$newCluster + ModelData.subset$LBTESTCD.TimeToSpin')
    DCS.peptidomics.results.SIRD  <- getModelForClusters(peptides, compCluster=3, refCluster="2,4,5,6", formula='ModelData.subset$%s~ModelData.subset$newCluster + ModelData.subset$LBTESTCD.TimeToSpin')
    DCS.peptidomics.results.MOD  <- getModelForClusters(peptides, compCluster=4, refCluster="2,3,5,6", formula='ModelData.subset$%s~ModelData.subset$newCluster + ModelData.subset$LBTESTCD.TimeToSpin')
    DCS.peptidomics.results.MARD  <- getModelForClusters(peptides, compCluster=5, refCluster="2,3,4,6", formula='ModelData.subset$%s~ModelData.subset$newCluster + ModelData.subset$LBTESTCD.TimeToSpin')
    DCS.peptidomics.results.MARDH  <- getModelForClusters(peptides, compCluster=6, refCluster="2,3,4,5", formula='ModelData.subset$%s~ModelData.subset$newCluster + ModelData.subset$LBTESTCD.TimeToSpin')

    save(DCS.peptidomics.results.SIDD, 
         DCS.peptidomics.results.SIRD, 
         DCS.peptidomics.results.MOD, 
         DCS.peptidomics.results.MARD, 
         DCS.peptidomics.results.MARDH, 
         file="./000.Final_Scripts/002.Clusters.Data/Peptides GLM_models across clusters_FDB_DCS_OneVsAll_adjusted.RData")

## Quick check, how many significant?

### SIDD

    DCS.peptidomics.results.SIDD[DCS.peptidomics.results.SIDD$p.adj <= 0.05,]
    table(DCS.peptidomics.results.SIDD$p.adj <= 0.05)

### SIRD

    DCS.peptidomics.results.SIRD[DCS.peptidomics.results.SIRD$p.adj <= 0.05,]
    table(DCS.peptidomics.results.SIRD$p.adj <= 0.05)

### MOD

    DCS.peptidomics.results.MOD[DCS.peptidomics.results.MOD$p.adj <= 0.05,]
    table(DCS.peptidomics.results.MOD$p.adj <= 0.05)

### MARD

    DCS.peptidomics.results.MARD[DCS.peptidomics.results.MARD$p.value <= 0.05,]
    table(DCS.peptidomics.results.MARD$p.adj <= 0.05)

### MARDH

    DCS.peptidomics.results.MARDH[DCS.peptidomics.results.MARDH$p.value <= 0.05,]
    table(DCS.peptidomics.results.MARDH$p.adj <= 0.05)
