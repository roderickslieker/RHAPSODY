## Load packages

    library(datashieldclient)
    library(dsBaseClient)
    library(dsSwissKnifeClient)
    library(survminer)
    library(survival)
    library(reshape2)
    library(viridis)

    #setwd("D:/001_Projects/001_RHAPSODY/")
    knitr::opts_knit$set(root.dir = "D:/001_Projects/001_RHAPSODY/")

## Login

    source("D:/001_Projects/001_RHAPSODY/001.Scripts/001 New utils/Utils_CoxModels.R")
    source("D:/001_Projects/001_RHAPSODY/001.Scripts/001 New utils/Utils_GB.R")
    source("D:/001_Projects/001_RHAPSODY/001.Scripts/001 New utils/Utils.R")

    datasource <- "opal"
    loginpath <- "logindata_template.txt"
    cohort <- "godarts"

    logins <- read.table(loginpath, header=T)
    logins <- logins[logins$server %in% c(cohort),]
    assign("opal",datashield.login(logins), envir = .GlobalEnv)

## Load data

    #Select required attributes 
    myAttributes <- c("SOMALOGIC","CLUSTER","SEX","AGE","BMI","HBA1C","HDL","LDL","TRIG","CPEPTIDE","TimeToSpin")
    myTransform <-c(3,rep(0, length(myAttributes)-1))
    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,cohort = "godarts",
                transformVector=myTransform)

# Run models

    peptides <- ds.colnames('ModelData')$godarts
    peptides <- peptides[grep("SL0", peptides)]

    # SIDD comparisons
    GoDARTS.pept.SIDD.adj  <- getModelForClusters(peptides, compCluster=2, refCluster="3,4,5,6", formula='ModelData.subset$%s~ModelData.subset$newCluster + ModelData.subset$LBTESTCD.TimeToSpin+ ModelData.subset$LBTESTCD.HBA1C')

    GoDARTS.pept.SIRD.adj  <- getModelForClusters(peptides, compCluster=3, refCluster="2,4,5,6", formula='ModelData.subset$%s~ModelData.subset$newCluster + ModelData.subset$LBTESTCD.TimeToSpin+ ModelData.subset$LBTESTCD.CPEPTIDE')

    GoDARTS.pept.MOD.adj  <- getModelForClusters(peptides, compCluster=4, refCluster="2,3,5,6", formula='ModelData.subset$%s~ModelData.subset$newCluster + ModelData.subset$LBTESTCD.TimeToSpin + ModelData.subset$VSTESTCD.BMI')

    GoDARTS.pept.MDH.adj   <- getModelForClusters(peptides, compCluster=6, refCluster="2,3,4,5", formula='ModelData.subset$%s~ModelData.subset$newCluster + ModelData.subset$LBTESTCD.TimeToSpin + ModelData.subset$LBTESTCD.HDL')

    GoDARTS.pept.MDH.adj_TRIG   <- getModelForClusters(peptides, compCluster=6, refCluster="2,3,4,5", formula='ModelData.subset$%s~ModelData.subset$newCluster + ModelData.subset$LBTESTCD.TimeToSpin + ModelData.subset$LBTESTCD.HDL + ModelData.subset$LBTESTCD.TRIG')

    save(GoDARTS.pept.SIDD.adj, 
         GoDARTS.pept.SIRD.adj, 
         GoDARTS.pept.MOD.adj, 
         GoDARTS.pept.MDH.adj, 
         GoDARTS.pept.MDH.adj_TRIG, 
         file="./000.Final_Scripts/002.Clusters.Data/Peptides GLM_models across clusters_FDB_GoDARTS_OneVsAll_adjusted.RData")
