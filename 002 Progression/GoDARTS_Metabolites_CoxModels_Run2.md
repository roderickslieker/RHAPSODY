## Load libraries and Utils scripts

    .libPaths("/mnt/d/000_Packages/Linux")
    library(datashieldclient)
    library(dsCDISCclient)
    library(survminer)
    library(survival)

    #Load utils
    source("./000.Final_Scripts/003.Utils/Utils_GB.R")
    source("./000.Final_Scripts/003.Utils/Utils_CoxModels.R")

## Login

    loginpath <- "./logindata_template.txt"
    #Shark Gerard
    #loginpath <- "/exports/molepi/users/gbouland/Login/logindata_template"
    cohort <- "godarts"
    datasource <- "opal"
    logIn(datasource, loginpath, cohort)

## Prepare Data

    datashield.assign(opal = opal, symbol = 'lb', value = 'rhapsody.LB')
    ds2.subset(symbol = 'lb',what = 'lb',row.filter = '-which(lb$LBRUNID == "1")')
    ds.dim("lb")

    ds2.subset('lb','lb',row.filter="-which(is.na(lb$LBORRES))")
    myAttributes <- c("METABOLOMICS","BMI","AGE","SEX","HBA1C","HDL","CPEPTIDE","DIABDUR","TIMEINSULIN","TIMEINSULINSTATUS","MEDICATION")
    myTransform <-c(3,0,0,0,0,0,0,0,0,0,0)

    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,transformVector=myTransform, plusOne = TRUE)
    ds2.subset(symbol = 'ModelData', what = 'ModelData', row.filter = 'which(!is.na(ModelData$LBTESTCD.TIMEINSULINSTATUS))')
    ds2.cut(x = 'LBTESTCD.HBA1C', new.name = 'HBA1C.CUT', df = 'ModelData', breaks = c(0,53,75,140), datasources = opal)
    ds.dim("ModelData")

## Prepare Formulas

    #MODEL 1 
    MetabolomicsFormula1 <- "survival::Surv(LBTESTCD.TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s + survival::strata(HBA1C.CUT) + age + sex + VSTESTCD.BMI"
    #MODEL 2
    MetabolomicsFormula2 <- "survival::Surv(LBTESTCD.TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s + survival::strata(HBA1C.CUT) + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE"
    #MODEL 3
    MetabolomicsFormula3 <- ("survival::Surv(LBTESTCD.TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s + survival::strata(HBA1C.CUT) + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE + DIABETESDURATION + glucoselowering")

## Prepare function call

    metabolites <- ds.colnames('ModelData')$godarts
    metabolites <- metabolites[which(metabolites == "LBTESTCD.AADA"):which(metabolites=="LBTESTCD.UDCA")]

## Calculate Cox models for model 1

HbA1c, age, sex, bmi

    MetabolomicsModel1godarts.repl <- CalculateCoxModel(vars = metabolites,
                                        formula = MetabolomicsFormula1,
                                        cohort = "godarts",
                                        omicType = "METABOLOMICS",
                                        data = "ModelData")

    MetabolomicsModel2godarts.repl <- CalculateCoxModel(vars = metabolites,
                                        formula = MetabolomicsFormula2,
                                        cohort = "godarts",
                                        omicType = "METABOLOMICS",
                                        data = "ModelData")

    MetabolomicsModel3godarts.repl <- CalculateCoxModel(vars = metabolites,
                                        formula = MetabolomicsFormula3,
                                        cohort = "godarts",
                                        omicType = "METABOLOMICS",
                                        data = "ModelData")

    save(MetabolomicsModel1godarts.repl, MetabolomicsModel2godarts.repl, MetabolomicsModel3godarts.repl, 
        file="./000.Final_Scripts/001.CoxModels.Data/GoDARTS_Metabolites_CoxModel_1_adjusted_repl.RData")

## Logout

    logOut(datasource = datasource)
