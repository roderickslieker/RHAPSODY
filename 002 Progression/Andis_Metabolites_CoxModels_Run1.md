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
    cohort <- "andis"
    datasource <- "opal"
    logIn(datasource, loginpath, cohort)

## Prepare Data

    datashield.assign(opal = opal, symbol = 'lb', value = 'rhapsody.LB')
    ds2.subset(symbol = 'lb',what = 'lb',row.filter = '-which(lb$LBRUNID == "Run1")')
    ds.dim("lb")

    ds2.subset('lb','lb',row.filter="-which(is.na(lb$LBORRES))")
    myAttributes <- c("METABOLOMICS","BMI","AGE","SEX","HBA1C","HDL","CPEP","DIABDUR","TIMEINSULIN","TIMEINSULINSTATUS","MEDICATION")
    myTransform <-c(3,0,0,0,0,0,0,0,0,0,0)

    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,transformVector=myTransform, plusOne = TRUE)
    ds2.subset(symbol = 'ModelData', what = 'ModelData', row.filter = 'which(!is.na(ModelData$LBTESTCD.TIMEINSULINSTATUS))')
    ds2.cut(x = 'LBTESTCD.HBA1C', new.name = 'HBA1C.CUT', df = 'ModelData', breaks = c(0,53,75,140), datasources = opal)
    ds.dim("ModelData")

## Prepare Formulas

    #MODEL 1 
    MetabolomicsFormula1 <- "survival::Surv(LBTESTCD.TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s + survival::strata(HBA1C.CUT) + age + sex + VSTESTCD.BMI"
    #MODEL 2
    MetabolomicsFormula2 <- "survival::Surv(LBTESTCD.TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s + survival::strata(HBA1C.CUT) + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEP "
    #MODEL 3
    MetabolomicsFormula3 <- ("survival::Surv(LBTESTCD.TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s + survival::strata(HBA1C.CUT) + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEP + DIABETESDURATION + glucoselowering")

## Prepare function call

    metabolites <- ds.colnames('ModelData')$andis
    metabolites <- metabolites[which(metabolites == "LBTESTCD.AADA"):which(metabolites=="LBTESTCD.Tyr")]

## Calculate Cox models for model 1

HbA1c, age, sex, bmi

    MetabolomicsModel1andis.disc <- CalculateCoxModel(vars = metabolites,
                                        formula = MetabolomicsFormula1,
                                        cohort = "andis",
                                        omicType = "METABOLOMICS",
                                        data = "ModelData")

    MetabolomicsModel2andis.disc <- CalculateCoxModel(vars = metabolites,
                                        formula = MetabolomicsFormula2,
                                        cohort = "andis",
                                        omicType = "METABOLOMICS",
                                        data = "ModelData")

    MetabolomicsModel3andis.disc <- CalculateCoxModel(vars = metabolites,
                                        formula = MetabolomicsFormula3,
                                        cohort = "andis",
                                        omicType = "METABOLOMICS",
                                        data = "ModelData")

    save(MetabolomicsModel1andis.disc, MetabolomicsModel2andis.disc, MetabolomicsModel3andis.disc, 
        file="./000.Final_Scripts/001.CoxModels.Data/ANDIS_Metabolites_CoxModel_1_adjusted_disc.RData")

## Logout

    logOut(datasource = datasource)
