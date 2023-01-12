## Load libraries and Utils scripts

    library(datashieldclient)
    library(dsCDISCclient)
    library(survminer)
    library(survival)

    #Load utils
    source("./000.Final_Scripts/003.Utils/Utils_GB.R")
    source("./000.Final_Scripts/003.Utils/Utils_CoxModels.R")

## Login

    loginpath <- "./logindata_template.txt"
    cohort <- "godarts"
    datasource <- "opal"
    logIn(datasource, loginpath, cohort)

## Prepare Data

    datashield.assign(opal = opal, symbol = 'lb', value = 'rhapsody.LB')
    ds2.subset(symbol = 'lb',what = 'lb',row.filter = '-which(lb$LBRUNID == "2")')
    ds.dim("lb")

    ds2.subset('lb','lb',row.filter="-which(is.na(lb$LBORRES))")
    myAttributes <- c("METABOLOMICS","BMI","AGE","SEX","HBA1C","HDL","CPEPTIDE","DIABDUR","TIMEINSULIN","TIMEINSULINSTATUS","MEDICATION")
    myTransform <-c(3,0,0,0,0,0,0,0,0,0,0)

    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,transformVector=myTransform, plusOne = TRUE)
    cleanUp(myAttributes)
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

    MetabolomicsModel1godarts.disc <- CalculateCoxModel(vars = metabolites,
                                        formula = MetabolomicsFormula1,
                                        cohort = "godarts",
                                        omicType = "METABOLOMICS",
                                        data = "ModelData")

    MetabolomicsModel2godarts.disc <- CalculateCoxModel(vars = metabolites,
                                        formula = MetabolomicsFormula2,
                                        cohort = "godarts",
                                        omicType = "METABOLOMICS",
                                        data = "ModelData")

    MetabolomicsModel3godarts.disc <- CalculateCoxModel(vars = metabolites,
                                        formula = MetabolomicsFormula3,
                                        cohort = "godarts",
                                        omicType = "METABOLOMICS",
                                        data = "ModelData")

    save(MetabolomicsModel1godarts.disc, MetabolomicsModel2godarts.disc, MetabolomicsModel3godarts.disc, 
        file="./000.Final_Scripts/001.CoxModels.Data/GoDARTS_Metabolites_CoxModel_1_adjusted_disc.RData")

## Logout

    logOut(datasource = datasource)
