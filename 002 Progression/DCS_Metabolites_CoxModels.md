## Load libraries and Utils scripts

    library(datashieldclient)
    library(dsCDISCclient)
    library(survminer)
    library(survival)

    source("./000.Final_Scripts/003.Utils/Utils.R")
    source("./000.Final_Scripts/003.Utils/Utils_GB.R")
    source("./000.Final_Scripts/003.Utils/Utils_CoxModels.R")

## Login

    #loginpath <- "./logindata_template.txt"
    #Shark Gerard
    loginpath <- "logindata_template.txt"
    cohort <- "dcs"
    datasource <- "opal"
    logIn(datasource, loginpath, cohort)

## Prepare Data

    myAttributes <- c("METABOLOMICS","BMI","AGE","SEX","HDL","HBA1C","CPEPTIDE","TIMETOINSULIN",
        "TIMEINSULINSTATUS","DIABDUR")
    myTransform <-c(3,0,0,0,0,0,0,0,0,0,0)
    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,transformVector=myTransform, plusOne = TRUE)

## Prepare Formulas

    #MODEL 1 
    MetabolomicsFormula1 <- "survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 LBTESTCD.HBA1C + age + sex + VSTESTCD.BMI"
    #MODEL 2
    MetabolomicsFormula2 <- "survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 LBTESTCD.HBA1C + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE "
    #MODEL 3
    MetabolomicsFormula3 <- ("survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 LBTESTCD.HBA1C + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE + DIABETESDURATION + glucoselowering")

## Prepare function call

    metabolites <- ds.colnames('ModelData')$dcs
    metabolites <- metabolites[which(metabolites == "LBTESTCD.AADA"):which(metabolites=="LBTESTCD.UDCA")]

## Calculate Cox models for model 1

HbA1c, age, sex, bmi

    MetabolomicsModel1DCS <- CalculateCoxModel(vars = metabolites,
                                        formula = MetabolomicsFormula1,
                                        cohort = "dcs",
                                        omicType = "METABOLOMICS",
                                        data = "ModelData")


    MetabolomicsModel1DCS$fdr <- p.adjust(MetabolomicsModel1DCS$p.val,method="fdr")

    save(MetabolomicsModel1DCS,file="./000.Final_Scripts/001.CoxModels.Data/DCS_Metabolites_CoxModel_1.RData")
    knitr::kable(head(MetabolomicsModel1DCS[order(MetabolomicsModel1DCS$fdr,decreasing = F),],10))

## Calculate Cox models for model 2

HbA1c, age, sex, bmi, hdl , cpeptide

    MetabolomicsModel2DCS <- CalculateCoxModel(vars = metabolites,
                                        formula = MetabolomicsFormula2,
                                        cohort = "dcs",
                                        omicType = "METABOLOMICS",
                                        data = "ModelData")


    MetabolomicsModel2DCS$fdr <- p.adjust(MetabolomicsModel2DCS$p.val,method="fdr")

    save(MetabolomicsModel2DCS,file="./000.Final_Scripts/001.CoxModels.Data/DCS_Metabolites_CoxModel_2.RData")
    knitr::kable(head(MetabolomicsModel2DCS[order(MetabolomicsModel2DCS$fdr,decreasing = F),],10))

## Calculate Cox models for model 3

HbA1c, age, sex, bmi, hdl , cpeptide,diabetes duration,glucose lowering
drugs

    MetabolomicsModel3DCS <- CalculateCoxModel(vars = metabolites,
                                        formula = MetabolomicsFormula3,
                                        cohort = "dcs",
                                        omicType = "METABOLOMICS",
                                        data = "ModelData")


    MetabolomicsModel3DCS$fdr <- p.adjust(MetabolomicsModel3DCS$p.val,method="fdr")

    save(MetabolomicsModel3DCS,file="./000.Final_Scripts/001.CoxModels.Data/DCS_Metabolites_CoxModel_3.RData")
    knitr::kable(head(MetabolomicsModel3DCS[order(MetabolomicsModel3DCS$fdr,decreasing = F),],10))

## Logout

    logOut(datasource)
