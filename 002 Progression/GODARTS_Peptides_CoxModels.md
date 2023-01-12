## Load libraries and Utils scripts

    library(datashieldclient)
    library(dsCDISCclient)
    library(survminer)
    library(survival)

    #Load utils
    source("./Utils/Utils_GB.R")
    source("./Utils/Utils_CoxModels.R")

## Login

    loginpath <- "./logindata_template.txt"
    cohort <- "godarts"
    datasource <- "opal"
    logIn(datasource, loginpath, cohort)

## Prepare Data

    #Select required attributes 
    myAttributes <- c("SOMALOGIC","MEDICATION","BMI","AGE","SEX","HDL","HBA1C","CPEPTIDE","TIMETOINSULIN",
        "TIMEINSULINSTATUS","DIABDUR","TimeToSpin")
    myTransform <-c(3,0,0,0,0,0,0,0,0,0,0,0)
    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,transformVector=myTransform)

## Prepare Formulas

    #MODEL 1 
    PeptidesFormula1 <- "survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 survival::strata(cut(LBTESTCD.HBA1C ,breaks=c(0,53,75,140))) + age + sex + VSTESTCD.BMI + LBTESTCD.TimeToSpin "
    #MODEL 2
    PeptidesFormula2 <- "survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 survival::strata(cut(LBTESTCD.HBA1C ,breaks=c(0,53,75,140))) + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE + LBTESTCD.TimeToSpin"
    #MODEL 3
    PeptidesFormula3 <- ("survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 survival::strata(cut(LBTESTCD.HBA1C ,breaks=c(0,53,75,140))) + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE + DIABETESDURATION + glucoselowering + LBTESTCD.TimeToSpin")

## Prepare arguments for function call

    load("./Data/Conversion_Peptides.RData")
    targets <- conversion;rm(conversion)
    peptides <- ds.colnames('ModelData')$godarts
    peptides <- peptides[which(peptides == "LBTESTCD.SL000001"):which(peptides=="LBTESTCD.SL021043")]

## Calculate Cox models for model 1

HbA1c, Age, Sex, BMI, TimeToSpin

    PeptidesModel1 <- CalculateCoxModel(vars = peptides[1:600],
                                        formula = PeptidesFormula1,
                                        targets = targets,cohort = "godarts",
                                        omicType = "SOMALOGIC",
                                        data = "ModelData")

    PeptidesModel1.2 <- CalculateCoxModel(vars = peptides[601:1195],
                                        formula = PeptidesFormula1,
                                        targets = targets,cohort = "godarts",
                                        omicType = "SOMALOGIC",
                                        data = "ModelData")

    PeptidesModel1GODARTS <- rbind(PeptidesModel1,PeptidesModel1.2)
    PeptidesModel1GODARTS$fdr <- p.adjust(PeptidesModel1GODARTS$p.val,method="fdr")

    save(PeptidesModel1GODARTS,file="SomalogicgodartsCoxModel1.RData")
    knitr::kable(head(PeptidesModel1GODARTS[order(PeptidesModel1GODARTS$fdr,decreasing = F),],10))

## Calculate Cox models for model 2

HbA1c, Age, Sex, BMI, HDL , CPeptide, Time to Spin

    PeptidesModel2 <- CalculateCoxModel(vars = peptides[1:600],
                                        formula = PeptidesFormula2,
                                        targets = targets,cohort = "godarts",
                                        omicType = "SOMALOGIC",
                                        data = "ModelData")

    PeptidesModel2.2 <- CalculateCoxModel(vars = peptides[601:1195],
                                        formula = PeptidesFormula2,
                                        targets = targets,cohort = "godarts",
                                        omicType = "SOMALOGIC",
                                        data = "ModelData")

    PeptidesModel2GODARTS <- rbind(PeptidesModel2,PeptidesModel2.2)
    PeptidesModel2GODARTS$fdr <- p.adjust(PeptidesModel2GODARTS$p.val,method="fdr")

    save(PeptidesModel2GODARTS,file="SomalogicgodartsCoxModel2.RData")
    knitr::kable(head(PeptidesModel2GODARTS[order(PeptidesModel2GODARTS$fdr,decreasing = F),],10))

## Calculate Cox models for model 3

HbA1c, Age, Sex, BMI, HDL , CPeptide, Diabetes Duration, Glucose
Lowering Drugs, Time to Spin

    PeptidesModel3 <- CalculateCoxModel(vars = peptides[1:600],
                                        formula = PeptidesFormula3,
                                        targets = targets,cohort = "godarts",
                                        omicType = "SOMALOGIC",
                                        data = "ModelData")

    PeptidesModel3.2 <- CalculateCoxModel(vars = peptides[601:1195],
                                        formula = PeptidesFormula3,
                                        targets = targets,cohort = "godarts",
                                        omicType = "SOMALOGIC",
                                        data = "ModelData")

    PeptidesModel3GODARTS <- rbind(PeptidesModel3,PeptidesModel3.2)
    PeptidesModel3GODARTS$fdr <- p.adjust(PeptidesModel3GODARTS$p.val,method="fdr")

    save(PeptidesModel3GODARTS,file="SomalogicgodartsCoxModel3.RData")
    knitr::kable(head(PeptidesModel3GODARTS[order(PeptidesModel3GODARTS$fdr,decreasing = F),],10))

## Logout

    logOut(datasource)
