## Load libraries and Utils scripts

    library(datashieldclient)
    library(dsCDISCclient)
    library(survminer)
    library(survival)

    #Load utils
    source("./000.Final_Scripts/003.Utils/Utils_GB.R")
    source("./000.Final_Scripts/003.Utils/Utils_CoxModels.R")

## Login

    loginpath <- "logindata_template.txt"
    cohort <- "godarts"
    datasource <- "opal"
    logIn(datasource, loginpath, cohort)

## Prepare Data

    #Select required attributes 
    myAttributes <- c("LIPIDOMICS","MEDICATION","BMI","AGE","SEX","HDL","HBA1C","CPEPTIDE","TIMETOINSULIN",
        "TIMEINSULINSTATUS","DIABDUR")
    myTransform <-c(3,0,0,0,0,0,0,0,0,0,0)
    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,transformVector=myTransform)

## Prepare Formulas

    #MODEL 1 
    LipoFormula1 <- "survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 survival::strata(cut(LBTESTCD.HBA1C ,breaks=c(0,53,75,140))) + age + sex + VSTESTCD.BMI "
    #MODEL 2
    LipoFormula2 <- "survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 survival::strata(cut(LBTESTCD.HBA1C ,breaks=c(0,53,75,140))) + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE "
    #MODEL 3
    LipoFormula3 <- ("survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 survival::strata(cut(LBTESTCD.HBA1C ,breaks=c(0,53,75,140))) + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE + DIABETESDURATION + glucoselowering")

## Prepare arguments for function call

    targets <- read.table("./Data/SwissLipids.csv", header=T, sep=",")
    lipids <- ds.colnames('ModelData')$godarts
    lipids <- lipids[which(lipids == "LBTESTCD.CE.14.0.0"):which(lipids=="LBTESTCD.TAG.58.8.0")]

## Calculate Cox models for model 1

HbA1c, Age, Sex, BMI

    LiposModel1GoDarts <- CalculateCoxModel(vars = lipids,
                                        formula = LipoFormula1,
                                        targets = targets,cohort = "godarts",
                                        omicType = "LIPIDOMICS",
                                        data = "ModelData")
    LiposModel1GoDarts$fdr <- p.adjust(LiposModel1GoDarts$p.val,method="fdr")

    save(LiposModel1GoDarts,file="lipotypegodartsCoxModel1.RData")
    knitr::kable(head(LiposModel1GoDarts[order(LiposModel1GoDarts$fdr,decreasing = F),],10))

## Calculate Cox models for model 2

HbA1c, Age, Sex, BMI, HDL , CPeptide

    LiposModel2GoDarts <- CalculateCoxModel(vars = lipids,
                                        formula = LipoFormula2,
                                        targets = targets,cohort = "godarts",
                                        omicType = "LIPIDOMICS",
                                        data = "ModelData")
    LiposModel2GoDarts$fdr <- p.adjust(LiposModel2GoDarts$p.val,method="fdr")

    save(LiposModel2GoDarts,file="lipotypegodartsCoxModel2.RData")
    knitr::kable(head(LiposModel2GoDarts[order(LiposModel2GoDarts$fdr,decreasing = F),],10))

## Calculate Cox models for model 3

HbA1c, Age, Sex, BMI, HDL , CPeptide, Diabetes Duration, Glucose
Lowering Drugs

    LiposModel3GoDarts <- CalculateCoxModel(vars = lipids,
                                        formula = LipoFormula3,
                                        targets = targets,cohort = "godarts",
                                        omicType = "LIPIDOMICS",
                                        data = "ModelData")
    LiposModel3GoDarts$fdr <- p.adjust(LiposModel3GoDarts$p.val,method="fdr")

    save(LiposModel3GoDarts,file="lipotypegodartsCoxModel3.RData")
    knitr::kable(head(LiposModel3GoDarts[order(LiposModel3GoDarts$fdr,decreasing = F),],10))

## Logout

    logOut(datasource)
