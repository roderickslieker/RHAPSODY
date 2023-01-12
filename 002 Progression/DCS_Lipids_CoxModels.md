## Load libraries and Utils scripts

    library(datashieldclient)
    library(dsCDISCclient)
    library(survminer)
    library(survival)

    #Load utils
    source("/home/gabouland/CoxUtils/Utils_GB.R")
    source("/home/gabouland/CoxUtils/Utils_CoxModels.R")

## Login

    logins <- read.table("/home/gabouland/CoxLogin/logindata_template", header=T)
    logins <- logins[logins$server %in% c("dcs"),]
    opal <- datashield.login(logins)
    opal

## Load Tables

    datashield.assign(opal, 'lb','rhapsody.LB')
    datashield.assign(opal, 'vs','rhapsody.VS')
    datashield.assign(opal, 'dm','rhapsody.DM')
    datashield.assign(opal, 'cm','rhapsody.CM')

## Prepare Data

    #Select required attributes 
    myAttributes <- c("LIPIDOMICS","MEDICATION","BMI","AGE","SEX","HDL","HBA1C","CPEPTIDE","TIMETOINSULIN",
        "TIMEINSULINSTATUS","DIABDUR")
    #Determine which attributes need to be transformed
    #   0 = none
    #   1 = log10(x+1)
    #   2 = scale(x)
    #   3 = scale(log10(x+1))
    myTransform <-c(3,0,0,0,0,0,0,0,0,0,0)
    #Prepare dataframe
    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,removeNA=FALSE,transformVector=myTransform,cohort="dcs")

## Prepare Formulas

    #MODEL 1 
    LipoFormula1 <- "survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 LBTESTCD.HBA1C + age + sex + VSTESTCD.BMI "
    #MODEL 2
    LipoFormula2 <- "survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 LBTESTCD.HBA1C + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE "
    #MODEL 3
    LipoFormula3 <- ("survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 LBTESTCD.HBA1C + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE + DIABETESDURATION + glucoselowering")

## Prepare arguments for function call

    targets <- read.table("/home/gabouland/CoxData/SwissLipids.csv", header=T, sep=",")
    lipids <- ds.colnames('ModelData')$dcs
    lipids <- lipids[which(lipids == "LBTESTCD.CE.14.0.0"):which(lipids=="LBTESTCD.TAG.58.8.0")]

## Calculate Cox models for model 1

HbA1c, Age, Sex, BMI

    LiposModel1dcs <- CalculateCoxModel(vars = lipids,
                                        formula = LipoFormula1,
                                        targets = targets,cohort = "dcs",
                                        omicType = "LIPIDOMICS",
                                        data = "ModelData")
    LiposModel1dcs$fdr <- p.adjust(LiposModel1dcs$p.val,method="fdr")

    save(LiposModel1dcs,file="lipotypedcsCoxModel1.RData")
    knitr::kable(head(LiposModel1dcs[order(LiposModel1dcs$fdr,decreasing = F),],10))

## Calculate Cox models for model 2

HbA1c, Age, Sex, BMI, HDL , CPeptide

    LiposModel2dcs <- CalculateCoxModel(vars = lipids,
                                        formula = LipoFormula2,
                                        targets = targets,cohort = "dcs",
                                        omicType = "LIPIDOMICS",
                                        data = "ModelData")
    LiposModel2dcs$fdr <- p.adjust(LiposModel2dcs$p.val,method="fdr")

    save(LiposModel2dcs,file="lipotypedcsCoxModel2.RData")
    knitr::kable(head(LiposModel2dcs[order(LiposModel2dcs$fdr,decreasing = F),],10))

## Calculate Cox models for model 3

HbA1c, Age, Sex, BMI, HDL , CPeptide, Diabetes Duration, Glucose
Lowering Drugs

    LiposModel3dcs <- CalculateCoxModel(vars = lipids,
                                        formula = LipoFormula3,
                                        targets = targets,cohort = "dcs",
                                        omicType = "LIPIDOMICS",
                                        data = "ModelData")
    LiposModel3dcs$fdr <- p.adjust(LiposModel3dcs$p.val,method="fdr")

    save(LiposModel3dcs,file="lipotypedcsCoxModel3.RData")
    knitr::kable(head(LiposModel3dcs[order(LiposModel3dcs$fdr,decreasing = F),],10))

## Logout

    datashield.logout(opal)
