## Load libraries and Utils scripts

    library(datashieldclient)
    library(dsCDISCclient)
    library(survminer)
    library(survival)

    #Load utils
    source("./000.Final_Scripts/003.Utils/Utils_GB.R")
    source("./000.Final_Scripts/003.Utils/Utils_CoxModels.R")

## Login

    logins <- read.table("./logindata_template.txt", header=T)
    logins <- logins[logins$server %in% c("dcs"),]
    opal <- datashield.login(logins)
    opal

## Prepare Data

    #Select required attributes 
    myAttributes <- c("SOMALOGIC","MEDICATION","BMI","AGE","SEX","HDL","HBA1C","CPEPTIDE","TIMETOINSULIN","TIMEINSULINSTATUS","DIABDUR","TimeToSpin","METFORMIN")
    #Determine which attributes need to be transformed
    #   0 = none
    #   1 = log10(x+1)
    #   2 = scale(x)
    #   3 = scale(log10(x+1))
    myTransform <-c(3,0,0,0,0,0,0,0,0,0,0,0,0)
    #Prepare dataframe
    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,transformVector=myTransform)

## Prepare Formulas

    #MODEL 1 
    PeptidesFormula1 <- "survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 LBTESTCD.HBA1C + age + sex + VSTESTCD.BMI + LBTESTCD.TimeToSpin"
    #MODEL 2
    PeptidesFormula2 <- "survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 LBTESTCD.HBA1C + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE + LBTESTCD.TimeToSpin "
    #MODEL 3
    PeptidesFormula3 <- ("survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ %s +
                                 LBTESTCD.HBA1C + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE + DIABETESDURATION + glucoselowering + LBTESTCD.TimeToSpin")

## Prepare

    load("/home/gabouland/CoxData/Conversion_Peptides.RData")
    targets <- conversion;rm(conversion)
    peptides <- ds.colnames('ModelData')$dcs
    peptides <- peptides[which(peptides == "LBTESTCD.SL000001"):which(peptides=="LBTESTCD.SL021043")]

## Calculate Cox models for model 1

HbA1c, Age, Sex, BMI, TimeToSpin

    PeptidesModel1 <- CalculateCoxModel(vars = peptides[1:600],
                                        formula = PeptidesFormula1,
                                        targets = targets,cohort = "dcs",
                                        omicType = "SOMALOGIC",
                                        data = "ModelData")

    PeptidesModel1.2 <- CalculateCoxModel(vars = peptides[601:1195],
                                        formula = PeptidesFormula1,
                                        targets = targets,cohort = "dcs",
                                        omicType = "SOMALOGIC",
                                        data = "ModelData")

    PeptidesModel1DCS <- rbind(PeptidesModel1,PeptidesModel1.2)
    PeptidesModel1DCS$fdr <- p.adjust(PeptidesModel1DCS$p.val,method="fdr")

    save(PeptidesModel1DCS,file="SomalogicDCSCoxModel1.RData")
    knitr::kable(head(PeptidesModel1DCS[order(PeptidesModel1DCS$fdr,decreasing = F),],10))

## Calculate Cox models for model 2

HbA1c, Age, Sex, BMI, HDL , CPeptide, Time to Spin

    PeptidesModel2 <- CalculateCoxModel(vars = peptides[1:600],
                                        formula = PeptidesFormula2,
                                        targets = targets,cohort = "dcs",
                                        omicType = "SOMALOGIC",
                                        data = "ModelData")

    PeptidesModel2.2 <- CalculateCoxModel(vars = peptides[601:1195],
                                        formula = PeptidesFormula2,
                                        targets = targets,cohort = "dcs",
                                        omicType = "SOMALOGIC",
                                        data = "ModelData")

    PeptidesModel2DCS <- rbind(PeptidesModel2,PeptidesModel2.2)
    PeptidesModel2DCS$fdr <- p.adjust(PeptidesModel2DCS$p.val,method="fdr")

    save(PeptidesModel2DCS,file="SomalogicDCSCoxModel2.RData")
    knitr::kable(head(PeptidesModel2DCS[order(PeptidesModel2DCS$fdr,decreasing = F),],10))

## Calculate Cox models for model 3

HbA1c, Age, Sex, BMI, HDL , CPeptide, Diabetes Duration, Glucose
Lowering Drugs, Time to Spin

    PeptidesModel3 <- CalculateCoxModel(vars = peptides[1:600],
                                        formula = PeptidesFormula3,
                                        targets = targets,cohort = "dcs",
                                        omicType = "SOMALOGIC",
                                        data = "ModelData")

    PeptidesModel3.2 <- CalculateCoxModel(vars = peptides[601:1195],
                                        formula = PeptidesFormula3,
                                        targets = targets,cohort = "dcs",
                                        omicType = "SOMALOGIC",
                                        data = "ModelData")

    PeptidesModel3DCS <- rbind(PeptidesModel3,PeptidesModel3.2)
    PeptidesModel3DCS$fdr <- p.adjust(PeptidesModel3DCS$p.val,method="fdr")

    save(PeptidesModel3DCS,file="SomalogicDCSCoxModel3.RData")
    knitr::kable(head(PeptidesModel3DCS[order(PeptidesModel3DCS$fdr,decreasing = F),],10))

## Sensitivity

    ds.replaceNA('ModelData$metformin', forNA = '0', datasources = opal)
    ds.dataframe('metformin.noNA', newobj = 'metformin.noNA')
    ds.cbind(c('ModelData', 'metformin.noNA'), newobj = 'ModelData2', datasources = opal) 
    #MODEL 3
    fit <- ds2.coxph(formula = survival::Surv(TIMEINSULIN, LBTESTCD.TIMEINSULINSTATUS) ~ LBTESTCD.SL003869 +
                                 LBTESTCD.HBA1C + age + sex + VSTESTCD.BMI + LBTESTCD.HDL + LBTESTCD.CPEPTIDE + DIABETESDURATION + glucoselowering + metformin.noNA + LBTESTCD.TimeToSpin,data="ModelData2", async=F)

## Logout
