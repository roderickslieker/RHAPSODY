## Load libraries and Utils scripts

    #library(datashieldclient)
    library(dsBaseClient)
    library(dsSwissKnifeClient)
    library(survminer)
    library(survival)
    library(reshape2)
    library(viridis)

    source("D:/001_Projects/033 Final scripts RHAPSODY/999 Utils/Utils_CoxModels.R")
    source("D:/001_Projects/033 Final scripts RHAPSODY/999 Utils/Utils_GB.R")
    source("D:/001_Projects/033 Final scripts RHAPSODY/999 Utils/Utils.R")


    datasource <- "opal"
    loginpath <- "logindata_template.txt"
    cohort <- "godarts"

    logins <- read.table(loginpath, header=T)
    logins <- logins[logins$server %in% c(cohort),]
    assign("opal",datashield.login(logins), envir = .GlobalEnv)

## Prepare Data

    myAttributes <- c("METABOLOMICS","LIPIDOMICS","SOMALOGIC","BMI","AGE","SEX","HDL","HBA1C","CPEPTIDE","TIMETOINSULIN", "TIMEINSULINSTATUS","DIABDUR")
    myTransform <-c(3,3,3,rep(0,9))
    prepareData(assign="ModelData",opal=opal,attributes=myAttributes,transformVector=myTransform, plusOne = TRUE, cohort = "godarts")

    dssCut(x = 'LBTESTCD.HBA1C', new.name = 'HBA1C.CUT', df = 'ModelData', breaks = c(0,53,75,140), datasources = opal)

## Prepare Formulas

    Formula1 <- "survival::Surv(as.numeric(TIMEINSULIN), LBTESTCD.TIMEINSULINSTATUS) ~ %s*VSTESTCD.BMI + survival::strata(HBA1C.CUT) + age + sex"

    Formula2 <- "survival::Surv(as.numeric(TIMEINSULIN), LBTESTCD.TIMEINSULINSTATUS) ~ %s*LBTESTCD.CPEPTIDE + survival::strata(HBA1C.CUT) + age + sex + VSTESTCD.BMI"

## Prepare function call

    vars <- c("Hcit","AADA","Ile","GCA","TCA",
    "TAG.50.1.0","TAG.46.1.0","TAG.46.2.0","TAG.48.1.0","TAG.51.1.0",
    "TAG.48.2.0","TAG.48.3.0","TAG.49.1.0","SM.42.2.2",
    "SL003869","SL005208","SL004152","SL012774","SL009045","SL002731",
    "SL000454","SL003733","SL004814","SL010384","SL018921")
    variables <- paste0("LBTESTCD.",vars)

## BMI interaction

    godarts.bmi.interaction <- lapply(variables, getCoxModelInteraction, form  = Formula1, cohort = "godarts", nameVar = "%s:VSTESTCD.BMI", data ="ModelData") %>% do.call(what=rbind)

## C-peptide interaction

    godarts.cpeptide.interaction <- lapply(variables, getCoxModelInteraction, form  = Formula2, cohort = "godarts", nameVar = "%s:LBTESTCD.CPEPTIDE", data ="ModelData") %>% do.call(what=rbind)

    save(godarts.bmi.interaction, godarts.cpeptide.interaction, file="Interaction_models_GoDarts.RData")

## Logout

    datashield.logout(opal)