dataStructure <- list()
#############################################################
###########################LB################################
#############################################################
dataStructure[['LIPIDOMICS']] <- c('lb','LBMETHOD','LIPIDOMICS')
dataStructure[['LIPIDOMICS_AGGREGATION']] <- c('lb','LBMETHOD', 'LIPIDOMICS_AGGREGATION')
dataStructure[['ELISA']] <- c('lb','LBMETHOD', 'ELISA')
dataStructure[['FATTY_ACIDS']] <- c('lb','LBMETHOD', 'FATTY_ACIDS')
dataStructure[['SOMALOGIC']] <- c('lb','LBMETHOD', 'SOMALOGIC')
dataStructure[['TimeToDecant']] <- c('lb','LBMETHOD', 'SOMALOGIC_COVAR','LBTESTCD',"TimeToDecant")
dataStructure[['TimeToSpin']] <- c('lb','LBMETHOD', 'SOMALOGIC_COVAR','LBTESTCD','TimeToSpin')
dataStructure[['CPEPTIDE']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','CPEPTIDE')
dataStructure[['MIRNA']] <- c('lb','LBMETHOD','MIRNA')
dataStructure[['METABOLOMICS']] <- c('lb','LBMETHOD', 'METABOLOMICS')
dataStructure[['TIMEINSULIN']]<- c('lb','LBMETHOD','ADDITIONAL','LBTESTCD','TIMEINSULIN')
dataStructure[['CLINICAL']] <- c('lb','LBMETHOD', 'CLINICAL')

dataStructure[['TIMEINSULINSTATUS']]<- c('lb','LBMETHOD','ADDITIONAL','LBTESTCD','TIMEINSULINSTATUS')
#dataStructure[['TIMEINSULIN ']]<- c('lb','LBMETHOD','ADDITIONAL','LBTESTCD','TIMEINSULIN ')
dataStructure[['CLUSTER']]<- c('lb','LBMETHOD','ADDITIONAL','LBTESTCD','CLUSTER')
dataStructure[['CREAT']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','CREAT')
dataStructure[['GLUC']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','GLUC')
dataStructure[['HBA1C']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','HBA1C')
dataStructure[['HDL']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','HDL')
dataStructure[['LDL']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','LDL')
dataStructure[['TRIG']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','TRIG')
dataStructure[['GADAB']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','GADAB')
dataStructure[['CPEP']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','CPEP')
dataStructure[['HOMA2B']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','HOMA2B')
dataStructure[['HOMA2IR']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','HOMA2IR')
dataStructure[['HOMA2S']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','HOMA2S')

#############################################################
###########################GD################################
#############################################################
dataStructure[['GWAS']] <- c('gd','LBMETHOD','GWAS')

#############################################################
###########################VS################################
#############################################################
dataStructure[['BMI']] <- c("vs","VSTESTCD","BMI")
dataStructure[['DIABP']] <- c("vs","VSTESTCD","DIABP")
dataStructure[['SYSBP']] <- c("vs","VSTESTCD","SYSBP")
dataStructure[['HEIGHT']] <- c("vs","VSTESTCD","HEIGHT")
dataStructure[['WEIGHT']] <- c("vs","VSTESTCD","WEIGHT")
#############################################################
###########################DM################################
#############################################################
dataStructure[['AGE']] <- c("dm")
dataStructure[['SEX']] <- c("dm")
#############################################################
###########################CM################################
#############################################################
dataStructure[['MEDICATION']] <- c("cm")
dataStructure[['MEDICATIONEXT']] <- c("cm")
dataStructure[['METFORMIN']] <- c("cm")
dataStructure[['TIMETOINSULIN']]<- c('lb')
dataStructure[['DIABDUR']] <- c("dm")
#############################################################


#' Prepares data based on selected attributes
#' @param assign New name of combined data object
#' @param opal Input is an opal object
#' @param attributes Desired columns
#' @return On the federated database, the desired data.frame is generated
#' @examples
#' prepareData("ModelData", opal, myAttributes)

dataframeSelector <- function(dataframeTag){
	print(dataframeTag)
  if(dataframeTag == "lb"){
    formulaDate <- "LBDTC"
    formulaOrres <- "LBORRES"
    formulaTest <- "LBTESTCD"
  }else if (dataframeTag == "vs"){
    formulaDate <- "VSDTC"
    formulaOrres <- "VSORRES"
    formulaTest <- "VSTESTCD"
  }else if (dataframeTag == "dm"){
    formulaDate <- "DMDTC"
    formulaOrres <- "none"
    formulaTest <- "none"
  }else if (dataframeTag == "gd"){
    formulaDate <- "LBDTC"
    formulaOrres <- "LBORRES"
    formulaTest <- "LBTESTCD"
  }else if(dataframeTag == "cm"){
  	formulaDate <- "CMSTDTC"
    formulaOrres <- "LBORRES"
    formulaTest <- "LBTESTCD"
  }
  return(c(formulaDate,formulaOrres,formulaTest))
}


diabetes.duration <- function(omicAttribute,opal){
	dssSubset(symbol = 'tmpdf', what = omicAttribute, col.filter = 'c("SUBJID","LBDTC")',datasources= opal)
	dssSubset(symbol = "DIABDUR", what = 'dm',row.filter = 'match(tmpdf$SUBJID, dm$SUBJID)',datasources= opal)
	dssSubset(symbol = "DIABDUR", what = 'DIABDUR',col.filter = 'c("RFSTDTC")',datasources= opal)
	ds.cbind(x = c("tmpdf","DIABDUR"), newobj = "DIABDUR",datasources = opal)
	dssAddDaysInterval(newobj = 'DIABDUR',df = 'DIABDUR',
		description.list =  list(DIABETESDURATION = list(start_date = 'RFSTDTC', end_date = 'LBDTC')), datasources = opal)

	dssSubset(symbol = 'neg.diab', what = 'DIABDUR', row.filter = "DIABDUR$DIABETESDURATION <= 0", datasources = opal)
	nrow.zero <- ds.dim('neg.diab', datasources = opal)[[1]][1]
  	if(nrow.zero != 0) # Because the date of sampling or diagnosis can be wrongly recorded.
  	{
    	dssSubset(symbol = 'non.neg.diab', what = 'DIABDUR', row.filter = "DIABDUR$DIABETESDURATION > 0", datasources = opal)

    	ds.vectorCalc(x = c('neg.diab$SUBJID','neg.diab$SUBJID'), calc = "-", newobj = 'DIABETESDURATION', datasources = opal)
    	dssSubset(symbol= 'neg.diab', what='neg.diab',col.filter='c("SUBJID","LBDTC","RFSTDTC")',datasources=opal)
    	ds.cbind(x = c("neg.diab","DIABETESDURATION"), newobj = "neg.diab",datasources = opal)
    	dssRbind(newobj = 'DIABDUR', x = 'non.neg.diab', y='neg.diab',datasources = opal)
    	dssSubset(symbol = "DIABDUR", what = 'DIABDUR',
    		row.filter = sprintf('match(%s$SUBJID, DIABDUR$SUBJID)',omicAttribute),datasources= opal)
  	}
  	#Organize data frame
  	ds.assign(toAssign = "DIABDUR$SUBJID", newobj = 'SUBJID', datasources = opal)
  	ds.assign(toAssign = "DIABDUR$LBDTC", newobj = 'DMDTC', datasources = opal)  	
  	ds.assign(toAssign = as.symbol("DIABDUR$DIABETESDURATION /365"), newobj = 'DIABETESDURATION', datasources = opal)
  	ds.dataFrame(x = c("SUBJID","DMDTC","DIABETESDURATION"), newobj = 'DIABDUR', datasources = opal)
}

make.row.FilterOmics.cm <- function(omics, delta)
{
     row.filter <- paste0('{
                    Reduce(c, lapply(',omics,'$SUBJID, function(sid)
                    {
                    # Select omics date
                    soma.sel <- unique(',omics,'[',omics,'$SUBJID %in% sid,"LBDTC"])[1]
                    sel.date <- as.Date(as.character(soma.sel))

                    # Substract 3 months period
                    sel.date.new <- sel.date - ',delta,'

                    # Restore CM date
                    idx <- which(cm$SUBJID == sid)
                    if(length(idx) == 0)
                    {
                         return(NA)
                         }else{
                         # Idx criterion
                         cm.sel <- cm[idx,]
                         # Date criterion
                         cm.sel$CMDTC <- as.Date(as.character(cm.sel$CMDTC))

                         idx.out <- idx[which(cm.sel$CMDTC <= sel.date.new)]
                         return(idx.out)                  
                     }
                    }))
                    }')
     return(row.filter)
}


medication.module.godarts <- function(omics, delta, opals)
{  
  if(is.null(omics))
  {
    stop("Please provide a object name")
  }
  if(is.null(delta))
  {
    stop("Please provide delta")
  }
  row.filter.soma <- make.row.FilterOmics.cm(omics =omics , delta= delta)
  
  #Subset Somalogic
  dssSubset(symbol = "medication", what = "cm", row.filter = row.filter.soma, datasources=opals)
  
  # Transform
  dssPivot("medication_wide", "medication", value.var = 'CMOCCUR', by.col="SUBJID",formula ="SUBJID ~ CMCAT", completeCases = FALSE,
            fun.aggregate="length", datasources = opals)
  
  
  # Glucose lowering drugs
  #Some magic to go from number of lines to 0/1
  ds.asNumeric("medication_wide$CMCAT.BLOOD.GLUCOSE.LOWERING.DRUGS..EXCL..INSULINS", "glucoselowering", datasources=opals)
  dssCut("glucoselowering", in.place=TRUE, breaks=c(-1,0,10),  datasources=opals)
  ds.recodeLevels("glucoselowering", newCategories=c(0,1), newobj = "glucoselowering", datasources = opals)
  ds.asCharacter("glucoselowering", "glucoselowering", datasources=opals)
  ds.asNumeric("glucoselowering", "glucoselowering", datasources=opals)
  
  # Glucose lowering drugs
  #Some magic to go from number of lines to 0/1
  ds.asNumeric("medication_wide$CMCAT.INSULINS.AND.ANALOGUES", "insulin", datasources=opals)
  dssCut("insulin", in.place=TRUE, breaks=c(-1,0,10),  datasources=opals)
  ds.recodeLevels("insulin", newCategories=c(0,1), newobj = "insulin", datasources = opals)
  ds.asCharacter("insulin", "insulin", datasources=opals)
  ds.asNumeric("insulin", "insulin", datasources=opals)  
  ds.dataFrame(x = c('medication_wide$SUBJID','glucoselowering','insulin'), newobj = 'medication', datasources = opals)
  dssSubset(symbol = 'omic.date', what = omics, row.filter = paste0('which(',omics,'$SUBJID %in% medication$SUBJID)'), datasources = opals)
  dssSubset(symbol = 'omic.date', what = 'omic.date', row.filter = 'match(medication$SUBJID,omic.date$SUBJID)', datasources = opals)
  ds.dataFrame(x = c('medication$SUBJID','omic.date$LBDTC','glucoselowering','insulin'), newobj = 'medication', datasources = opals)
  
  # Zero table
  dssSubset(symbol = "no.med", what = omics, row.filter = paste0('-which(',omics,'$SUBJID %in% medication$SUBJID)'), datasources=opals)
  dssSubset(symbol = "no.med", what = "no.med", col.filter = 'c("SUBJID","LBDTC")', datasources = opals)
  ds.vectorCalc(x=c("no.med$SUBJID","no.med$SUBJID"),calc="-",newobj="glucoselowering",datasources=opals)
  ds.vectorCalc(x=c("no.med$SUBJID","no.med$SUBJID"),calc="-",newobj="insulin",datasources=opals)
  ds.cbind(x = c("no.med","glucoselowering","insulin"), newobj = "no.med",datasources = opals)
  dssRbind('MEDICATION', x = 'medication', y = 'no.med', datasources=opals)
  dssSubset(symbol = 'MEDICATION', what = 'MEDICATION', row.filter = paste0('match(',omics,'$SUBJID, MEDICATION$SUBJID))'), datasources = opals)
  ds.assign(toAssign = "MEDICATION$LBDTC", newobj = 'CMSTDTC', datasources = opal)
  ds.dataFrame(x = c('MEDICATION$SUBJID','CMSTDTC','MEDICATION$glucoselowering','MEDICATION$insulin'), newobj = 'MEDICATION', datasources = opals)
}


medication.module.dcs <- function(object=NULL, opal)
{
  if(is.null(object))
  {
    stop("Please provide a object name")
  }

  # Rowfilter
  ds.assign(toAssign = object, newobj = "whiteobj", datasources=opal)
  rowfilter.cm <- paste0("paste0(cm$SUBJID, '.', cm$CMSTDTC) %in% paste0(",
    object,"$SUBJID,'.',", object, "$LBDTC)")

  dssSubset(symbol = "medication",  what = "cm", row.filter =  rowfilter.cm, datasources=opal)
  
  # Transform
  dssPivot("medication_wide", "medication", value.var = 'CMOCCUR', by.col="SUBJID",formula ="SUBJID ~ CMCAT", completeCases = FALSE,
            fun.aggregate="length", datasources = opal)
  # Glucose lowering drugs
  #Some magic to go from number of lines to 0/1
  ds.asNumeric("medication_wide$CMCAT.BLOOD.GLUCOSE.LOWERING.DRUGS..EXCL..INSULINS", "glucoselowering", datasources=opal)
  dssCut("glucoselowering", in.place=TRUE, breaks=c(-1,0,10),  datasources=opal)
  ds.recodeLevels("glucoselowering", newCategories=c(0,1), newobj = "glucoselowering", datasources = opal)
  ds.asCharacter("glucoselowering", "glucoselowering", datasources=opal)
  ds.asNumeric("glucoselowering", "glucoselowering", datasources=opal)

  # Glucose lowering drugs
  #Some magic to go from number of lines to 0/1
  ds.asNumeric("medication_wide$CMCAT.INSULINS.AND.ANALOGUES", "insulin", datasources=opal)
  dssCut("insulin", in.place=TRUE, breaks=c(-1,0,10),  datasources=opal)
  ds.recodeLevels("insulin", newCategories=c(0,1), newobj = "insulin", datasources = opal)
  ds.asCharacter("insulin", "insulin", datasources=opal)
  ds.asNumeric("insulin", "insulin", datasources=opal)
  ds.assign(toAssign = "whiteobj$LBDTC", newobj = "CMSTDTC", datasources=opal)
  ds.dataFrame(x = c('medication_wide$SUBJID',"CMSTDTC",'glucoselowering','insulin'),   newobj = 'MEDICATION', datasources = opal)
}


transformations <- function(Attribute, transform, opal){
	if(transform == 1){
		datashield.assign(opal, Attribute, as.symbol(paste0("log(",Attribute,"+1, 10)")))
	}
	else if(transform == 2){
		dssScale(symbol = Attribute, what = Attribute, datasources = opal,type="split")
	}
	else if(transform == 3){		
		datashield.assign(opal, Attribute, as.symbol(paste0("log(",Attribute,"+1, 10)")))
		dssScale(symbol = Attribute, what = Attribute, datasources = opal, type="split")
	}
}

calc.age <- function(omicAttribute, opal){	
	dssSubset(symbol = 'AGE', what = 'dm', row.filter = paste0('match(',omicAttribute,'$SUBJID, dm$SUBJID)'), datasources = opal)
	ds.assign(toAssign = "AGE$BRTHDTC", newobj = 'BIRTH', datasources = opal)
	ds.assign(toAssign = paste0(omicAttribute,'$LBDTC'), newobj = 'OMICS.DATE', datasources = opal)
	ds.dataFrame(x = c("BIRTH","OMICS.DATE"), newobj = 'AGE', datasources = opal)
	dssAddDaysInterval(newobj = 'AGE',df = 'AGE',description.list =  list(AGE = list(start_date = 'BIRTH', end_date = 'OMICS.DATE')), datasources = opal)
  	ds.assign(toAssign = as.symbol("AGE$AGE / 365"), newobj = 'age', datasources = opal)
  	ds.assign(toAssign = paste0(omicAttribute,'$LBDTC'), newobj = 'DMDTC', datasources = opal) 
	ds.dataFrame(x = c(paste0(omicAttribute,'$SUBJID'),"DMDTC","age"), newobj = 'AGE', datasources = opal)	
}

get.sex <- function(omicAttribute,opal){
	dssSubset(symbol = 'SEX', what = 'dm', row.filter = paste0('match(',omicAttribute,'$SUBJID, dm$SUBJID)'), datasources = opal)
	ds.assign(toAssign = "SEX$SEX", newobj = 'sex', datasources = opal)
	ds.assign(toAssign = paste0(omicAttribute,'$LBDTC'), newobj = 'DMDTC', datasources = opal)
	ds.dataFrame(x = c(paste0(omicAttribute,'$SUBJID'),"DMDTC","sex"), newobj = 'SEX', datasources = opal)
}

godarts.insulinstatus <- function(omicAttribute, opal){
	  dssSubset(symbol = "TIMEINSULINSTATUS", what = 'lb', row.filter = sprintf('%s == "%s"',"LBTESTCD","TIMEINSULINSTATUS"), datasources = opal)      
    dssPivot("TIMEINSULINSTATUS", "TIMEINSULINSTATUS", by.col="SUBJID",value.var = "LBORRES", formula =sprintf("SUBJID + %s  ~ %s","LBDTC","LBTESTCD"),
                completeCases = FALSE,fun.aggregate="mean", datasources = opal)
    #PATIENTS WITH INSULIN STATUS DEFINED
    dssSubset(symbol = 'TIMEINSULINSTATUS.with', what = 'TIMEINSULINSTATUS', row.filter = paste0('na.omit(match(',omicAttribute,'$SUBJID, TIMEINSULINSTATUS$SUBJID))'), datasources = opal)
    dssSubset(symbol = "TIMEINSULINSTATUS.with", what = "TIMEINSULINSTATUS.with", col.filter = 'c("SUBJID","LBDTC")', datasources = opal)
    ds.vectorCalc(x = c('TIMEINSULINSTATUS.with$SUBJID', 'TIMEINSULINSTATUS.with$SUBJID'), calc = '/', newobj = 'LBTESTCD.TIMEINSULINSTATUS', datasources = opal)
    ds.cbind(x = c("TIMEINSULINSTATUS.with","LBTESTCD.TIMEINSULINSTATUS"), newobj = "TIMEINSULINSTATUS.with",datasources = opal)    
    #PARIENTS WITHOUT INSULIN STATUS DEFINED    
    dssSubset(symbol = "TIMEINSULINSTATUS.without", what = omicAttribute, row.filter = paste0('-which(',omicAttribute,'$SUBJID %in% TIMEINSULINSTATUS.with$SUBJID)'), datasources = opal)
    dssSubset(symbol = "TIMEINSULINSTATUS.without", what = "TIMEINSULINSTATUS.without", col.filter = 'c("SUBJID","LBDTC")', datasources = opal)
    ds.vectorCalc(x = c('TIMEINSULINSTATUS.without$SUBJID','TIMEINSULINSTATUS.without$SUBJID'), calc = "-", newobj = 'LBTESTCD.TIMEINSULINSTATUS', datasources = opal)
    ds.cbind(x = c("TIMEINSULINSTATUS.without","LBTESTCD.TIMEINSULINSTATUS"), newobj = "TIMEINSULINSTATUS.without",datasources = opal)
    dssRbind(newobj = 'TIMEINSULINSTATUS', x = 'TIMEINSULINSTATUS.with', y='TIMEINSULINSTATUS.without',datasources = opal)
    dssSubset(symbol = 'TIMEINSULINSTATUS', what = 'TIMEINSULINSTATUS', row.filter = paste0('match(',omicAttribute,'$SUBJID, TIMEINSULINSTATUS$SUBJID))'), datasources = opal)
}



loadTables <- function(attributes,opal){
  tablesList <- list(
    "lb" = "rhapsody.LB",
    "vs" = "rhapsody.VS",
    "dm" = "rhapsody.DM",
    "cm" = "rhapsody.CM",
    "gd" = "rhapsody.GD")
  tables <- unique(unlist(lapply(attributes,function(x){dataStructure[[x]][1]})))
  
  currentObjects <- datashield.symbols(conns = opal)[[1]]#$objects.found
  
  #currentObjects <- ds.ls2(opal)[[names(opal)[1]]]
  for(index in seq_along(tables)){
    table <- tables[index]
    if(table %in% currentObjects){
      message(paste0("\nThe ",table," is already available"))
    }else{
      message(paste0("\nLoading the ", table, " table into environment"))
      datashield.assign(opal,table,tablesList[[table]])
    }
  }
}

prepareData <- function(assign, opal, attributes, removeNA=FALSE, transformVector,cohort, plusOne=FALSE){
  loadTables(attributes,opal)
  
  
  #############################################################
  #SUBSET ALL SELECTED DATAFRAMESss
  #############################################################
  for(index in seq_along(attributes)){
    attr <- attributes[index]    		
    cat(attr)
    arg <- dataStructure[[attributes[index]]]
    formulaValues <- dataframeSelector(arg[1])		
    #When data is one level deep
    #attr <- ifelse(length(dataStructure[[attr]]) == 2,paste0(attr,"DF") ,attr)
    omicAttribute <- attributes[1]
    #DIABETES DURATION
    if(attr == "DIABDUR"){    	
    	diabetes.duration(omicAttribute,opal)    	
    #MEDICATION DATA			 	   	
    }else if (attr == "MEDICATION" & cohort == "dcs"){
    	medication.module.dcs(omicAttribute,opal) 
    }else if(attr == "MEDICATION" & cohort == "godarts"){
    	medication.module.godarts(omicAttribute, 90, opal)        	   	   
    #TIME TO INSULIN 	 	   	
    }else if (attr == "TIMETOINSULIN"){
    	omic.to.insulin(omicAttribute,opal = opal,cohort)    	   	
    }else if(attr == "AGE"){
    	calc.age(omicAttribute,opal)
    }else if(attr == "SEX"){
    	get.sex(omicAttribute,opal)
    #}else if(attr == "TIMEINSULINSTATUS" & cohort == "godarts"){
  #  	godarts.insulinstatus(omicAttribute,opal)
    }else if(attr == "METABOLOMICS" & cohort == "godarts"){
      METABOLS <- "c('AOHB','AADA','SDMA.ADMA','Ala','AzeA','BOHB','CA','CDCA','Cit','Crea','DCA','GBB','GCA','GDCA','GCDCA','Gln','Glu','Gly','GUDCA','Hcit','Ile','IndS','Kynu','Leu','NMNA','Phe','Taurine','TCA','TDCA','TCDCA','Trp','TUDCA','Tyr','UDCA')"
      dssSubset(symbol = 'METABOLOMICS', what = "lb", row.filter = sprintf('%s %%in%% %s',"LBTESTCD",METABOLS), datasources = opal)
      dssPivot('METABOLOMICS', 'METABOLOMICS', by.col="SUBJID",value.var = 'LBORRES', formula =sprintf("SUBJID + %s  ~ %s",'LBDTC','LBTESTCD'),
                completeCases = FALSE,fun.aggregate="mean", datasources = opal)
    }
    ###EXAMPLE###
    #dataStructure[['SEX']] <- c("dm","SEX")
    else if(length(arg) == 2){      
      dssSubset(symbol = attr, what = arg[1], col.filter = sprintf('c("SUBJID", "DMDTC", "%s")',arg[2]) , datasources = opal)
    ###EXAMPLE###
    #dataStructure[['LIPIDOMICS']] <- c('lb','LBMETHOD','LIPIDOMICS')
    }
    else if(length(arg) == 3){
      #Subset data
      dssSubset(symbol = attr, what = arg[1], row.filter = sprintf('%s == "%s"',arg[2],arg[3]), datasources = opal)
      #Transpose
      dssPivot(attr, attr, value.var = formulaValues[2], by.col="SUBJID", formula =sprintf("SUBJID + %s  ~ %s",formulaValues[1],formulaValues[3]),
                completeCases = FALSE, fun.aggregate="mean", datasources = opal)
      
    #When data is two levels deep
    ###EXAMPLE###
    #dataStructure[['CREAT']]<- c('lb','LBMETHOD','CLINICAL','LBTESTCD','CREAT')
    }
    else if(length(arg) == 5){
      #Subset first level
      dssSubset(symbol = attr, what = arg[1], row.filter = sprintf('%s == "%s"',arg[2],arg[3]), datasources = opal)
      #Subset second level
      dssSubset(symbol = attr, what = attr, row.filter = sprintf('%s == "%s"',arg[4],arg[5]), datasources = opal)
      #Transpose      
      if(attr == "HBA1C"){
      	#"GODARTS IS REVERSED"
      	LBORRESHB <- "mmol/mol"#ifelse(cohort == "godarts","%","mmol/mol")      	
      	dssSubset(symbol = attr, what = attr, row.filter = paste0('LBORRESU == "',LBORRESHB,'"'), datasources = opal)
      }
      dssPivot(attr, attr, by.col="SUBJID", value.var = formulaValues[2], formula =sprintf("SUBJID + %s  ~ %s",formulaValues[1],formulaValues[3]),
                completeCases = FALSE,	fun.aggregate="mean", datasources = opal)
    }  

  }
	
  #############################################################
  #MATCH ALL SELECTED ATTRIBUTES BASED ON SUBJID and closest date to OMICS date
  #############################################################
  message("Matching data...")
  matchData(opal, attributes)
  
  
  message("Removing SUBJIDs and DATES")
  ds.assign(toAssign = paste0(omicAttribute,"$SUBJID"),newobj = 'SUBJID', datasources = opal)
  removeColumns(opal,attributes)
   
  TransformData(attributes,transformVector,plusOne, opal)
  newAttr <-c("SUBJID",attributes)
  ds.cbind(x = newAttr, newobj=assign, datasources = opal)
  
}



#############################################################
#MATCH ALL SELECTED ATTRIBUTES BASED ON SUBJID and closest date to OMICS date
#############################################################
matchData <- function(opal, attributes){
  omicAttribute <- attributes[1]
  remainAttributes <- attributes[-1]
  for(idx in seq_along(remainAttributes)){
    attr <- remainAttributes[idx]
    arg <- dataStructure[[attr]]
    formulaValues <- dataframeSelector(arg[1])
    filterID <- paste0(attr, "$SUBJID %in% ",omicAttribute,"$SUBJID")
    dssSubset(symbol = attr, what = attr, row.filter = filterID, datasources = opal)
    
    
    if(attr == "CLUSTER")
    {
      filterDate <- sprintf('apply(%s,1,function(x){
      if(!x["SUBJID"] %%in%% %s$SUBJID){return(NA)}
      else{
        tmpDF <- %s[%s$SUBJID==x["SUBJID"],]
        rownames(tmpDF)
      }
      
      })',omicAttribute,attr,attr,attr,formulaValues[1],formulaValues[1])
    }else{
      filterDate <- sprintf('apply(%s,1,function(x){
      if(!x["SUBJID"] %%in%% %s$SUBJID){return(NA)}
      else{
        tmpDF <- %s[%s$SUBJID==x["SUBJID"],]
        tmpdate1 <- as.Date(as.character(tmpDF$%s))
        tmpdate2 <- as.Date(as.character(x["LBDTC"]))
        idx <- which.min(abs(tmpdate1 - tmpdate2))[1]
        rownames(tmpDF[idx,])
      }
      
      })',omicAttribute,attr,attr,attr,formulaValues[1],formulaValues[1])
    }
   
    
    
    dssSubset(symbol = attr, what = attr, row.filter = filterDate, datasources = opal)
    
    #Final match
    filterFinal <- sprintf('match(%s$SUBJID, %s$SUBJID)', omicAttribute, attr)
    
    dssSubset(symbol = attr, what = attr, row.filter = filterFinal, datasources = opal)
  }
}



#############################################################
#REMOVE SUBJIDS AND LBDTC COLUMNS
#############################################################
removeColumns <- function(opal, attributes){
  for(it in seq_along(attributes)){
    attr <- attributes[it]
    arg <- dataStructure[[attributes[it]]]
    formulaValues <- dataframeSelector(arg[1])
    colFilter <- paste0('-which(colnames(',attr,') %in% c("SUBJID","',formulaValues[1],'"))')
    dssSubset(symbol = attr, what = attr, col.filter = colFilter, datasources = opal)
  }
}

TransformData <- function(attributes,transformVector, plusOne, opal){
  for(iTr in seq_along(attributes)){
    attr <- attributes[iTr]
    transformation <- transformVector[iTr]
    transformations(Attribute = attr,transformation,plusOne, opal)
  }
}

transformations <- function(Attribute, transform,plusOne, opal){
  if(plusOne){
    form <- as.symbol(paste0("log(",Attribute,"+1, 10)"))
  }else{
    form <- as.symbol(paste0("log(",Attribute,", 10)"))
  }
  
  if(transform == 1){
    opal::datashield.assign(opal, Attribute, form)
  }
  else if(transform == 2){
    dssScale(symbol = Attribute, what = Attribute, datasources = opal, type="split")
  }
  else if(transform == 3){
    datashield.assign(opal, Attribute, form)
    dssScale(symbol = Attribute, what = Attribute, datasources = opal, type="split")
  }
}


logIn <- function(datasource, loginPath, cohort, destfile= getwd()){
  logins <- read.table(loginPath, header=T)
  logins <- logins[logins$server %in% c(cohort),]
  assign(datasource,datashield.login(logins), envir = .GlobalEnv)
  #print(get(datasource))
  newline <- paste(c(as.character(Sys.time()),
                     get(datasource)[[cohort]]$username,
                     get(datasource)[[cohort]]$rid), collapse = ",")
  write(newline,file=paste0(destfile,"/datasources.tmp"),append=TRUE)
}




omic.to.insulin <- function(omicAttribute,opal,cohort){
  if(cohort=="dcs"){
    dssSubset(symbol = 'TIMETOINSULIN', what = 'lb', row.filter = 'LBTESTCD == "TIMEINSULIN"', datasources = opal)
    dssSubset(symbol = 'TIMETOINSULIN', what = 'TIMETOINSULIN', row.filter = paste0('match(',omicAttribute,'$SUBJID, TIMETOINSULIN$SUBJID)'), datasources = opal)
    ds.assign(toAssign = "TIMETOINSULIN$LBDTC", newobj = 'INSULIN.DATE', datasources = opal)
    ds.assign(toAssign = paste0(omicAttribute,'$LBDTC'), newobj = 'OMICS.DATE', datasources = opal)
    ds.assign(toAssign = paste0(omicAttribute,'$SUBJID'), newobj = 'SUBJID', datasources = opal)
    ds.asCharacter("OMICS.DATE","OMICS.DATE", datasources = opal)
    ds.asCharacter("INSULIN.DATE","INSULIN.DATE", datasources = opal)
    ds.dataFrame(x = c("SUBJID","OMICS.DATE","INSULIN.DATE"), newobj = 'DATES', datasources = opal)
    dssAddDaysInterval(newobj = 'DATES',df = 'DATES',description.list =  list(LBTESTCD.TIMEINSULIN.OMICS = list(start_date = 'OMICS.DATE', end_date = 'INSULIN.DATE')), datasources = opal)
    dssSubset(symbol = 'TIMETOINSULIN', what = 'DATES', row.filter = "DATES$LBTESTCD.TIMEINSULIN.OMICS > 0", datasources = opal)
    dssSubset(symbol = omicAttribute, what = omicAttribute, row.filter = paste0('match(TIMETOINSULIN$SUBJID,',omicAttribute,'$SUBJID)'), datasources = opal)
    #Organize data frame
    ds.assign(toAssign = "TIMETOINSULIN$SUBJID", newobj = 'SUBJID', datasources = opal)
    ds.assign(toAssign = "TIMETOINSULIN$OMICS.DATE", newobj = 'LBDTC', datasources = opal)
    ds.assign(toAssign = as.symbol("TIMETOINSULIN$LBTESTCD.TIMEINSULIN.OMICS / 365"), newobj = 'TIMEINSULIN', datasources = opal)
    ds.asCharacter("TIMEINSULIN","TIMEINSULIN", datasources = opal)
    ds.dataFrame(x = c("SUBJID","LBDTC","TIMEINSULIN"), newobj = 'TIMETOINSULIN', datasources = opal)
  }
  else if(cohort == "godarts"){
    dssSubset(symbol = 'TIMETOINSULIN', what = 'lb', row.filter = 'LBTESTCD == "TIMEINSULIN"', datasources = opal)
    dssSubset(symbol = 'TIMETOINSULIN', what = 'TIMETOINSULIN', row.filter = paste0('match(',omicAttribute,'$SUBJID, TIMETOINSULIN$SUBJID)'), datasources = opal)
    ds.assign(toAssign = "TIMETOINSULIN$LBDTC", newobj = 'INSULIN.DATE', datasources = opal)
    ds.assign(toAssign = paste0(omicAttribute,'$LBDTC'), newobj = 'OMICS.DATE', datasources = opal)
    ds.assign(toAssign = paste0(omicAttribute,'$SUBJID'), newobj = 'SUBJID', datasources = opal)
    ds.dataFrame(x = c("SUBJID","OMICS.DATE","INSULIN.DATE"), newobj = 'DATES', datasources = opal)
    dssAddDaysInterval(newobj = 'DATES',df = 'DATES',description.list =  list(LBTESTCD.TIMEINSULIN.OMICS = list(start_date = 'OMICS.DATE', end_date = 'INSULIN.DATE')), datasources = opal)
    dssSubset(symbol = 'TIMETOINSULIN', what = 'DATES', row.filter = "DATES$LBTESTCD.TIMEINSULIN.OMICS > 0", datasources = opal)
    dssSubset(symbol = omicAttribute, what = omicAttribute, row.filter = paste0('match(TIMETOINSULIN$SUBJID,',omicAttribute,'$SUBJID)'), datasources = opal)
    #Organize data frame
    ds.assign(toAssign = "TIMETOINSULIN$SUBJID", newobj = 'SUBJID', datasources = opal)
    ds.assign(toAssign = "TIMETOINSULIN$OMICS.DATE", newobj = 'LBDTC', datasources = opal)
    ds.assign(toAssign = as.symbol("TIMETOINSULIN$LBTESTCD.TIMEINSULIN.OMICS / 365"), newobj = 'TIMEINSULIN', datasources = opal)
    ds.dataFrame(x = c("SUBJID","LBDTC","TIMEINSULIN"), newobj = 'TIMETOINSULIN', datasources = opal)
  }
  else if(cohort == "andis"){
    dssSubset(symbol = 'TIMETOINSULIN', what = 'lb', row.filter = 'LBTESTCD == "TIMEINSULIN"', datasources = opal)
    dssSubset(symbol = 'TIMETOINSULIN', what = 'TIMETOINSULIN', row.filter = paste0('match(',omicAttribute,'$SUBJID, TIMETOINSULIN$SUBJID)'), datasources = opal)
    ds.assign(toAssign = "TIMETOINSULIN$LBDTC", newobj = 'INSULIN.DATE', datasources = opal)
    ds.assign(toAssign = paste0(omicAttribute,'$LBDTC'), newobj = 'OMICS.DATE', datasources = opal)
    ds.assign(toAssign = paste0(omicAttribute,'$SUBJID'), newobj = 'SUBJID', datasources = opal)
    ds.dataFrame(x = c("SUBJID","OMICS.DATE","INSULIN.DATE"), newobj = 'DATES', datasources = opal)
    dssAddDaysInterval(newobj = 'DATES',df = 'DATES',description.list =  list(LBTESTCD.TIMEINSULIN.OMICS = list(start_date = 'OMICS.DATE', end_date = 'INSULIN.DATE')), datasources = opal)
    dssSubset(symbol = 'TIMETOINSULIN', what = 'DATES', row.filter = "DATES$LBTESTCD.TIMEINSULIN.OMICS > 0", datasources = opal)
    dssSubset(symbol = omicAttribute, what = omicAttribute, row.filter = paste0('match(TIMETOINSULIN$SUBJID,',omicAttribute,'$SUBJID)'), datasources = opal)
    #Organize data frame
    ds.assign(toAssign = "TIMETOINSULIN$SUBJID", newobj = 'SUBJID', datasources = opal)
    ds.assign(toAssign = "TIMETOINSULIN$OMICS.DATE", newobj = 'LBDTC', datasources = opal)
    ds.assign(toAssign = as.symbol("TIMETOINSULIN$LBTESTCD.TIMEINSULIN.OMICS / 365"), newobj = 'TIMEINSULIN', datasources = opal)
    ds.dataFrame(x = c("SUBJID","LBDTC","TIMEINSULIN"), newobj = 'TIMETOINSULIN', datasources = opal)
  }
  
}