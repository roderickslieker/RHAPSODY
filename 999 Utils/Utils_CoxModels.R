#' @param vars These are the attributes that are tested with a cox model
#' @param formula A formula as character
#' @param targets A dataframe with to convert the attributes to full names (available on federated node)
#' @param cohort This to get the right cohort back
#' @param omicType Type of omics that need to be tested
#' @param data Name of the datashield dataframe
CalculateCoxModel <- function(vars, formula, targets = NULL, cohort, omicType, data){
	if(omicType == "MIRNA"){		
		ModelOut <- getMiRNAModelData(vars = vars, formula = formula, targets = targets, cohort = cohort, data = data)
	}
	else if(omicType == "METABOLOMICS"){
		ModelOut <- getMetaboliteModelData(vars = vars, formula = formula, cohort = cohort, data = data)
	}
	else if(omicType == "LIPIDOMICS"){
		ModelOut <- getLipidModelData(vars = vars, formula = formula, targets = targets, cohort = cohort, data = data)
	}
	else if (omicType == "SOMALOGIC"){
		ModelOut <- getPeptideModelData(vars = vars, formula = formula, targets = targets, cohort = cohort, data = data)
	}
	else if (omicType == "GWAS"){
		ModelOut <- getSNPModelData(vars = vars, formula = formula, targets = targets, cohort = cohort, data = data)
	}
	else{
		message("\nOmic type not available for cox models.\nOptions are:\n\tMIRNA\n\tLIPIDOMICS\n\tMETABOLOMICS\n\tSOMALOGIC\n\tGWAS\n")
	}
	return(ModelOut)
}


#' Function to retrieve model from FDB, fucntion expectes data called ModelData on the FDB, with vars time to spin (LBTESTCD.TimeToSpin), LBTEST.TIMEINSULIN and LBTEST.TIMEINSULINSTATUS.
#' @param peptide This is the swiss peptide ID, so for example LBTESTCD.SL000001
#' @param cohort This to get the right cohort back
getCoxModel <- function(peptide, form, cohort="dcs", nameVar = "%s", data){
  cat(peptide,'\n')
  form <- sprintf(form, peptide)
  form <- as.formula(form)

  #Fit GLM
  cox.res <- tryCatch(suppressWarnings(dssCoxph(formula = form, data = data, async = FALSE, datasources = opal)), error=function(e) return("Error"))

  if(cox.res[1] == "Error")
  {
    cat(sprintf("Error model did not work for var %s", peptide))
    peptide2 <- gsub("LBTESTCD\\.","",peptide)
    out <- data.frame(variable=peptide2, hr=NA, se.hr=NA, lower=NA, upper=NA, p.val=NA, n.case=NA, n.event=NA)

  }else{
    #Select cohort
    cox.res.sel <- cox.res[[cohort]]

    #Create out table
    n.case <- cox.res.sel$model$n
    n.event <- cox.res.sel$model$nevent
    sm <- summary(cox.res.sel$model)

    #New log name
    peptide.log <- sprintf(nameVar, peptide)

    #Fetch
    p.val <- sm$coefficients[peptide.log,"Pr(>|z|)"]
    lower <- sm$conf.int[peptide.log,"lower .95"]
    upper <- sm$conf.int[peptide.log,"upper .95"]
    hr <- sm$conf.int[peptide.log,"exp(coef)"]
    se.hr <- sm$coefficients[peptide.log,"se(coef)"]

    peptide2 <- gsub("LBTESTCD\\.","",peptide)
    out <- data.frame(variable=peptide2, hr, se.hr, lower, upper, p.val, n.case, n.event)

    rm(cox.res, form,peptide2, hr, se.hr, lower, upper, p.val, n.case, n.event)

  }

  #out
  return(out)
}

#' Function to run a cox model on the federated node.
#' @param vars These are the lipids that are tested with a cox model on ModelData
#' @param formula A formula as character
#' @param targets A dataframe with to convert the lipids to full names (available on federated node)
#' @param cohort A string defining whether dcs or godarts is called
#' @param data Name of the datashield dataframe
getPeptideModelData <- function(vars, formula, targets,cohort, data){
	if(is.null(targets)){
		stop("Please supply target file")
	}
  	res <- lapply(vars, getCoxModel,form = formula,cohort=cohort, nameVar="%s", data=data)
  	results <- do.call(rbind, res)
  	targets.s <- targets[match(results$variable, targets$SomaId),]
  	results <- data.frame(targets.s[,c(4:6)],results)
  	results.somalogic <- results[order(results$p.val, decreasing=F),]
  	return(results.somalogic)
}

getMiRNAModelData <- function(vars, formula, targets, cohort, data){
	if(is.null(targets)){
		stop("Please supply target file")
	}
	res <- lapply(vars, getCoxModel,form = formula,cohort=cohort, nameVar = "%s", data=data)
  	results <- do.call(rbind, res)
  	results <- results[order(results$p.val, decreasing=F),]
  	return(results)
}


getMetaboliteModelData <- function(vars, formula, cohort, data){
	res <- lapply(vars, getCoxModel,form = formula,cohort=cohort, nameVar="%s", data=data)
	results <- do.call(rbind, res)
	results.metabolites <- results[order(results$p.val, decreasing=F),]
	return(results.metabolites)
}

getLipidModelData <- function(vars, formula, targets, cohort, data){
	if(is.null(targets)){
		stop("Please supply target file")
	}
  	res <- lapply(vars, getCoxModel,form = formula,cohort=cohort, nameVar = "%s", data=data)
  	results <- do.call(rbind, res)
  	targets <- prepare.targets(targets)
  	targets.s <- targets[match(results$variable, targets$id),]
  	results <- data.frame(targets.s[,c(2,3)],results)
  	results.lipidomics <- results[order(results$p.val, decreasing=F),]
  	return(results.lipidomics)
}

getSNPModelData <- function(vars, formula, targets, cohort,data){
	if(is.null(targets)){
		stop("Please supply target file")
	}
  	res <- lapply(vars, getCoxModel,form = formula,cohort=cohort, nameVar = '%s', data=data)
  	results <- do.call(rbind, res)
  	results <- results[order(results$p.val, decreasing=F),]
  	return(results)
}

prepare.targets <- function(targets)
{
  targets$id <- gsub(" ","\\.",targets$feature)
  targets$id <- gsub(";",".",targets$id)
  targets$id <- gsub(":",".",targets$id)
  targets$id <- gsub("/",".",targets$id)
  targets$id <- gsub("-",".",targets$id)

  return(targets)

}


CombineCoxThree <- function(Data01, Data02, Data03, ColVar, studyLabels){
  #Check if dataframes have equal number of rows
  if(!nrow(Data01) == nrow(Data02)){
    stop("Dataframes are not of equal length")
  }
  #Check if variables are equal and in same order
  if(!table(as.character(Data01[,ColVar]) == as.character(Data02[,ColVar]))==nrow(Data01)){
    stop("Variables of dataframes are not in same order")
  }
  if(!table(as.character(Data01[,ColVar]) == as.character(Data03[,ColVar]))==nrow(Data01)){
    stop("Variables of dataframes are not in same order")
  }
  
  combined <- data.frame()
  variables <-as.character(Data01[,ColVar])
  
  for(i in seq_along(variables)){
    #cat(i)
    var <- variables[i]
    tmp.merge <- rbind(Data01[Data01[,ColVar] == var,1:8],
                       Data02[Data02[,ColVar] == var,1:8],
                       Data03[Data03[,ColVar] == var,1:8])
    
    
    tmp.combined <- metagen(TE = log(tmp.merge$hr),
                            seTE = tmp.merge$se.hr,
                            studlab=studyLabels,
                            comb.fixed = TRUE,
                            comb.random = TRUE,
                            sm="HR")
    cis <- ci(TE = tmp.combined$TE.random, seTE = tmp.combined$seTE.random, level=0.95, df=NULL, null.effect = 0)
    summ <- data.frame(var,
                       #Effect.fixed = tmp.combined$TE.fixed,
                       #SE.fixed = tmp.combined$seTE.fixed,
                       #Zval.fixed = tmp.combined$zval.fixed,
                       #Pval.fixed = tmp.combined$pval.fixed,
                       Effect.random = exp(tmp.combined$TE.random),
                       #SE.random = tmp.combined$seTE.random,
                       lower = exp(cis$lower),
                       upper = exp(cis$upper),
                       Zval.random = tmp.combined$zval.random,
                       Pval.random  = tmp.combined$pval.random,
                       I2 = summary(tmp.combined)$I2,
                       Het = pchisq(tmp.combined$Q,1, lower.tail=F))
    
    combined <- rbind(combined,summ)
  }
  #combined <- combined[order(combined$Pval.random, decreasing=F),]
  combined$fdr.random <- p.adjust(combined$Pval.random,method="fdr")
  return(combined)
}




CombineCox <- function(Data01, Data02, ColVar, studyLabels){
  #Check if dataframes have equal number of rows
  if(!nrow(Data01) == nrow(Data02)){
    stop("Dataframes are not of equal length")
  }
  #Check if variables are equal and in same order
  if(!table(as.character(Data01[,ColVar]) == as.character(Data02[,ColVar]))==nrow(Data01)){
    stop("Variables of dataframes are not in same order")
  }
  
  
  combined <- data.frame()
  variables <-as.character(Data01[,ColVar])
  
  for(i in seq_along(variables)){
    #cat(i)
    var <- variables[i]
    tmp.merge <- rbind(Data01[Data01[,ColVar] == var,1:8],
                       Data02[Data02[,ColVar] == var,1:8])
    
    tmp.combined <- metagen(TE = log(tmp.merge$hr),
                            seTE = tmp.merge$se.hr,
                            studlab=studyLabels,
                            comb.fixed = TRUE,
                            comb.random = TRUE,
                            sm="HR")
    cis <- ci(TE = tmp.combined$TE.random, seTE = tmp.combined$seTE.random, level=0.95, df=NULL, null.effect = 0)
    summ <- data.frame(var,
                       #Effect.fixed = tmp.combined$TE.fixed,
                       #SE.fixed = tmp.combined$seTE.fixed,
                       #Zval.fixed = tmp.combined$zval.fixed,
                       #Pval.fixed = tmp.combined$pval.fixed,
                       Effect.random = exp(tmp.combined$TE.random),
                       #SE.random = tmp.combined$seTE.random,
                       lower = exp(cis$lower),
                       upper = exp(cis$upper),
                       Zval.random = tmp.combined$zval.random,
                       Pval.random  = tmp.combined$pval.random,
                       I2 = summary(tmp.combined)$I2,
                       Het = pchisq(tmp.combined$Q,1, lower.tail=F))
    
    combined <- rbind(combined,summ)
  }
  #combined <- combined[order(combined$Pval.random, decreasing=F),]
  combined$fdr.random <- p.adjust(combined$Pval.random,method="fdr")
  return(combined)
}




makeMetaSignMatrix <- function(N1,N2,N3,N4,adjust=NULL)
{
  #Assign names
  for(k in 1:4)
  {
    temp <- get(get(sprintf("N%s", k)))
    assign(x = sprintf("Data%s", k), value = temp)
    rm(temp)
  }
  
  
  if(!is.null(adjust))
  {
    Data1$p.adj <- p.adjust(Data1$Pval.random.adj, method = adjust)
    Data2$p.adj <- p.adjust(Data2$Pval.random.adj, method = adjust)
    Data3$p.adj <- p.adjust(Data3$Pval.random.adj, method = adjust)
    Data4$p.adj <- p.adjust(Data4$Pval.random.adj, method = adjust)
  }
  
  sign.matrix <- matrix(rep(0), ncol=4, nrow=nrow(Data1))
  
  #rename
  rename <- function(x){x <- gsub("meta.","",x); x <- gsub("\\.","/",x)}
  new.colnames <- do.call(c, lapply(c(N1,N2,N3,N4), rename))
  colnames(sign.matrix) <- new.colnames
  rownames(sign.matrix) <- Data1$var
  
  Data1.sign <- Data1[Data1$p.adj <= 0.05,]
  Data2.sign <- Data2[Data2$p.adj <= 0.05,]
  Data3.sign <- Data3[Data3$p.adj <= 0.05,]
  Data4.sign <- Data4[Data4$p.adj <= 0.05,]
  
  nrow(Data1.sign)
  nrow(Data2.sign)
  nrow(Data3.sign)
  nrow(Data4.sign)
  
  sign.matrix[match(Data1.sign$var, rownames(sign.matrix)), new.colnames[1]] <- 1
  sign.matrix[match(Data2.sign$var, rownames(sign.matrix)), new.colnames[2]] <- 1
  sign.matrix[match(Data3.sign$var, rownames(sign.matrix)), new.colnames[3]] <- 1
  sign.matrix[match(Data4.sign$var, rownames(sign.matrix)), new.colnames[4]] <- 1
  return(sign.matrix)
}
