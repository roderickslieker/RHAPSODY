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
  cox.res <- tryCatch(suppressWarnings(ds2.coxph(formula = form, data = data, async = FALSE, datasources = opal)), error=function(e) return("Error"))

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