#' Get meta-analyzed data for two cluster groups
#' @param group1 This is the name of the first cluster
#' @param group2 This is the name of the second cluster
#' @export
getMetaData <- function(group1, group2)
{
  data1 <- get(sprintf("results.%s.%s", group1, group2))
  data2 <- get(sprintf("results.%s.%s.GoDARTS", group1, group2))

  vars <- unique(data1$var)

  .getForVar <- function(var.sel)
  {
    #k=1
    var.sel <- as.character(var.sel)

    #Subset
    dx <- rbind(data1[match(var.sel, data1$var),],
                data2[match(var.sel, data2$var),])

    #Combine
    fit <- metagen(TE = dx$Estimate, seTE = dx$Std..Error, studlab = dx$study, comb.fixed = TRUE, comb.random = TRUE)
    out <- data.frame(var.sel,
                      Effect.fixed = fit$TE.fixed, SE.fixed = fit$seTE.fixed,
                      Zval.fixed = fit$zval.fixed, Pval.fixed = fit$pval.fixed,
                      Effect.random = fit$TE.random, SE.random = fit$seTE.random,
                      Zval.random = fit$zval.random, Pval.random  = fit$pval.random,
                      I2 = summary(fit)$I2$TE, Het = pchisq(fit$Q,1, lower.tail=F))


    out
  }
  res.all <- lapply(vars, .getForVar)
  outData <- do.call(rbind, res.all)
  #Sort pval
  outData <- outData[order(outData$Pval.random, decreasing=F),]

  #ADj pval
  outData$Pval.random.adj <- p.adjust(outData$Pval.random, 'bonferroni')
  print(table(outData$Pval.random.adj <= 0.05))

  #Return to disk and envir
  write.table(outData[outData$Pval.random.adj <= 0.05,], file=sprintf('./003.Figures.combined/Lipidomics meta %s %s.csv', group1, group2), sep=",", quote=F, row.names=F)
  return(outData)
}

#' Specific function for lipid data to make a plot of the mean values of the clusters
#' @param var Lipid to be plotted, requires mean.LipidData.DCS and meanLipidData.godarts
#' @export
makeMetaPlot <- function(var)
{
  var2 <- paste0("LBTESTCD.",var)
  d1 <- rbind(mean.LipidData.DCS[mean.LipidData.DCS$lip %in% var2,],
              mean.LipidData.godarts[mean.LipidData.godarts$lip %in% var2,])

  colnames(d1)

  d1m <- melt(d1, id.vars = c('lip','Mean','cluster','Cohort'))
  d1m$cluster <- ordered(d1m$cluster,
                         levels = c(2,3,4,5,6),
                         labels = c("2/SIDD", "3/SIRD", "4/MOD","5/MARD","6/MARD-HDL"))
  d1m$cluster2 <- factor(as.character(d1m$cluster))

  p1 <- ggplot(d1m, aes(x=cluster2, y=value, col=Cohort, fill=cluster2))+
    geom_boxplot()+
    #facet_grid(~Cohort)+
    scale_fill_manual(values = values)+
    ggtitle(var)+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
    ylab("Levels (log10 pmol)")+
    xlab("Cluster")+
    scale_y_continuous(trans="log10")+
    scale_color_manual(values = c("#333333","#333333"))+
    theme(legend.position = "none")
  return(p1)
}

get.plot.hr <- function(dataDCS, dataGODARTS, title)
{
  pd <- data.frame(var = dataDCS[,1], Beta.DCS = dataDCS$hr, Beta.GoDARTS = dataGODARTS$hr)
  pd$Sign <- ifelse(dataDCS$p.adj <= 0.05 & dataGODARTS$p.adj <= 0.05, "Sign", "NonSign")
  #cov <- conversion[match(dataDCS$var, conversion$SomaId),]

  pd$label <- ifelse(pd$Sign == "Sign", cov$Target, NA)

  minx <- round(min(c(pd$Beta.DCS, pd$Beta.GoDARTS)), 1)
  maxx <- round(max(c(pd$Beta.DCS, pd$Beta.GoDARTS)), 1)



  P <- ggplot(pd, aes(x=Beta.DCS, y=Beta.GoDARTS, label=label, col=Sign))+
    geom_point()+
    geom_text()+
    scale_colour_manual(values = c("#A4A4A4","#21908CFF"))+
    xlab("Estimate DCS")+
    ylab("Estimate GoDARTS")+
    ggtitle(title)+
    geom_smooth(method=lm, col='black')+
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0)+
    theme(legend.position = 'none')+
    scale_x_continuous(trans='log2')+
    scale_y_continuous(trans='log2')
  #ymin(minx,maxx)+
  #xlim(minx,maxx)
  return(P)
}


getMetaDataCox <- function(Model)
{
  data1 <- get(sprintf("%s.Data.DCS", Model))
  data2 <- get(sprintf("%s.Data.GoDARTS", Model))

  #add study
  data1$study <- "DCS"
  data2$study <- "GoDARTS"

  vars <- unique(data1$var)

  .getForVar <- function(var.sel)
  {
    #k=1
    var.sel <- as.character(var.sel)

    #Subset
    dx <- rbind(data1[match(var.sel, data1$var),],
                data2[match(var.sel, data2$var),])

    #Combine
    fit <- metagen(TE = log(dx$hr), seTE = dx$se.hr, studlab = dx$study, comb.fixed = TRUE, comb.random = TRUE, sm="HR")
    out <- data.frame(var.sel,
                      Effect.fixed = fit$TE.fixed, SE.fixed = fit$seTE.fixed,
                      Zval.fixed = fit$zval.fixed, Pval.fixed = fit$pval.fixed,
                      Effect.random = fit$TE.random, SE.random = fit$seTE.random,
                      Zval.random = fit$zval.random, Pval.random  = fit$pval.random,
                      I2 = summary(fit)$I2$TE, Het = pchisq(fit$Q,1, lower.tail=F))


    out
  }
  res.all <- lapply(vars, .getForVar)
  outData <- do.call(rbind, res.all)
  #Sort pval
  outData <- outData[order(outData$Pval.random, decreasing=F),]

  #ADj pval
  outData$Pval.random.adj <- p.adjust(outData$Pval.random, 'fdr')
  print(table(outData$Pval.random.adj <= 0.05))
  return(outData)
}







#' Function to retrieve model from FDB, fucntion expectes data called ModelData on the FDB, with vars time to spin (LBTESTCD.TimeToSpin), LBTEST.TIMEINSULIN and LBTEST.TIMEINSULINSTATUS.
#' @param peptide This is the swiss peptide ID, so for example LBTESTCD.SL000001
#' @param cohort This to get the right cohort back
getCoxModel <- function(peptide, form, cohort="dcs", nameVar = "scale(log10(%s))")
{
  cat(peptide,'\n')
  form <- sprintf(form, peptide)
  form <- as.formula(form)

  #Fit GLM
  cox.res <- tryCatch(suppressWarnings(ds2.coxph(formula = form, data = 'ModelData', async = FALSE, datasources = opal)), error=function(e) return("Error"))

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



prepare.targets <- function(targets)
{
  targets$id <- gsub(" ","\\.",targets$feature)
  targets$id <- gsub(";",".",targets$id)
  targets$id <- gsub(":",".",targets$id)
  targets$id <- gsub("/",".",targets$id)
  targets$id <- gsub("-",".",targets$id)

  return(targets)

}

.getModel <- function(var.sel, formula)
{
  cat(var.sel, "\n")
  formula2 <- sprintf(formula, var.sel)
  
  #Fit GLM
  ####################################################################################
  fit <- suppressMessages(tryCatch(ds.glm(formula = formula2, family = 'gaussian', datasources = opal, viewIter=FALSE), error=function(e) return("Error")))
  ####################################################################################
  
  if(fit[1] == "Error")
  {
    cat(sprintf("Error model did not work for lipid %s", var.sel))
    var.sel.normal <- gsub("LBTESTCD.","",var.sel)
    out <- data.frame(var = var.sel.normal, Estimate = NA, 'Std..Error' = NA, z.value = NA, p.value = NA, low0.95CI=NA, high0.95CI=NA)
    
  }else{
    var.sel.normal <- gsub("LBTESTCD.","",var.sel)
    out <- data.frame(var = var.sel.normal, t(fit$coefficients['ModelData.subset$newCluster',]))
    
  }
  
  
  #Push out
  return(out)
}


#' Function to retrieve linear model from FDB
#' @param vars Lipids or peptids
#' @param cohort This to get the right cohort back
getModelForClusters <- function(vars, formula = 'ModelData.subset$%s~ModelData.subset$newCluster', compCluster, refCluster=6, family='gaussian')
{
  #Subset data for clusters
  rowFilter <- paste0("ModelData$LBTESTCD.CLUSTER %in% c(", compCluster,",",refCluster, ")")
  dssSubset(symbol = 'ModelData.subset', what = 'ModelData', row.filter = rowFilter, datasources = opal)
  
  ds.asFactor('ModelData.subset$LBTESTCD.CLUSTER', newobj = 'CLUSTER.SUBSET',datasources = opal)
  lvls <- ds.levels("CLUSTER.SUBSET", datasources = opal)[[1]]
  lvls <- lvls$Levels
  
  #check for number of lines input and recode if necessary
  if(length(lvls) > 2)
  {
    lvls <- ifelse(lvls == compCluster,"0","1")
    ds.recodeLevels(x = "CLUSTER.SUBSET", newCategories = lvls, newobj = "CLUSTER.SUBSET.RECODE",datasources = opal)
  }else{
    ds.recodeLevels(x = 'CLUSTER.SUBSET', newCategories= c(1,0), newobj = 'CLUSTER.SUBSET.RECODE',datasources = opal)
  }
  
  
  ds.asNumeric(x = 'CLUSTER.SUBSET.RECODE', newobj = 'newCluster', datasources = opal)
  ds.dataFrame(x = 'newCluster', datasources = opal)
  ds.cbind(x = c('ModelData.subset','newCluster'), newobj = 'ModelData.subset', datasources = opal)
  
  #Assign variables (much faster than first assigning to DF)
  #ds.asFactor('ModelData.subset$SEX', newobj = 'SEX')
  #ds.assign('ModelData.subset$LBTESTCD.TimeToSpin', newobj = 'TimeToSpin')
  
  #Run model
  res <- lapply(vars, .getModel, formula=formula)
  rmx <- which(sapply(1:length(res), FUN=function(x){sum(is.na(res[[x]]))})!=0)
  if(length(rmx) !=0)
  {
    res <- res[-rmx]
  }
  results <- do.call(rbind, res)
  
  results$comparison <- sprintf('C%s/C%s', refCluster, compCluster)
  results$p.adj <- p.adjust(results$p.value, method = "bonferroni")
  return(results)
}









makeSignMatrix <- function(Data1,Data2,Data3,Data4,adjust=NULL)
{
  if(!is.null(adjust))
  {
    Data1$p.adj <- p.adjust(Data1$p.adj, method = adjust)
    Data2$p.adj <- p.adjust(Data2$p.adj, method = adjust)
    Data3$p.adj <- p.adjust(Data3$p.adj, method = adjust)
    Data4$p.adj <- p.adjust(Data4$p.adj, method = adjust)
  }
  sign.matrix <- matrix(rep(0), ncol=4, nrow=nrow(Data1))
  colnames(sign.matrix) <- c("SIDD","SIRD","MOD","MARD")
  rownames(sign.matrix) <- Data1$var

  Data1.sign <- Data1[Data1$p.adj <= 0.05,]
  Data2.sign <- Data2[Data2$p.adj <= 0.05,]
  Data3.sign <- Data3[Data3$p.adj <= 0.05,]
  Data4.sign <- Data4[Data4$p.adj <= 0.05,]

  nrow(Data1.sign)
  nrow(Data2.sign)
  nrow(Data3.sign)
  nrow(Data4.sign)

  sign.matrix[match(Data1.sign$var, rownames(sign.matrix)), "SIDD"] <- 1
  sign.matrix[match(Data2.sign$var, rownames(sign.matrix)), "SIRD"] <- 1
  sign.matrix[match(Data3.sign$var, rownames(sign.matrix)), "MOD"] <- 1
  sign.matrix[match(Data4.sign$var, rownames(sign.matrix)), "MARD"] <- 1
  return(sign.matrix)
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

getEnrichment <- function(signData, lipidClassData,name)
{
  #n00
  nxx <- table(lipidClassData$Class)
  lcd <- lipidClassData[match(signData$var.sel, lipidClassData$FDBName),"Class"]
  n11 <- table(lcd)
  nss <- nrow(signData)

  sign.table <- data.frame(class = names(n11), n11 = as.numeric(n11), n01=as.numeric(nss-n11), n10=as.numeric(nxx-n11))
  sign.table$n00 <- nrow(lipidClassData) - sign.table$n10 - sign.table$n01 - sign.table$n11
  sign.table$perc <- sign.table$n11 / nss

  getOR <- function(x, sign.table)
  {
    fit <- oddsratio(sign.table$n11[x],sign.table$n01[x],sign.table$n10[x],sign.table$n00[x])
    out <- data.frame(var = as.character(sign.table$class[x]), OR = fit$estimate, CI.lower = fit$conf.int[[1]], CI.upper = fit$conf.int[[2]], P.value = fit$p.value)
    return(out)
  }

  add.OR <- do.call(rbind, lapply(1:nrow(sign.table), getOR, sign.table=sign.table))


  return.out <- data.frame(sign.table, add.OR[,-1])
  return.out$name <- name
  return(return.out)
}


getTable <- function(direction, SIDD, SIRD, MOD, MARD)
{
  SIDD <- SIDD[SIDD$Pval.random.adj <= 0.05,]
  SIRD <- SIRD[SIRD$Pval.random.adj <= 0.05,]
  MOD <- MOD[MOD$Pval.random.adj <= 0.05,]
  MARD <- MARD[MARD$Pval.random.adj <= 0.05,]

  if(direction == "negative")
  {
    SIDD.n <- SIDD[SIDD$Effect.random <= 0 ,]
    SIRD.n <- SIRD[SIRD$Effect.random <= 0,]
    MOD.n  <- MOD[MOD$Effect.random <= 0,]
    MARD.n <- MARD[MARD$Effect.random <= 0,]

    enrichments <- rbind(
      getEnrichment(signData=SIDD.n, lipidClassData=lipidClass.subset, name="SIDD"),
      getEnrichment(signData=SIRD.n, lipidClassData=lipidClass.subset, name="SIRD"),
      getEnrichment(signData=MOD.n, lipidClassData=lipidClass.subset, name="MOD"),
      getEnrichment(signData=MARD.n, lipidClassData=lipidClass.subset, name="MARD"))
    enrichments$name <- factor(enrichments$name, levels=c("SIDD","SIRD","MOD","MARD"))
    enrichments$direction <- direction
  }else{
    SIDD.p <- SIDD[SIDD$Effect.random >= 0,]
    SIRD.p <- SIRD[SIRD$Effect.random >= 0,]
    MOD.p  <- MOD[MOD$Effect.random >= 0,]
    MARD.p <- MARD[MARD$Effect.random >= 0,]

    enrichments <- rbind(
      getEnrichment(signData=SIDD.p, lipidClassData=lipidClass.subset, name="SIDD"),
      getEnrichment(signData=SIRD.p, lipidClassData=lipidClass.subset, name="SIRD"),
      getEnrichment(signData=MOD.p, lipidClassData=lipidClass.subset, name="MOD"),
      getEnrichment(signData=MARD.p, lipidClassData=lipidClass.subset, name="MARD"))
    enrichments$name <- factor(enrichments$name, levels=c("SIDD","SIRD","MOD","MARD"))
    enrichments$direction <- direction
  }


  enrichments$ymin <- ifelse(enrichments$OR <= 1, enrichments$CI.lower, enrichments$OR)
  enrichments$ymax <- ifelse(enrichments$OR <= 1, enrichments$OR, enrichments$CI.upper)
  return(enrichments)
}




getMiRNAModelData <- function(vars, formula, targets, cohort)
{
  res <- lapply(vars, getCoxModel,
                form = formula,
                cohort=cohort, nameVar = "scale(log10(%s + 1))")
  results <- do.call(rbind, res)
  results <- results[order(results$p.val, decreasing=F),]
  return(results)
}

#'
getMetaboliteModelData <- function(vars, formula, cohort)
{
  res <- lapply(vars, getCoxModel,
                form = formula,
                cohort=cohort, nameVar="%s")
  results <- do.call(rbind, res)
  #targets <- prepare.targets(targets)

  #targets.s <- targets[match(results$variable, targets$id),]
  #results <- data.frame(targets.s[,c(2,3)],results)
  results.metabolites <- results[order(results$p.val, decreasing=F),]
  return(results.metabolites)
}

#' Function to run a cox model on the lipid data on the federated node.
#' @param vars These are the lipids that are tested with a cox model on ModelData
#' @param formula A formula as character
#' @param targets A dataframe with to convert the lipids to full names (available on federated node)
#' @param cohort A string defining whether dcs or godarts is called
getLipidModelData <- function(vars, formula, targets, cohort)
{
  res <- lapply(vars, getCoxModel,
                form = formula,
                cohort=cohort, nameVar = "%s")
  results <- do.call(rbind, res)
  targets <- prepare.targets(targets)

  targets.s <- targets[match(results$variable, targets$id),]
  results <- data.frame(targets.s[,c(2,3)],results)
  results.lipidomics <- results[order(results$p.val, decreasing=F),]
  return(results.lipidomics)
}


#' Function to run a cox model on the lipid data on the federated node.
#' @param vars These are the lipids that are tested with a cox model on ModelData
#' @param formula A formula as character
#' @param targets A dataframe with to convert the lipids to full names (available on federated node)
#' @param cohort A string defining whether dcs or godarts is called
getSNPModelData <- function(vars, formula, targets, cohort)
{
  res <- lapply(vars, getCoxModel,
                form = formula,
                cohort=cohort, nameVar = '%s')
  results <- do.call(rbind, res)
  results <- results[order(results$p.val, decreasing=F),]
  return(results)
}


#' Function to run a cox model on the federated node.
#' @param vars These are the lipids that are tested with a cox model on ModelData
#' @param formula A formula as character
#' @param targets A dataframe with to convert the lipids to full names (available on federated node)
#' @param cohort A string defining whether dcs or godarts is called
getPeptideModelData <- function(vars, formula, targets,cohort, nameVar = "scale(log10(%s))")
{
  res <- lapply(vars, getCoxModel,
                form = formula,
                cohort=cohort, nameVar=nameVar)
  results <- do.call(rbind, res)

  targets.s <- targets[match(results$variable, targets$SomaId),]
  results <- data.frame(targets.s[,c(4:6)],results)
  results.somalogic <- results[order(results$p.val, decreasing=F),]

  return(results.somalogic)
}


getAndSavePlot <- function(var, meanData)
{
  #meanData = mean.Metabolomics
  #var = "Kynu"
  p1 <- ggplot(meanData[meanData$FDBName %in% var,], aes(x=factor(cluster), y=value, fill=factor(cluster)))+
    geom_boxplot()+
    scale_fill_manual(values = values)+
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), legend.position = 'none')+
    xlab("Cluster")+
    scale_x_discrete(labels=c("2/SIDD","3/SIRD","4/MOD","5/MARD","6/MARD high HDL"))+
    ylab("Levels")

  pdf(sprintf("./003.Figures/Figure of %s.pdf",var), width=2.5, height=3)
  print(p1)
  dev.off()
  return(p1)
}

generateDSTable <- function()
{
  #Sex 
  SEX <- ds.table1D('ModelData$SEX', datasources = opal, type = "split")
  
  .getVar <- function(var, newVar)
  {
    AGE.dcs <- ds.mean(paste0('ModelData$',var), datasources = opal["dcs"], type = "split")
    AGE.godarts <- ds.mean(paste0('ModelData$',var), datasources = opal["godarts"], type = "split")
    out <- data.frame(Variable = newVar, DCS = round(as.numeric(AGE.dcs),2), GoDARTS = round(as.numeric(AGE.godarts),2))
    out
  }
  
  countsN <- data.frame(Variable = c("N",
                                     "%Males"),
                        DCS = c(SEX$counts$dcs["Total",],
                                SEX$percentage$dcs["M",]),
                        GoDARTS = c(SEX$counts$godarts["Total",],
                                    SEX$percentage$godarts["M",])
  )
  
  
  clinVars <- rbind(
    .getVar("AGE","Age"),
    .getVar("VSTESTCD.BMI","BMI"),
    .getVar("LBTESTCD.HDL","HDL"),
    .getVar("LBTESTCD.HBA1C","HbA1c"),
    .getVar("LBTESTCD.CPEPTIDE","C-peptide"),
    .getVar("DIABETESDURATION","Diabetes duration") 
  )
  
  Gluc.dcs <- ds.table1D('ModelData$glucoselowering', datasources = opal["dcs"], type = "split")
  Gluc.godarts <- ds.table1D('ModelData$glucoselowering', datasources = opal["godarts"], type = "split")
  
  countsGluc <- data.frame(Variable = "Glucose lowering",
                           DCS = Gluc.dcs$percentage$dcs["1",],
                           GoDARTS = Gluc.godarts$percentage$godarts["1",])
  out <- rbind(countsN, clinVars, countsGluc)
  out
}