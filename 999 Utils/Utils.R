get.plot <- function(dataDCS, datagodarts, title, cohortLabels, ymin, ymax, xmin, xmax)
{
  pd <- data.frame(Beta.DCS = dataDCS$Estimate, Beta.godarts = datagodarts$Estimate)
  pd$Sign <- ifelse(dataDCS$p_adjusted <= 0.05 & datagodarts$p_adjusted <= 0.05, "Sign", "NonSign")
  
  pd$label <- ifelse(pd$Sign == "Sign", NA, NA)
  
  minx <- round(min(c(pd$Beta.DCS, pd$Beta.godarts)), 1)
  maxx <- round(max(c(pd$Beta.DCS, pd$Beta.godarts)), 1)
  
  P <- ggplot(pd, aes(x=Beta.DCS, y=Beta.godarts, label=label))+
    geom_point()+
    geom_text()+
    #scale_colour_manual(values = viridis::viridis_pal()(length(unique(cols))))+
    xlab(paste("Estimate",cohortLabels[1]))+
    ylab(paste("Estimate",cohortLabels[2]))+
    ggtitle(title)+
    geom_smooth(method=lm, col='black')+
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0)+
    xlim(xmin, xmax)+ylim(ymin, ymax)

  return(P)
}


CombineGLM <- function(Data01, Data02, Data03, ColVar, studyLabels){
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
    var <- variables[i]
    tmp.merge <- rbind(Data01[Data01[,ColVar] == var,],
                       Data02[Data02[,ColVar] == var,],
                       Data03[Data03[,ColVar] == var,])
    

    tmp.combined <- metagen(TE = tmp.merge$Estimate, seTE = tmp.merge$Std..Error, studlab = studyLabels, comb.fixed = TRUE, comb.random = TRUE)
    cis <- ci(TE = tmp.combined$TE.random, seTE = tmp.combined$seTE.random, level=0.95, df=NULL, null.effect = 0)
    summ <- data.frame(var,
                      Effect.random = tmp.combined$TE.random,
                      lower = cis$lower,
                      upper = cis$upper,
                      Zval.random = tmp.combined$zval.random,
                      Pval.random  = tmp.combined$pval.random,
                      I2 = summary(tmp.combined)$I2$TE,
                      Het = pchisq(tmp.combined$Q,1, lower.tail=F))
    
    combined <- rbind(combined,summ)
  }
  combined <- combined[order(combined$Pval.random, decreasing=F),]
  combined$fdr.random <- p.adjust(combined$Pval.random,method="fdr")
  return(combined)
}


updateData <- function(data,varNames){
  data <- data[match(varNames, data$var),]
  data$p_adjusted <- p.adjust(data$p.value,'fdr')
  return(data)
}
updateData2 <- function(data,varNames){
  data <- data[match(varNames, data$var),]
  #data$p_adjusted <- p_adjustedust(data$p.value,'fdr')
  return(data)
}


cols <- c("#132B41","#009AC6","#E1B70A","#B81E3D","#878786")

getPlot <- function(cohorts, var, title=NULL, breaks=NULL){
  dataframes <- lapply(cohorts,function(co){get(co)[get(co)$lip == var,]})
  pd <- do.call(rbind,dataframes)
  pdx <- pd[pd$variable %in% c("25%","75%","50%"),]
  pdx.t <- pdx[pdx$variable %in% "25%",]
  pdx.t <- pdx.t[,-(c(ncol(pdx.t)-1):ncol(pdx.t))]
  pdx.t$v25 <- pdx[pdx$variable %in% "25%","value"]
  pdx.t$v50 <- pdx[pdx$variable %in% "50%","value"]
  pdx.t$v75 <- pdx[pdx$variable %in% "75%","value"]
  
  conv <- data.frame(id = 2:6,
                     name = c("1/SIDD","2/SIRD","3/MOD","4/MD","5/MDH"))
  pdx.t$Cluster <- conv[match(pdx.t$cluster, conv$id),"name"]
  
  pdx.t$cohort2 <- toupper(substr(pdx.t$cohort, 1, 1))
  pdx.t$cohort <- factor(pdx.t$cohort, levels=c("dcs","godarts","andis"))
  p1 <- ggplot(pdx.t, aes(x=factor(Cluster), y=v50))+
    geom_point(size=2, position=position_dodge2(width=  0.4), aes(col=Cluster, group=cohort))+
    geom_errorbar(aes(x=Cluster, ymin=v25, ymax=v75, col=Cluster, group=cohort), width=0, position=position_dodge(width = 0.4))+
    scale_y_continuous(trans="log10")+#, limits = c(min(pdx.t$v25)/1.5,max(pdx.t$v75)*1.5))+
    scale_colour_manual(values = cols)+
    #geom_text(aes(x=Cluster, y=v75*1.04, label=cohort2), position=position_dodge2(width=0.4))+
    ylab("Levels")+
    xlab("Cluster")+
    facet_wrap(~cohort, scale="free_y")+
    theme(legend.position = "none")
  
  if(!is.null(title))
  {
    p1 <- p1+ggtitle(title)
  }
  if(!is.null(breaks))
  {
    p1 <- p1+scale_y_continuous(trans="log10",breaks=breaks)
  }
  
  pdf(sprintf("../003.Figures.combined/Plot_of_%s.pdf", var), width=8, height=3)
  print(p1)
  dev.off()
  return(p1)
}


CombineGLMTwo <- function(Data01, Data02, ColVar, studyLabels){
  #Check if dataframes have equal number of rows
  if(!nrow(Data01) == nrow(Data02)){
    stop("Dataframes are not of equal length")
  }
  #Check if variables are equal and in same order
  #if(!sum(as.character(Data01[,ColVar]) == as.character(Data02[,ColVar]))==nrow(Data01)){
  #  stop("Variables of dataframes are not in same order")
  #}
  combined <- data.frame()
  variables <-as.character(intersect(na.omit(Data01[,ColVar]), na.omit(Data02[,ColVar])))
  
  for(i in seq_along(variables)){
    #cat(i)
    var <- variables[i]
    tmp.merge <- rbind(Data01[which(Data01[,ColVar] == var),],
                       Data02[which(Data02[,ColVar] == var),])
    
    
    tmp.combined <- metagen(TE = tmp.merge$Estimate, seTE = tmp.merge$Std..Error, studlab = studyLabels, comb.fixed = TRUE, comb.random = TRUE)
    cis <- ci(TE = tmp.combined$TE.random, seTE = tmp.combined$seTE.random, level=0.95, df=NULL, null.effect = 0)
    summ <- data.frame(var,
                       Effect.random = tmp.combined$TE.random,
                       lower = cis$lower,
                       upper = cis$upper,
                       Zval.random = tmp.combined$zval.random,
                       Pval.random  = tmp.combined$pval.random,
                       I2 = summary(tmp.combined)$I2$TE,
                       Het = pchisq(tmp.combined$Q,1, lower.tail=F))
    
    combined <- rbind(combined,summ)
  }
  combined$fdr.random <- p.adjust(combined$Pval.random,method="fdr")
  combined <- combined[order(combined$Pval.random, decreasing=F),]
  return(combined)
}