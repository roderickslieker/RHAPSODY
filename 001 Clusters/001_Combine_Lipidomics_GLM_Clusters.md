    library(meta)
    library(reshape2)
    library(ggplot2)
    library(tidyverse)
    source("./003.Utils/forestplot.R")

    load("./002.Clusters.Data/Lipids GLM_models across clusters_FDB_ANDIS_OneVsAll.RData")
    load("./002.Clusters.Data/Lipids GLM_models across clusters_FDB_DCS_OneVsAll.RData")
    load("./002.Clusters.Data/Lipids GLM_models across clusters_FDB_GoDARTS_OneVsAll.RData")

    updateData <- function(data,varNames){
      data <- data[match(varNames, data$var),]
      data$p.adj <- p.adjust(data$p.value,'fdr')
      return(data)
    }

    tmp <- c(as.character(ANDIS.lipidomics.results.SIDD$var),
             as.character(DCS.lipidomics.results.SIDD$var),
             as.character(GoDARTS.lipidomics.results.SIDD$var))
    lipids <- names(table(tmp)[table(tmp) == 3])

    lipids.classes.Model_1 <- lipids.classes[match(lipids,lipids.classes$FDBName),]

    #Update DCS
    results.SIDD.DCS  <- updateData(DCS.lipidomics.results.SIDD,lipids)
    results.SIRD.DCS  <- updateData(DCS.lipidomics.results.SIRD,lipids)
    results.MOD.DCS   <- updateData(DCS.lipidomics.results.MOD,lipids)
    results.MARD.DCS  <- updateData(DCS.lipidomics.results.MARD,lipids)
    results.MARDH.DCS <- updateData(DCS.lipidomics.results.MARDH,lipids)

    #Update GoDARTS
    results.SIDD.GoDARTS  <- updateData(GoDARTS.lipidomics.results.SIDD,lipids)
    results.SIRD.GoDARTS  <- updateData(GoDARTS.lipidomics.results.SIRD,lipids)
    results.MOD.GoDARTS   <- updateData(GoDARTS.lipidomics.results.MOD,lipids)
    results.MARD.GoDARTS  <- updateData(GoDARTS.lipidomics.results.MARD,lipids)
    results.MARDH.GoDARTS <- updateData(GoDARTS.lipidomics.results.MARDH,lipids)

    #Update Andis
    results.SIDD.ANDIS  <- updateData(ANDIS.lipidomics.results.SIDD,lipids)
    results.SIRD.ANDIS  <- updateData(ANDIS.lipidomics.results.SIRD,lipids)
    results.MOD.ANDIS   <- updateData(ANDIS.lipidomics.results.MOD,lipids)
    results.MARD.ANDIS  <- updateData(ANDIS.lipidomics.results.MARD,lipids)
    results.MARDH.ANDIS <- updateData(ANDIS.lipidomics.results.MARDH,lipids)

    head(results.SIDD.DCS)
    head(results.SIDD.GoDARTS)


    pd <- data.frame(Cluster = factor(rep(c("SIDD","SIRD","MOD","MD","MDH"), each=nrow(results.SIDD.DCS)), levels=c("SIDD","SIRD","MOD","MD","MDH")), 
                     DCS = c(results.SIDD.DCS$Estimate, results.SIRD.DCS$Estimate, results.MOD.DCS$Estimate, results.MARD.DCS$Estimate, results.MARDH.DCS$Estimate),
                     GoDARTS = c(results.SIDD.GoDARTS$Estimate, results.SIRD.GoDARTS$Estimate, results.MOD.GoDARTS$Estimate, results.MARD.GoDARTS$Estimate, results.MARDH.GoDARTS$Estimate),
                     ANDIS = c(results.SIDD.ANDIS$Estimate, results.SIRD.ANDIS$Estimate, results.MOD.ANDIS$Estimate, results.MARD.ANDIS$Estimate, results.MARDH.ANDIS$Estimate))
                     
    a1 <- rbind(
      ggplotGrob(
        ggplot(pd, aes(x=DCS, y=GoDARTS))+
        geom_point()+
        facet_grid(~Cluster)+
        geom_vline(xintercept = 0)+
        geom_hline(yintercept = 0)+
        geom_smooth(method=lm)
      ),ggplotGrob(
        ggplot(pd, aes(x=DCS, y=ANDIS))+
        geom_point()+
        facet_grid(~Cluster)+
        geom_vline(xintercept = 0)+
        geom_hline(yintercept = 0)+
        geom_smooth(method=lm)
      ),ggplotGrob(
        ggplot(pd, aes(x=GoDARTS, y=ANDIS))+
        geom_point()+
        facet_grid(~Cluster)+
        geom_vline(xintercept = 0)+
        geom_hline(yintercept = 0)+
        geom_smooth(method=lm)
    ), size="first")

    pdf("Comparison of effect size_lipidomics.pdf", width=10, height=10)
    plot(a1)
    dev.off()

## Combine

    meta.SIDD.lipidomics <- CombineGLM(results.SIDD.DCS,results.SIDD.GoDARTS,results.SIDD.ANDIS,ColVar = "var",studyLabels = c("dcs","godarts","andis"))
    meta.SIRD.lipidomics <- CombineGLM(results.SIRD.DCS,results.SIRD.GoDARTS,results.SIRD.ANDIS,ColVar = "var",studyLabels = c("dcs","godarts","andis"))
    meta.MOD.lipidomics <- CombineGLM(results.MOD.DCS,results.MOD.GoDARTS,results.MOD.ANDIS,ColVar = "var",studyLabels = c("dcs","godarts","andis"))
    meta.MARD.lipidomics <- CombineGLM(results.MARD.DCS,results.MARD.GoDARTS,results.MARD.ANDIS,ColVar = "var",studyLabels = c("dcs","godarts","andis"))
    meta.MARDH.lipidomics <- CombineGLM(results.MARDH.DCS,results.MARDH.GoDARTS,results.MARDH.ANDIS,ColVar = "var",studyLabels = c("dcs","godarts","andis"))

## SIDD

    knitr::kable(meta.SIDD.lipidomics[meta.SIDD.lipidomics$fdr.random < 0.05,])
    table(meta.SIDD.lipidomics$fdr.random <= 0.05)
    table(meta.SIDD.lipidomics[meta.SIDD.lipidomics$Effect.random < 0,]$fdr.random <= 0.05)

    getPlot(cohorts = c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),var = "SM.34.2.2", title = "SM.34.2.2")
        
    #print(getPlot(c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),"PC.O.16.1.0.18.1.0"))

    #knitr::include_graphics(path = 'Z:/Roderick/001 Projects/Rhapsody/001.1 Federated analyses/002.Clusters.Combine/SM_34_2_2_SIDD_forest.png')
    knitr::kable(meta.SIRD.lipidomics[meta.SIRD.lipidomics$fdr.random < 0.05,])
    get.plot(results.SIDD.DCS,results.SIDD.GoDARTS,"Estimates SIDD",c("DCS","GoDARTS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)
    get.plot(results.SIDD.DCS,results.SIDD.ANDIS,"Estimates SIDD",c("DCS","ANDIS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)
    get.plot(results.SIDD.GoDARTS,results.SIDD.ANDIS,"Estimates SIDD",c("GoDARTS","ANDIS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)

## SIRD

    table(meta.SIRD.lipidomics$fdr.random < 0.05)
    table(meta.SIRD.lipidomics[meta.SIRD.lipidomics$Effect.random < 0,]$fdr.random <= 0.05)

    getPlot(c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),"PC.O.16.1.0.18.1.0")
    getPlot(c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),"PC.O.16.1.0.18.2.0")
    getPlot(c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),"PC.O.16.1.0.18.1.0")

    getPlot(cohorts = c("DCS_Melt","GoDARTS_Melt","Andis_Melt"), breaks = seq(1,2,by=0.2), var = "PC.O.16.0.0.18.1.0", title = "PC O-16:0;0/18:1;0")

    getPlot(cohorts = c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),  var = "TAG.51.3.0", title = "TAG 51:3;0")

    getPlot(cohorts = c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),  var = "CE.16.1.0", title = "CE.16.1.0")





    get.plot(results.SIRD.DCS,resualts.SIRD.ANDIS,"Estimates SIRD",c("DCS","ANDIS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)
    get.plot(results.SIRD.GoDARTS,results.SIRD.ANDIS,"Estimates SIRD",c("GoDARTS","ANDIS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)

## MOD

    knitr::kable(meta.MOD.lipidomics[meta.MOD.lipidomics$fdr.random < 0.05,])
    table(meta.MOD.lipidomics$fdr.random < 0.05)
    table(meta.MOD.lipidomics[meta.MOD.lipidomics$Effect.random < 0,]$fdr.random <= 0.05)

    getPlot(cohorts = c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),var = "TAG.52.4.0", title = "TAG.52.4.0")
    getPlot(cohorts = c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),var = "TAG.51.3.0", title = "TAG.51.3.0")

    getPlot(c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),"TAG.52.4.0")
    get.plot(results.MOD.DCS,results.MOD.GoDARTS,"Estimates MOD",c("DCS","GoDARTS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)
    get.plot(results.MOD.DCS,results.MOD.ANDIS,"Estimates MOD",c("DCS","ANDIS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)
    get.plot(results.MOD.GoDARTS,results.MOD.ANDIS,"Estimates MOD",c("GoDARTS","ANDIS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)

## MARD

    knitr::kable(meta.MARD.lipidomics[meta.MARD.lipidomics$fdr.random < 0.05,])
    table(meta.MARD.lipidomics$fdr.random < 0.05)

    get.plot(results.MARD.DCS,results.MARD.GoDARTS,"Estimates MARD",c("DCS","GoDARTS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)
    get.plot(results.MARD.DCS,results.MARD.ANDIS,"Estimates MARD",c("DCS","ANDIS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)
    get.plot(results.MARD.GoDARTS,results.MARD.ANDIS,"Estimates MARD",c("GoDARTS","ANDIS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)

\#\#MARDH

    knitr::kable(meta.MARDH.lipidomics[meta.MARDH.lipidomics$fdr.random < 0.05,])
    meta.MARDH.lipidomics.sign <- meta.MARDH.lipidomics[meta.MARDH.lipidomics$fdr.random <= 0.05,]
    table(meta.MARDH.lipidomics.sign$fdr.random < 0.05)
    table(meta.MARDH.lipidomics.sign$Effect.random*-1 <= 0.05)

    getPlot(c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),"TAG.51.1.0")
    getPlot(c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),"PC.O.16.1.0.18.1.0")

    get.plot(results.MARDH.DCS,results.MARDH.GoDARTS,"Estimates MARDH",c("DCS","GoDARTS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)
    get.plot(results.MARDH.DCS,results.MARDH.ANDIS,"Estimates MARDH",c("DCS","ANDIS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)
    get.plot(results.MARDH.GoDARTS,results.MARDH.ANDIS,"Estimates MARDH",c("GoDARTS","ANDIS"),-1,1,-1,1,col=lipids.classes.Model_1$Class)

# Circos

    library(OmicCircos)
    lipids.classes <- read.table("../002.Data/SwissLipids.csv", sep=",", header=T)
    lipids.classes <- lipids.classes[lipids.classes$FDBName %in% lipids,]

    meta.SIDD.lipidomics.o <- meta.SIDD.lipidomics[match(lipids.classes$FDBName,meta.SIDD.lipidomics$var),]
    meta.SIRD.lipidomics.o <- meta.SIRD.lipidomics[match(lipids.classes$FDBName,meta.SIRD.lipidomics$var),]
    meta.MOD.lipidomics.o <- meta.MOD.lipidomics[match(lipids.classes$FDBName,meta.MOD.lipidomics$var),]
    meta.MARD.lipidomics.o <- meta.MARD.lipidomics[match(lipids.classes$FDBName,meta.MARD.lipidomics$var),]
    meta.MARDH.lipidomics.o <- meta.MARDH.lipidomics[match(lipids.classes$FDBName,meta.MARDH.lipidomics$var),]

    #New df
    pd <- data.frame(
                seg.name = lipids.classes$Class,
                seg.Start = seq(1,nrow(lipids.classes) * 3, by=3),
                seg.End = seq(3,nrow(lipids.classes) * 3, by=3),
                Name = lipids.classes$Class,
                Description = lipids.classes$Class)


    segv <- data.frame(
                seg.name = lipids.classes$Class,
                seg.po = seq(2,nrow(lipids.classes)*3, by=3),
                Gene = meta.SIDD.lipidomics.o$var,
                SIDD = meta.SIDD.lipidomics.o$Effect.random * -1,
                SIRD = meta.SIRD.lipidomics.o$Effect.random * -1,
                MOD = meta.MOD.lipidomics.o$Effect.random * -1,
                MARD = meta.MARD.lipidomics.o$Effect.random * -1,
                MARDH = meta.MARDH.lipidomics.o$Effect.random * -1
    )

    head(segv)



    source("./003.Utils/OmicCircos.R")

    options(stringsAsFactors = F)
    # load the OmicCircos-package
    db<-segAnglePo(pd,seg=unique(pd$seg.name));

    list.in <- table(pd$seg.name)
    list.in <- list.in[db[,1]]
    getDegree <- function(i, list.in){
      if(i==1 )
      {
       start <- 1 
      }else{
       start <- cumsum(list.in)[i-1]*2.011696+(i-1)*2
      }
      end <- cumsum(list.in)[i]*2.011696+(i-1)*2.011696
      data.frame(name = names(list.in)[i], i, start, end)
    }
    degrees <- do.call(rbind, lapply(1:9, getDegree, list.in))
    degrees <- degrees[match(db[,1],degrees$name),]
    db[,2] <- degrees$start + 270
    db[,3] <- degrees$end + 270

    pdf("./002.Clusters.Combine/Heatmap plot lipids.pdf", width=7, height=6.75)
    par(mar=c(2,2,2,2));
    plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="",main="");
    circos2(R=400,type="chr",cir=db,print.chr.lab=TRUE,W=4,scale=FALSE,col=viridis::viridis_pal()(9));
    circos2(R=250,cir=db,W=150,mapping=segv,col.v=4,type="heatmap2",lwd=1, B=F, col.bar=T, col.bar.po="bottomright", col=viridis::viridis_pal()(4))
    circos2(R=410,cir=db,W=20,mapping=segv,type="label",side="out",cex=0.4, col="black")
    dev.off()

# Enrichments

    lipidClass <- read.table("./002.Data/SwissLipids.csv", sep=",", header=T)
    lipidClass.subset <- lipidClass[lipidClass$FDBName %in% meta.SIDD.lipidomics$var,]



    pd <- data.frame(table(lipidClass.subset$Class))
    pd$Var1 <- as.character(pd$Var1) %>% factor(levels = pd$Var1[order(pd$Freq, decreasing = T)])

    pd$fraction = pd$Freq / sum(pd$Freq)
    pd = pd[order(pd$fraction,decreasing = T), ]
    pd$ymax = cumsum(pd$fraction)
    pd$ymin = c(0, head(pd$ymax, n=-1))

    pd$Var1 <- as.character(pd$Var1) %>% factor(levels = pd$Var1[order(pd$Freq, decreasing = T)])

    library(tidyverse)

    #Dataframw
    px <- ggplot(pd, aes(fill=Var1, ymax=ymax, ymin=ymin, xmax=4, xmin=3))+
      geom_rect() +
      coord_polar(theta="y") +
      xlim(c(0, 4)) +
      #theme_sniper() +
      theme(panel.grid=element_blank(),
            axis.title = element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            axis.line = element_blank()) +
      geom_text(aes(label=sprintf("%s\n(N=%s,%%=%s)",pd$Var1, pd$Freq, round(pd$fraction*100,2)), x=3.5,
                    y=pd$ymax-((pd$ymax-pd$ymin)/2)),
                col="black", check_overlap = FALSE)+
      scale_fill_manual(values = viridis::viridis_pal()(9))

    pdf("Figure Sx.pdf")
    px
    dev.off()

    enrichments <- rbind(getTable2("SIDD"),getTable2("SIRD"),getTable2("MOD"),getTable2("MARDH"))
    #Add labels
    enrichments$label <- paste0(round(enrichments$perc*100,1),"%")
    enrichments$label <- ifelse(enrichments$label == '0%', NA, enrichments$label)

    enrichments$ymin <- ifelse(enrichments$OR <= 1, enrichments$CI.lower, enrichments$OR)
    enrichments$ymax <- ifelse(enrichments$OR <= 1, enrichments$OR, enrichments$CI.upper)

    ymax <- max(enrichments$ymax[!is.infinite(enrichments$ymax)])
    enrichments2 <- data.frame(enrichments, colsplit(enrichments$name, "\\.", LETTERS[1:2]))
    enrichments2$B  <- factor(enrichments2$B, c("down","up"))
    enrichments2$A  <- factor(enrichments2$A, c("SIDD","SIRD","MOD","MARDH"))


    p1 <- ggplot(enrichments2, aes(x=class, y=OR, fill=class, label=label))+
      geom_col()+
      facet_grid(B~A)+
      scale_y_continuous(trans="log2", limits=c(0.01,ymax), breaks=c(0.01,0.25,1,100,1000))+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      geom_hline(yintercept = 1)+
      ylab("Odds ratio (log2)")+
      geom_text()+
      geom_errorbar(aes(ymin=ymin, ymax=ymax, col=class, width=0), alpha=0.5)+
      scale_fill_manual(values = viridis::viridis_pal()(9))+
      scale_colour_manual(values = viridis::viridis_pal()(9))

    pdf("./002.Clusters.Combine//Lipids_enrichment_classes_clusters.pdf", height=6, width=10)
    print(p1)
    dev.off()

    rio::export(enrichments2, file="./002.Clusters.Combine/Enrichment_Lipids.xlsx")


    p2 <- ggplot(enrichments2, aes(x=class, y=perc*100, fill=class, label=label))+
      geom_col()+
      facet_grid(B~A)+
      #scale_y_continuous(trans="log2", limits=c(0.01,ymax), breaks=c(0.01,0.25,1,100,1000,5000))+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      geom_hline(yintercept = 1)+
      ylab("Percentage in group")+
      geom_text(aes(y=(perc*100)+5))+
      #geom_errorbar(aes(ymin=ymin, ymax=ymax, col=class, width=0), alpha=0.5)+
      scale_fill_manual(values = viridis::viridis_pal()(9))+
      scale_colour_manual(values = viridis::viridis_pal()(9))


    pdf("./002.Clusters.Combine//Lipids_percentage_classes_clusters.pdf", height=6, width=10)
    print(p2)
    dev.off()

    enrichments

    rio::export(enrichments, file="./002.Clusters.Combine/Enrichments.xlsx")

    head(meta.MARDH.lipidomics)
    head(meta.SIRD.lipidomics, 20)
    getPlot(cohorts = c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),var = "TAG.48.0.0", title = "TAG.48.0.0")
    getPlot(cohorts = c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),var = "TAG.51.3.0", title = "TAG.51.3.0")

    getPlot(cohorts = c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),var = "PC.O.16.0.0.18.1.0", title = "PC.O.16.0.0.18.1.0", breaks=c(1.5, 2, 2.5,3,3.5))




    getPlot(cohorts = c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),var = "TAG.48.0.0", title = "TAG 48:0;0")
    getPlot(cohorts = c("DCS_Melt","GoDARTS_Melt","Andis_Melt"),var = "TAG.56.5.0", title = "TAG 56:5;0")
