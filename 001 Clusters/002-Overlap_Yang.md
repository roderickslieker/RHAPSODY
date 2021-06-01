## Packages

    library(magrittr)
    library(ggplot2)

    load("D:/001_Projects/033_RHAPSODY_Shiny/www/Cluster_data.RData")
    meta <- rio::import("D:/001_Projects/001_RHAPSODY/002.Data/MetaData.txt")
    corr <- rio::import("Table S2.xlsx")
    hr.cvd <- rio::import("Incident_CVD.xlsx")
    hr.cvd$Protein <- gsub("\r\n","", hr.cvd$Protein)

    pd <- protein.data[protein.data$FDR <= 0.05,]
    table(pd$Cluster)

    ## 
    ## SIDD SIRD  MOD   MD  MDH 
    ##    8  367  261   61  270

    clusters <- c("SIDD","SIRD","MOD","MD","MDH")


    corr.data <- lapply(clusters, function(cluster, pd, meta, corr){
      cat(cluster)
      pdx <- pd[pd$Cluster %in% cluster,]
      metax <- meta[match(pdx$Variable, meta$SomaId),]
      corrx <- corr[corr$Protein %in% metax$TargetFullName,]
      metax2 <- meta[match(corrx$Protein, meta$TargetFullName),]
      corrx$SomaId <- metax2$SomaId
      pdx2 <- pdx[match(corrx$SomaId, pdx$Variable),]
      corrx$Effect <- pdx2$Effect
      corrx$Direction <- ifelse(corrx$Effect <= 0,"Downregulated","Upregulated")
      
      corrx$Cluster <- cluster
      #corrx$Direction <- direction
      corrx
    }, pd=pd, meta=meta, corr=corr) %>% do.call(what = rbind)

    ## SIDDSIRDMODMDMDH

    colnames(corr.data)[3] <- "Correlation"
    corr.data$Cluster <-factor(corr.data$Cluster, levels=clusters)
    corr.data$Target <- meta[match(corr.data$SomaId, meta$SomaId),"Target"]
    corr.data$label <- ifelse(abs(corr.data$Correlation) >= 0.6, corr.data$Target,NA)

    ggplot(corr.data, aes(x=Cluster, y=Correlation, fill=Direction, label=label))+
      geom_boxplot(alpha=0.5, width=0.5)+
      ggbeeswarm::geom_quasirandom(alpha=0.5)+
      geom_hline(yintercept = 0)+
      scale_fill_manual(values = c("#009AC7","#900D09BB"))+
      ylab("Correlation with eGFR (Spearman rho)")+
      facet_grid(~Direction)+
      ggrepel::geom_text_repel()

    ## Warning: Removed 660 rows containing missing values (geom_text_repel).

![](Overlap_Yang_files/figure-markdown_strict/unnamed-chunk-4-1.png)

## Export

    pdf("Correlation data with eGFR_Yang paper.pdf", width=12, height=5)
    ggplot(corr.data, aes(x=Cluster, y=Correlation, fill=Direction, label=label))+
      geom_boxplot(alpha=0.5, width=0.5)+
      ggbeeswarm::geom_quasirandom(alpha=0.5)+
      geom_hline(yintercept = 0)+
      scale_fill_manual(values = c("#009AC7","#900D09BB"))+
      ylab("Correlation with eGFR (Spearman rho)")+
      facet_grid(~Direction)+
      ggrepel::geom_text_repel()

    ## Warning: Removed 660 rows containing missing values (geom_text_repel).

    dev.off()

    ## png 
    ##   2

# CVD

    head(hr.cvd)

    ##                                                             Protein UniProt ID
    ## 1                     15-hydroxyprostaglandindehydrogenase [NAD(+)]     P15428
    ## 2                  6-phosphogluconatedehydrogenase, decarboxylating     P52209
    ## 3 A disintegrin and metalloproteinase with thrombospondin motifs 13     Q76LX8
    ## 4                                              Adapter molecule crk     P46108
    ## 5                                                            Afamin     P43652
    ## 6                                   Alcohol dehydrogenase [NADP(+)]     P14550
    ##     HR  p-value Bonferroni corrected p- value HR\r\ncorrected for eGFR
    ## 1 0.72 1.15e-09                      1.21e-06                     0.79
    ## 2 0.80 3.44e-06                      3.63e-03                     0.92
    ## 3 0.77 1.22e-09                      1.29e-06                     0.87
    ## 4 1.37 2.20e-16                      2.20e-16                     1.17
    ## 5 0.80 1.89e-09                      1.99e-06                     0.87
    ## 6 1.24 3.45e-06                      3.63e-03                     1.13
    ##   p-value (HR corrected\r\nfor eGFR)
    ## 1                          6.430e-06
    ## 2                          6.900e-02
    ## 3                          5.530e-03
    ## 4                          7.316e-03
    ## 5                          1.620e-03
    ## 6                          1.280e-02
    ##   p-value\r\n(HR corrected for eGFR and Bonferroni corrected)
    ## 1                                                     0.00678
    ## 2                                                     1.00000
    ## 3                                                     1.00000
    ## 4                                                     1.00000
    ## 5                                                     1.00000
    ## 6                                                     1.00000

    cvd.data <- lapply(clusters, function(cluster, pd, meta, hr.cvd){
      cat(cluster)
      pdx <- pd[pd$Cluster %in% cluster,]
      metax <- meta[match(pdx$Variable, meta$SomaId),]
      hr.cvdx <- hr.cvd[hr.cvd$Protein %in% metax$TargetFullName,]
      metax2 <- meta[match(hr.cvdx$Protein, meta$TargetFullName),]
      hr.cvdx$SomaId <- metax2$SomaId
      pdx2 <- pdx[match(hr.cvdx$SomaId, pdx$Variable),]
      hr.cvdx$Effect <- pdx2$Effect
      hr.cvdx$Direction <- ifelse(hr.cvdx$Effect <= 0,"Downregulated","Upregulated")
      
      hr.cvdx$Cluster <- cluster
      #corrx$Direction <- direction
      hr.cvdx
    }, pd=pd, meta=meta, hr.cvd=hr.cvd) %>% do.call(what = rbind)

    ## SIDDSIRDMODMDMDH

    colnames(cvd.data) <- gsub("\r\n","",colnames(cvd.data))
    cvd.data$Cluster <-factor(cvd.data$Cluster, levels=clusters)
    cvd.data$Target <- meta[match(cvd.data$SomaId, meta$SomaId),"Target"]
    cvd.data$label <- ifelse(abs(log2(cvd.data$HR)) >= log2(1.55), cvd.data$Target,NA)

    px <- ggplot(cvd.data, aes(x=Cluster, y=HR, fill=Direction, label=label))+
      geom_boxplot(alpha=0.5, width=0.5)+
      ggbeeswarm::geom_quasirandom(alpha=0.5)+
      geom_hline(yintercept = 1)+
      scale_fill_manual(values = c("#009AC7","#900D09BB"))+
      ylab("Hazard ratio (log2)")+
      facet_grid(~Direction)+
      scale_y_continuous(trans="log2", breaks=c(0.5,.75,1,1.25,1.5))+
      ggrepel::geom_text_repel()

    plot(px)

    ## Warning: Removed 205 rows containing missing values (geom_text_repel).

![](Overlap_Yang_files/figure-markdown_strict/unnamed-chunk-7-1.png)

    pdf("CVD HR_Yang paper.pdf", width=12, height=5)
    px

    ## Warning: Removed 205 rows containing missing values (geom_text_repel).

    dev.off()

    ## png 
    ##   2

    rio::export(corr.data, file="eGFR_correlation_Yang.xlsx")
    rio::export(cvd.data, file="CVD_data_Yang.xlsx")
