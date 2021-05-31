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

    pdf("Comparison of effect size_metabolomics.pdf", width=10, height=10)
    plot(a1)
    dev.off()

## SIDD

    getForestPlotGLM(Data01 = results.SIDD.DCS,
                     Data02 = results.SIDD.GoDARTS,
                     Data03 = results.SIDD.ANDIS,
                     Meta = meta.SIDD.metabolomics,
                     variable = "Tyr",
                     omic = "metabolomics")
    getForestPlotGLM(Data01 = results.SIDD.DCS,
                     Data02 = results.SIDD.GoDARTS,
                     Data03 = results.SIDD.ANDIS,
                     Meta = meta.SIDD.metabolomics,
                     variable = "SDMA.ADMA",
                     omic = "metabolomics")

    getForestPlotGLM(Data01 = results.SIDD.DCS,
                     Data02 = results.SIDD.GoDARTS,
                     Data03 = results.SIDD.ANDIS,
                     Meta = meta.SIDD.metabolomics,
                     variable = "SDMA.ADMA",
                     omic = "metabolomics")



    pdf("../003.Figures.combined/Forest plot Tyr.pdf", width=10, height=3)
    getForestPlotGLM(Data01 = results.SIDD.DCS,
                     Data02 = results.SIDD.GoDARTS,
                     Data03 = results.SIDD.ANDIS,
                     Meta = meta.SIDD.metabolomics,
                     variable = "Tyr",
                     omic = "metabolomics")
    dev.off()

## SIRD

    getForestPlotGLM(Data01 = results.SIRD.DCS,
                     Data02 = results.SIRD.GoDARTS,
                     Data03 = results.SIRD.ANDIS,
                     Meta = meta.SIRD.metabolomics,
                     variable = "Tyr",
                     omic = "metabolomics")
    getForestPlotGLM(Data01 = results.SIRD.DCS,
                     Data02 = results.SIRD.GoDARTS,
                     Data03 = results.SIRD.ANDIS,
                     Meta = meta.SIRD.metabolomics,
                     variable = "Leu",
                     omic = "metabolomics")

    getForestPlotGLM(Data01 = results.SIRD.DCS,
                     Data02 = results.SIRD.GoDARTS,
                     Data03 = results.SIRD.ANDIS,
                     Meta = meta.SIRD.metabolomics,
                     variable = "Ile",
                     omic = "metabolomics")

    getForestPlotGLM(Data01 = results.SIRD.DCS,
                     Data02 = results.SIRD.GoDARTS,
                     Data03 = results.SIRD.ANDIS,
                     Meta = meta.SIRD.metabolomics,
                     variable = "AADA",
                     omic = "metabolomics")

## MOD

    getForestPlotGLM(Data01 = results.MOD.DCS,
                     Data02 = results.MOD.GoDARTS,
                     Data03 = results.MOD.ANDIS,
                     Meta = meta.MOD.metabolomics,
                     variable = "Ile",
                     omic = "metabolomics")

## MARD

## MARDH

    getForestPlotGLM(Data01 = results.MARDH.DCS,
                     Data02 = results.MARDH.GoDARTS,
                     Data03 = results.MARDH.ANDIS,
                     Meta = meta.MARDH.metabolomics,
                     variable = "Ile",
                     omic = "metabolomics")

    getForestPlotGLM(Data01 = results.MARDH.DCS,
                     Data02 = results.MARDH.GoDARTS,
                     Data03 = results.MARDH.ANDIS,
                     Meta = meta.MARDH.metabolomics,
                     variable = "Taurine",
                     omic = "metabolomics")

    save(meta.SIDD.metabolomics,file="Metabolomics_Clusters_Meta_SIDD.RData")
    save(meta.SIRD.metabolomics,file="Metabolomics_Clusters_Meta_SIRD.RData")
    save(meta.MOD.metabolomics,file="Metabolomics_Clusters_Meta_MOD.RData")
    save(meta.MARD.metabolomics,file="Metabolomics_Clusters_Meta_MARD.RData")
    save(meta.MARDH.metabolomics,file="Metabolomics_Clusters_Meta_MARDH.RData")


    meta.SIDD.metabolomics$group <- "SIDD"
    meta.SIRD.metabolomics$group <- "SIRD"
    meta.MOD.metabolomics$group <- "MOD"
    meta.MARD.metabolomics$group <- "MD"
    meta.MARDH.metabolomics$group <- "MDH"


    out <- rbind(
      meta.SIDD.metabolomics[meta.SIDD.metabolomics$fdr.random <= 0.05,],
      meta.SIRD.metabolomics[meta.SIRD.metabolomics$fdr.random <= 0.05,],
      meta.MOD.metabolomics[meta.MOD.metabolomics$fdr.random <= 0.05,],
      meta.MARD.metabolomics[meta.MARD.metabolomics$fdr.random <= 0.05,],
      meta.MARDH.metabolomics[meta.MARDH.metabolomics$fdr.random <= 0.05,]
      )

    rio::export(out, file="Metabolomics clusters out.xlsx")

## Get heatmap of meta-estimates.

    vars <- metabolites

     #Make sign matrix
    sign.matrix <- matrix(rep(0), ncol=5, nrow=length(vars))
    colnames(sign.matrix) <- c("SIDD","SIRD","MOD","MARD","MARDH")
    rownames(sign.matrix) <- vars


    sign.matrix[match(meta.SIDD.metabolomics$var, rownames(sign.matrix)), 1] <- 1
    sign.matrix[match(meta.SIRD.metabolomics$var, rownames(sign.matrix)), 2] <- 1
    sign.matrix[match(meta.MOD.metabolomics$var, rownames(sign.matrix)), 3] <- 1
    sign.matrix[match(meta.MARD.metabolomics$var, rownames(sign.matrix)), 4] <- 1
    sign.matrix[match(meta.MARDH.metabolomics$var, rownames(sign.matrix)), 5] <- 1
    estimates <- data.frame(
              SIDD = meta.SIDD.metabolomics[match(vars, meta.SIDD.metabolomics$var),"Effect.random"],
              SIRD = meta.SIRD.metabolomics[match(vars, meta.SIRD.metabolomics$var),"Effect.random"],
              MOD = meta.MOD.metabolomics[match(vars, meta.MOD.metabolomics$var),"Effect.random"],
              MARD = meta.MARD.metabolomics[match(vars, meta.MARD.metabolomics$var),"Effect.random"],
              MARDH = meta.MARDH.metabolomics[match(vars, meta.MARDH.metabolomics$var),"Effect.random"])


    rownames(estimates) <- vars

    # Estimates are inverted 
    estimates <- estimates*-1

    # Functions for correlation
    hclust2 <- function(x, method="average", ...)
        hclust(x, method=method, ...)
    dist2 <- function(x, ...)
        as.dist(1-cor(t(x)))
      


    aheatmap(data.matrix(estimates), hclustfun = hclust2, distfun =dist2, Colv=NA,  
             annRow = sign.matrix[,-4], 
             annColors = c("#4DB3E6","#00B399","#E69900","#8B1A4F"), 
             color = colorRampPalette(colors =  c("#030189","white","#FD0000"))(50),
             breaks = seq(-1,1,length.out = 51))
