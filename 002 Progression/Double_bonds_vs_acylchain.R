setwd("D:/001_Projects/001_RHAPSODY/")

df <- rio::import("./Supplementary Data 1.xlsx",skip=3)

head(df)


df.tag <- df[grep("TAG", df$Lipid),]
head(df.tag)

library(magrittr) 

dfx <- data.frame(Lipid = df.tag$Lipid, 
                  HR = df.tag$HR...2,
                  Lower = df.tag$Lower...3,
                  Upper = df.tag$Upper...4)

dfx$Acyl <- (reshape2::colsplit(df.tag$Lipid, " ", LETTERS[1:2])[,2] %>% reshape2::colsplit(pattern = ";", names = c("Acyl", "rest")))[,1]

temp <- (reshape2::colsplit(df.tag$Lipid, " ", LETTERS[1:2])[,2] %>% reshape2::colsplit(pattern = ";", names = c("Acyl", "rest")))[,2] %>% reshape2::colsplit(pattern = ":", names = LETTERS[1:2])
dfx$Double <- temp$A
dfx$Pos <- temp$B


library(ggplot2)
library(patchwork)

px <- ggplot(dfx, aes(x=Acyl, y=HR))+
  geom_point()+
  geom_smooth(method = lm)+
  ylab("Hazard ratio")+
  xlab("Acyl chain length")+
  ggtitle("Acyl chain length")+
  scale_y_continuous(trans="log10", breaks=seq(1,1.2, by=0.05))+
ggplot(dfx, aes(x=Double, y=HR))+
  geom_point()+
  geom_smooth(method = lm)+
  ylab("Hazard ratio")+
  xlab("Number of double bonds")+
  ggtitle("Number of double bonds")+
  scale_y_continuous(trans="log2", breaks=seq(1,1.2, by=0.05))+
  ggplot(dfx, aes(x=Acyl, y=Double, col=HR, group=Lipid))+
  geom_point(size=2)+
  geom_smooth(method = lm)+
  ylab("Number of double bonds")+
  xlab("Acyl chain length")+
  ggtitle("Acyl chain length vs # double bonds")+
  scale_colour_gradientn(colours = MetBrewer::met.brewer("Cross", 10))


pdf("Double bonds vs acyl chains.pdf", width=12, height=4)
print(px)
dev.off()