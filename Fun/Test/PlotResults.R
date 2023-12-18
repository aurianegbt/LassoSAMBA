Res=Reduce(rbind,lapply(1:50,function(i){get(load(paste0("Fun/Test/Results/Res_",i,".RData")))}))

library(ggplot2)
library(ggpubr)
library(dplyr)
source("~/Travail/00_Theme.R")

# Comparaison paramètre 

ggplot(Res,aes(x=Method,y=FDR,color=as.factor(corrInY)))+
  scale_color_manual(name="Parameters correlation",labels=c("False","True"),values=cbPalette[2:3])+
  geom_boxplot() +  
  ggtitle("FDR comparison between Lasso and Elastic Net selection", subtitle="Among 50 simulated model.")

ggplot(Res,aes(x=Method,y=FNR,color=as.factor(corrInY)))+
  scale_color_manual(name="Parameters correlation",labels=c("False","True"),values=cbPalette[2:3])+  
  geom_boxplot() +  
  ggtitle("FNR comparison between Lasso and Elastic Net selection", subtitle="Among 50 simulated model.")

# COmparaison covariables 
ggplot(Res,aes(x=Method,y=FDR,color=as.factor(corrInX)))+
  scale_color_manual(name="Covariates correlation",labels=c("None","Low","High"),values=cbPalette[c(4,3,2)])+
  geom_boxplot() +  
  ggtitle("FDR comparison between Lasso and Elastic Net selection", subtitle="Among 50 simulated model.")

ggplot(Res,aes(x=Method,y=FNR,color=as.factor(corrInX)))+
  scale_color_manual(name="Covariates correlation",labels=c("None","Low","High"),values=cbPalette[c(4,3,2)])+
  geom_boxplot() +  
  ggtitle("FNR comparison between Lasso and Elastic Net selection,\nwhen no correlation is considered between the parameters.", subtitle="Among 50 simulated model.")

# Both
ggplot(Res,aes(x=Method,y=FDR,color=as.factor(corrInX),fill=as.factor(corrInY)))+
  scale_color_manual(name="Covariates correlation",labels=c("None","Low","High"),values=cbPalette[c(4,3,2)])+
  scale_fill_manual(name="Parameters correlation",labels=c("False","True"),values=cbPalette[2:3])+  
  geom_boxplot(alpha=0.4,lwd=1) +  
  ggtitle("FDR comparison between Lasso and Elastic Net selection", subtitle="Among 50 simulated model.")

ggplot(Res,aes(x=Method,y=FNR,color=as.factor(corrInX),fill=as.factor(corrInY)))+
  scale_color_manual(name="Covariates correlation",labels=c("None","Low","High"),values=cbPalette[c(4,3,2)])+
  scale_fill_manual(name="Parameters correlation",labels=c("False","True"),values=cbPalette[2:3])+  
  geom_boxplot(alpha=0.4,lwd=1) +  
  ggtitle("FNR comparison between Lasso and Elastic Net selection", subtitle="Among 50 simulated model.")


## Cas réaliste : high correlation for X and moderated for Y 
ggplot(Res[Res$corrInX==3 & Res$corrInY==2,],aes(x=Method,y=FDR))+
  geom_violin() + 
  ggtitle("FDR comparison between Lasso and Elastic Net selection\nwhen high correlated covariates and correlated parameters .", subtitle="Among 50 simulated model.")

ggplot(Res[Res$corrInX==3 & Res$corrInY==2,],aes(x=Method,y=FNR))+
  geom_violin() + 
  ggtitle("FDR comparison between Lasso and Elastic Net selection\nwhen high correlated covariates and correlated parameters .", subtitle="Among 50 simulated model.")

