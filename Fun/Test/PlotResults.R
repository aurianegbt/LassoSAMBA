Res=Reduce(rbind,lapply(1:50,function(i){get(load(paste0("Fun/Test/ResultsNew/Res_",i,".RData")))}))

library(ggplot2)
library(ggpubr)
library(dplyr)
source("~/Travail/00_Theme.R")


## Comparaison avec et sans blanchiment quand aucune corrélation 

ggplot(Res[Res$corrInY==1,],aes(x=Method,y=FDR,color=as.factor(Wh)))+
  geom_boxplot() +  
  ggtitle("FDR comparison between Lasso and Elastic Net selection,\nwhen no correlation is considered between the parameters.", subtitle="Among 50 simulated model.")+
  labs(color = "Whitening step")

ggplot(Res[Res$corrInY==2,],aes(x=Method,y=FDR,color=as.factor(Wh)))+
  geom_boxplot()+ 
  ggtitle("FDR comparison between Lasso and Elastic Net selection,\nwhen correlation is considered between the parameters.", subtitle="Among 50 simulated model.")+
  labs(color = "Whitening step")


ggplot(Res[Res$corrInY==1,],aes(x=Method,y=FNR,color=as.factor(Wh)))+
  geom_boxplot() +  
  ggtitle("FNR comparison between Lasso and Elastic Net selection,\nwhen no correlation is considered between the parameters.", subtitle="Among 50 simulated model.")+
  labs(color = "Whitening step")



ggplot(Res[Res$corrInY==2,],aes(x=Method,y=FNR,color=as.factor(Wh)))+
  geom_boxplot()+ 
  ggtitle("FNR comparison between Lasso and Elastic Net selection,\nwhen correlation is considered between the parameters.", subtitle="Among 50 simulated model.")+
  labs(color = "Whitening step")

## Dans le cas où on blanchi, résultat des méthodes peu importe le modèle entran t: 
methKeep = c("EN","Lasso","Lasso2","Lasso2Crit","LassoCrit")
ggplot(Res[Res$Method %in% methKeep,],aes(x=Method,y=FNR,color=as.factor(corrInX),fill=Wh))+
  geom_violin(lwd=1,alpha=0.2) + 
  ggtitle("FDR comparison between Lasso and Elastic Net selection.", subtitle="Among 50 simulated model.")+
  scale_color_discrete(name = "Correlated covariates", labels = c("None", "Moderated","High"))

ggplot(Res[Res$Method %in% methKeep,],aes(x=Method,y=FNR,color=as.factor(corrInX),fill=Wh))+
  geom_violin(lwd=1,alpha=0.2) + 
  ggtitle("FNR comparison between Lasso and Elastic Net selection.", subtitle="Among 50 simulated model.")+
  scale_color_discrete(name = "Correlated covariates", labels = c("None", "Moderated","High"))

## Cas réaliste : high correlation for X and moderated for Y 
ggplot(Res[Res$corrInX==3 & Res$corrInY==2 & Res$Method %in% methKeep,],aes(x=Method,y=FDR,color=Wh))+
  geom_violin() + 
  ggtitle("FDR comparison between Lasso and Elastic Net selection\nwhen high correlated covariates and correlated parameters .", subtitle="Among 50 simulated model.")

ggplot(Res[Res$corrInX==3 & Res$corrInY==2 & Res$Method %in% methKeep,],aes(x=Method,y=FNR,color=Wh))+
  geom_violin() + 
  ggtitle("FDR comparison between Lasso and Elastic Net selection\nwhen high correlated covariates and correlated parameters .", subtitle="Among 50 simulated model.")

