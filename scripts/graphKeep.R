library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggpattern)
library(webshot2)
library(ggh4x)
library(flextable)
library(grid)
library(gtable)
library(gridExtra)
library(scales)
sapply(list.files("scripts/resultsFun/graphsFun",full.names = T),FUN=function(d){source(d,echo=F)})
source("scripts/resultsFun/gatherResults.R")
source("scripts/resultsFun/getResults.R")
source("scripts/resultsFun/graphs.R")

gatherResults(project=c("Pasin","GaussianPasin","Naveau"))
getResults(project=c("Pasin","GaussianPasin","Naveau"))

pN0 = graphsGenerate(project="Naveau",buildMethod = c("stepAIC","lassoFDP10"),JPEG = T,PNG=T)
pP = graphsGenerate(project="Pasin",buildMethod = c("stepAIC","lassoFDP10"),JPEG = T,PNG=T)
pG = graphsGenerate(project="GaussianPasin",buildMethod = c("stepAIC","lassoFDP10"),JPEG = T,PNG=T)

pN = graphsGenerate(project=c("Naveau"),buildMethod = c("stepAIC","lassoFDP10","SAEMVS"),JPEG = T,PNG=T)

pN_all = graphsGenerate(project="Naveau",buildMethod = c("stepAIC","lassoFDP5","lassoFDP10","lassoFDP20"),JPEG = T,PNG=T)
pP_all = graphsGenerate(project="Pasin",buildMethod = c("stepAIC","lassoFDP5","lassoFDP10","lassoFDP20"),JPEG = T,PNG=T)
pG_all = graphsGenerate(project="GaussianPasin",buildMethod = c("stepAIC","lassoFDP5","lassoFDP10","lassoFDP20"),JPEG = T,PNG=T)

## Generate Final graphs : 
PNG <- JPEG <- TRUE

# 
ggarrange(pN$StatsComp+theme(axis.title.x = element_blank(),
                             plot.background = element_rect(linewidth=0.7,color="black")),
          pG$StatsComp+theme(axis.title.x = element_blank(),
                             plot.background = element_rect(linewidth=0.7,color="black")),
          pP$StatsComp+theme(axis.title.x = element_blank(),
                             plot.background = element_rect(linewidth=0.7,color="black")),
          nrow=1,common.legend=TRUE,legend="bottom",labels=c("A","B","C"),widths=c(1.3,1,1))

if(PNG){
  ggsave("outputs/figures/finalFigures/ComparisonStats.png",
         height = 2400, width =   6000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/ComparisonStats.jpeg",
         height = 2400, width =   6000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/ComparisonStats.eps",
       height=3,width=9,device=cairo_ps)

# 
ggarrange(pN$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          pG$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          pP$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          ncol=1,labels=c("A","B","C"),heights = c(1.4,1,1))

if(PNG){
  ggsave("outputs/figures/finalFigures/NumberSelectionParameter.png",
         height = 6000, width =   3500,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/NumberSelectionParameter.jpeg",
         height = 6000, width =   3500,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/NumberSelectionParameter.eps",
       height=11,width=6,device=cairo_ps)

# 
ggarrange(pN_all$StatsComp+theme(axis.title.x = element_blank(),
                                 plot.background = element_rect(linewidth=0.7,color="black")),
          pG_all$StatsComp+theme(axis.title.x = element_blank(),
                                 plot.background = element_rect(linewidth=0.7,color="black")),
          pP_all$StatsComp+theme(axis.title.x = element_blank(),
                                 plot.background = element_rect(linewidth=0.7,color="black")),
          ncol=1,labels = c("A","B","C"),common.legend = TRUE,legend="bottom")

if(PNG){
  ggsave("outputs/figures/finalFigures/ComparisonStats_all.png",
         height = 6000, width =   3000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/ComparisonStats_all.jpeg",
         height = 6000, width =   3000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/ComparisonStats_all.eps",
       height=9,width=4.5,device=cairo_ps)

# 
ggarrange(pN0$LLComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          pG$LLComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          pP$LLComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          nrow=1,labels=c("A","B","C"),common.legend = TRUE,legend="bottom")

if(PNG){
  ggsave("outputs/figures/finalFigures/LL.png",
         height = 2400, width =   6000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/LL.jpeg",
         height = 2400, width =   6000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/LL.eps",
       height=3,width=9,device=cairo_ps)

ggarrange(
  ggarrange(pN0$TimeComp$time+theme(axis.title.x = element_blank()),pN0$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black")),
  ggarrange(pG$TimeComp$time+theme(axis.title.x = element_blank()),pG$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black")),
  ggarrange(pP$TimeComp$time+theme(axis.title.x = element_blank()),pP$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black")),
  labels=c("A","B","C"),nrow=1,label.y=0,vjust=-1)

if(PNG){
  ggsave("outputs/figures/finalFigures/Time.png",
         height = 3000, width =   8000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/Time.jpeg",
         height = 3000, width =   8000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/Time.eps",
       height=4.5,width=14,device=cairo_ps)
