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

gatherResults(project=c("Pasin","GaussianPasin","Naveau","GaussianPasin1000"))
getResults(project=c("Pasin","GaussianPasin","Naveau","GaussianPasin1000"))

# graphsGenerate(project="Naveau",buildMethod = c("stepAIC","lassoBIC","lassoFDP10"),JPEG = T,PNG=T)
# graphsGenerate(project="Naveau",buildMethod = c("stepAIC","elasticnet5","lassoFDP10"),JPEG = T,PNG=T)
# graphsGenerate(project="Pasin",buildMethod = c("stepAIC","lassoBIC","lassoFDP10"),JPEG = T,PNG=T)
# graphsGenerate(project="Pasin",buildMethod = c("stepAIC","elasticnet5","lassoFDP10"),JPEG = T,PNG=T)
# graphsGenerate(project="GaussianPasin",buildMethod = c("stepAIC","lassoBIC","lassoFDP10"),JPEG = T,PNG=T)
# graphsGenerate(project="GaussianPasin",buildMethod = c("stepAIC","elasticnet5","lassoFDP10"),JPEG = T,PNG=T)
# graphsGenerate(project="GaussianPasin",buildMethod = paste0("elasticnet",1:9),JPEG = T,PNG=T)
# graphsGenerate(project="GaussianPasin",buildMethod = c("stepAIC","lassoSSrepFDP10","lassoFDP10"),JPEG = T,PNG=T)
# graphsGenerate(project="GaussianPasin1000",buildMethod = c("lassoFDP10"),JPEG = T,PNG=T)


pN0 = graphsGenerate(project=c("Naveau"),buildMethod = c("stepAIC","lassoFDP10"),JPEG = T,PNG=T)
pN = graphsGenerate(project=c("Naveau"),buildMethod = c("stepAIC","lassoFDP10","SAEMVS"),JPEG = T,PNG=T)
pP = graphsGenerate(project="Pasin",buildMethod = c("stepAIC","lassoFDP10"),JPEG = T,PNG=T)
pG = graphsGenerate(project="GaussianPasin1000",buildMethod = "lassoFDP10",JPEG = T,PNG=T)
pG2 = graphsGenerate(project="GaussianPasin",buildMethod = c("stepAIC","lassoFDP10"),JPEG = T,PNG=T)

pN_all = graphsGenerate(project="Naveau",buildMethod = c("stepAIC","lassoFDP5","lassoFDP10","lassoFDP20"),JPEG = T,PNG=T)
pP_all = graphsGenerate(project="Pasin",buildMethod = c("stepAIC","lassoFDP5","lassoFDP10","lassoFDP20"),JPEG = T,PNG=T)
pG_all = graphsGenerate(project="GaussianPasin",buildMethod = c("stepAIC","lassoFDP5","lassoFDP10","lassoFDP20"),JPEG = T,PNG=T)

## Generate Final graphs : 
PNG <- TRUE
JPEG <- TRUE

# 
ggarrange(pN$StatsComp+theme(axis.title.x = element_blank(),
                             plot.background = element_rect(linewidth=0.7,color="black")),
          pG$StatsComp+theme(axis.title.x = element_blank(),
                             plot.background = element_rect(linewidth=0.7,color="black")),
          pP$StatsComp+theme(axis.title.x = element_blank(),
                             plot.background = element_rect(linewidth=0.7,color="black")),
          nrow=1,common.legend=TRUE,legend="bottom",labels=c("A","B","C"),widths=c(1.2,0.65,1))

if(PNG){
  ggsave("outputs/figures/finalFigures/Figure2_colored.png",
         height = 2400, width =   6000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/Figure2_colored.jpeg",
         height = 2400, width =   6000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/Figure2_colored.eps",
       height=3,width=9,device=cairo_ps)

# 
ggarrange(pN$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          pG$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          pP$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          ncol=1,labels=c("A","B","C"),heights = c(1.3,0.6,1))

if(PNG){
  ggsave("outputs/figures/finalFigures/Figure3_colored.png",
         height = 5250, width =   5000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/Figure3_colored.jpeg",
         height = 5250, width =   5000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/Figure3_colored.eps",
       width = 8.33, height = 8.75,units="in",device=cairo_ps)


# Supplementary Figure - additional graphs  -------------------------------
pN_all = graphsGenerate(project="Naveau",buildMethod = c("stepAIC","lassoFDP5","lassoFDP10","lassoFDP20"),JPEG = T,PNG=T)
pP_all = graphsGenerate(project="Pasin",buildMethod = c("stepAIC","lassoFDP5","lassoFDP10","lassoFDP20"),JPEG = T,PNG=T)
pG_all = graphsGenerate(project="GaussianPasin",buildMethod = c("stepAIC","lassoFDP5","lassoFDP10","lassoFDP20"),JPEG = T,PNG=T)

ggarrange(pN_all$StatsComp+theme(axis.title.x = element_blank(),
                                 plot.background = element_rect(linewidth=0.7,color="black")),
          pG_all$StatsComp+theme(axis.title.x = element_blank(),
                                 plot.background = element_rect(linewidth=0.7,color="black")),
          pP_all$StatsComp+theme(axis.title.x = element_blank(),
                                 plot.background = element_rect(linewidth=0.7,color="black")),
          ncol=1,nrow=3,labels = c("A","B","C"),common.legend = TRUE,legend="bottom")

if(PNG){
  ggsave("outputs/figures/finalFigures/SuppFigure23.png",
         height = 6000, width =   3000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/SuppFigure23.jpeg",
         height = 6000, width =   3000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/SuppFigure23.eps",
       height=9,width=6.5,device=cairo_ps)

# 
ggarrange(pN0$LLComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          pG2$LLComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          pP$LLComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          nrow=1,labels=c("A","B","C"),common.legend = TRUE,legend="bottom")

if(PNG){
  ggsave("outputs/figures/finalFigures/SuppFigure24.png",
         height = 2400, width =   6000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/SuppFigure24.jpeg",
         height = 2400, width =   6000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/SuppFigure24.eps",
       height=3,width=9,device=cairo_ps)

ggarrange(
  ggarrange(pN$TimeComp$time+theme(axis.title.x = element_blank()),pN$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black")),
  ggarrange(pG2$TimeComp$time+theme(axis.title.x = element_blank()),pG2$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black")),
  ggarrange(pP$TimeComp$time+theme(axis.title.x = element_blank()),pP$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black")),
  labels=c("A","B","C"),nrow=1,label.y=0,vjust=-1)

if(PNG){
  ggsave("outputs/figures/finalFigures/SuppFigure25.png",
         height = 3000, width =   9000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/SuppFigure25.jpeg",
         height = 3000, width =   9000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/SuppFigure25.eps",
       height=4.5,width=14,device=cairo_ps)


# LassoBIC ----------------------------------------------------------------

pN = graphsGenerate(project=c("Naveau"),buildMethod = c("lassoFDP10","lassoBIC"),JPEG = T,PNG=T)
pP = graphsGenerate(project="Pasin",buildMethod = c("lassoFDP10","lassoBIC"),JPEG = T,PNG=T)
pG = graphsGenerate(project="GaussianPasin",buildMethod = c("lassoFDP10","lassoBIC"),JPEG = T,PNG=T)

# 
ggarrange(pN$StatsComp+theme(axis.title.x = element_blank(),
                             plot.background = element_rect(linewidth=0.7,color="black")),
          pG$StatsComp+theme(axis.title.x = element_blank(),
                             plot.background = element_rect(linewidth=0.7,color="black")),
          pP$StatsComp+theme(axis.title.x = element_blank(),
                             plot.background = element_rect(linewidth=0.7,color="black")),
          nrow=1,common.legend=TRUE,legend="bottom",labels=c("A","B","C"),widths=c(1,1,1))

if(PNG){
  ggsave("outputs/figures/finalFigures/SuppFigure1.png",
         height = 2400, width =   6000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/SuppFigure1.jpeg",
         height = 2400, width =   6000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/SuppFigure1.eps",
       height=3,width=9,device=cairo_ps)

# 
ggarrange(pN$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black"),strip.text = element_text(size = 8)),
          pG$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black"),strip.text = element_text(size = 8)),
          pP$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black"),strip.text = element_text(size = 8)),
          ncol=1,labels=c("A","B","C"),heights = c(1,1,1))

if(PNG){
  ggsave("outputs/figures/finalFigures/SuppFigure2.png",
         height = 5525, width =   5000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/SuppFigure2.jpeg",
         height = 5525, width =   5000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggarrange(pN$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black"),strip.text = element_text(size = 10)),
          pG$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black"),strip.text = element_text(size = 10)),
          pP$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black"),strip.text = element_text(size = 10)),
          ncol=1,labels=c("A","B","C"),heights = c(1,1,1))
ggsave("outputs/figures/finalFigures/SuppFigure2.eps",
       width = 8.33, height = 8.75,units="in",device=cairo_ps)

# 
ggarrange(
  ggarrange(pN$TimeComp$time+theme(axis.title.x = element_blank()),pN$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black")),
  ggarrange(pG$TimeComp$time+theme(axis.title.x = element_blank()),pG$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black")),
  ggarrange(pP$TimeComp$time+theme(axis.title.x = element_blank()),pP$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black")),
  labels=c("A","B","C"),nrow=1,label.y=0,vjust=-1)

if(PNG){
  ggsave("outputs/figures/finalFigures/SuppFigure3.png",
         height = 3000, width =   8000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/SuppFigure3.jpeg",
         height = 3000, width =   8000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/SuppFigure3.eps",
       height=3,width=12,device=cairo_ps)


# lassoSSrep --------------------------------------------------------------

pG = graphsGenerate(project="GaussianPasin",buildMethod = c("lassoFDP10","lassoSSrepFDP10"),JPEG = T,PNG=T)

if(PNG){
  ggsave(plot = pG$StatsComp,"outputs/figures/finalFigures/SuppFigure4.png",
         height = 2400, width =   2000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave(plot=pG$StatsComp,"outputs/figures/finalFigures/SuppFigure4.jpeg",
         height = 2400, width =   2000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave(plot=pG$StatsComp,"outputs/figures/finalFigures/SuppFigure4.eps",
       height=3,width=3,device=cairo_ps)

# 

if(PNG){
  ggsave(plot=pG$ParComp, "outputs/figures/finalFigures/SuppFigure5.png",
         height = 2000, width =   5000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave(plot=pG$ParComp,"outputs/figures/finalFigures/SuppFigure5.jpeg",
         height = 2000, width =   5000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave(plot=pG$ParComp,"outputs/figures/finalFigures/SuppFigure5.eps",
       width = 8.33, height = 3.5,units="in",device=cairo_ps)

# 
ggarrange(pG$TimeComp$time+theme(axis.title.x = element_blank()),pG$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black"))

if(PNG){
  ggsave("outputs/figures/finalFigures/SuppFigure6.png",
         height = 3000, width =   2500,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/SuppFigure6.jpeg",
         height = 3000, width =   2500,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/SuppFigure6.eps",
       height=4,width=5,device=cairo_ps)


# ElasticNet --------------------------------------------------------------

pN = graphsGenerate(project=c("Naveau"),buildMethod = c("lassoFDP10","elasticnet5"),JPEG = T,PNG=T)
pP = graphsGenerate(project="Pasin",buildMethod = c("lassoFDP10","elasticnet5"),JPEG = T,PNG=T)
pG = graphsGenerate(project="GaussianPasin",buildMethod = c("lassoFDP10","elasticnet5"),JPEG = T,PNG=T)

pG_all = graphsGenerate(project="GaussianPasin",buildMethod = paste0("elasticnet",1:9),JPEG = T,PNG=T)

# 
ggarrange(pN$StatsComp+theme(axis.title.x = element_blank(),
                             plot.background = element_rect(linewidth=0.7,color="black")),
          pG$StatsComp+theme(axis.title.x = element_blank(),
                             plot.background = element_rect(linewidth=0.7,color="black")),
          pP$StatsComp+theme(axis.title.x = element_blank(),
                             plot.background = element_rect(linewidth=0.7,color="black")),
          nrow=1,common.legend=TRUE,legend="bottom",labels=c("A","B","C"),widths=c(1,1,1))

if(PNG){
  ggsave("outputs/figures/finalFigures/SuppFigure7.png",
         height = 2400, width =   6000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/SuppFigure7.jpeg",
         height = 2400, width =   6000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/SuppFigure7.eps",
       height=3,width=9,device=cairo_ps)

# 
ggarrange(pN$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          pG$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          pP$ParComp+theme(plot.background = element_rect(linewidth=0.7,color="black")),
          ncol=1,labels=c("A","B","C"),heights = c(1,1,1))

if(PNG){
  ggsave("outputs/figures/finalFigures/SuppFigure8.png",
         height = 5525, width =   5000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/SuppFigure8.jpeg",
         height = 5525, width =   5000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/SuppFigure8.eps",
       width = 8.33, height = 8.75,units="in",device=cairo_ps)

# 
ggarrange(
  ggarrange(pN$TimeComp$time+theme(axis.title.x = element_blank()),pN$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black")),
  ggarrange(pG$TimeComp$time+theme(axis.title.x = element_blank()),pG$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black")),
  ggarrange(pP$TimeComp$time+theme(axis.title.x = element_blank()),pP$TimeComp$iter,ncol=1)+theme(plot.background = element_rect(linewidth=0.7,color="black")),
  labels=c("A","B","C"),nrow=1,label.y=0,vjust=-1)

if(PNG){
  ggsave("outputs/figures/finalFigures/SuppFigure9.png",
         height = 3000, width =   8000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave("outputs/figures/finalFigures/SuppFigure9.jpeg",
         height = 3000, width =   8000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave("outputs/figures/finalFigures/SuppFigure9.eps",
       height=3,width=12,device=cairo_ps)


if(PNG){
  ggsave(plot=pG_all$ParComp,"outputs/figures/finalFigures/SuppFigure10.png",
         height = 5525, width =   5000,dpi=600, units = "px", bg='transparent',device=grDevices::png)
}
if(JPEG){
  ggsave(plot=pG_all$ParComp,"outputs/figures/finalFigures/SuppFigure10.jpeg",
         height = 5525, width =   5000,dpi=600, units = "px",device=grDevices::jpeg)
}
ggsave(plot=pG_all$ParComp,"outputs/figures/finalFigures/SuppFigure10.eps",
       width = 8.33, height = 8.75,units="in",device=cairo_ps)
