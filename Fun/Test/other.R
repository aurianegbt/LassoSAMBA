library(caret)
library(dplyr)

design = read.csv("design_train.csv") 
genes = read.csv("expr_train.csv") 

genesKeep = genes %>% filter(sample_id %in% design$sample_id)

data = cbind.data.frame(MFC=design[,"MFC",drop=FALSE],design[,c("age","sex")],genesKeep[genesKeep$sample_id==design$sample_id,-1])[,c("MFC","age","sex","CYP4F3","EIF1AY","FGF3","GRK3","H1.3" ,"H3C1","HERC5","KDM5D","LRRC2","MATN4","MEP1A","MLXIP","PFKFB1","PRKY","PUDP","RUFY1","TNIP3","TP53BP2","UTY" ,"UVRAG","XPA","ZBTB14","ZFY","ZNF207","ZNF75D")]

data$sex <- as.numeric(data$sex)

knnmodel = knnreg(x=data[,-1],y=data[,1])

designTEST = read.csv("design_test.csv")
genesTEST = read.csv("expr_test.csv") 
genesKeepTEST = genesTEST %>% filter(sample_id %in% designTEST$sample_id)

dataTEST = cbind(designTEST[,c("age","sex")],genesKeepTEST[genesKeepTEST$sample_id==designTEST$sample_id,-1]) [,c("age","sex","CYP4F3","EIF1AY","FGF3","GRK3","H1.3" ,"H3C1","HERC5","KDM5D","LRRC2","MATN4","MEP1A","MLXIP","PFKFB1","PRKY","PUDP","RUFY1","TNIP3","TP53BP2","UTY" ,"UVRAG","XPA","ZBTB14","ZFY","ZNF207","ZNF75D")]

dataTEST$sex <- as.numeric(dataTEST$sex)

prediction = predict(knnmodel,dataTEST)

resultPREDICT = cbind.data.frame(participant_id=designTEST$participant_id,pred=as.numeric(prediction))

resultPREDICT <- resultPREDICT[!(duplicated(resultPREDICT$participant_id)),]
resultPREDICT[resultPREDICT$pred<0,] <- 0
write.csv(resultPREDICT,"awfulpredKNN.csv",row.names = FALSE)
