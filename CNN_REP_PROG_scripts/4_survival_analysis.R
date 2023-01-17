
# A partir des features obtenu au CNN2, on recherche les features pronostique
# 13/10/2022

library(stringr)
library(survival)
library(ggplot2)
library(survminer)
library('gtsummary')
library('FactoMineR')
library(maxstat)
require(dplyr)
library(factoextra)
library(glmnet)

# mettre en forme les données cnn pour l'analyse de survie qui suit
########################################## CNN DENSENET PROG REP ###############################################
#densenet_v6_julieDatagen
#densenet_julieDatagen
features_path <- '/work/shared/ptbc/CNN_Pancreas_V2/results/CP/macenko_prog_rep/Feature_extration_prog_rep/without_color_aug/'
predictions_path <- '/work/shared/ptbc/CNN_Pancreas_V2/results/CP/macenko_prog_rep/Predictions_normalVStumor/' 
myFiles <- list.files(features_path)
length(myFiles)
myTable_mean <- patient <- c()
#test <- read.csv(paste0(features_path,myFiles[1]), header = F)


for (i in 1:length(myFiles)){
  print(i)
  nom_patient <- gsub("Features_","",myFiles[i]) #enleve Features du debut
  prediction <- read.csv(paste0(predictions_path,nom_patient),header = F) #myPrediction_
  prediction$pred_class <- ifelse(prediction$V2>0.5,"Normal","Tumor")
  nom_patient <- tools::file_path_sans_ext(nom_patient) #get namefile without extension (.csv)
  #nom_patient <- str_sub(nom_patient,1,12) # if cohor is TCGA
  patient <- c(patient, nom_patient) #concatene avec le vecteur patient
  myTable <- read.csv(paste0(features_path,myFiles[i]), header = F) #lire les features
  rownames(myTable) <- myTable$V1
  myTable <- myTable[prediction$V1[which(prediction$pred_class== "Tumor"&prediction$V3>0.95)],] #select les tuiles tumoral>95% de ma table des features
  myTable_mean <- rbind(myTable_mean,apply(myTable[2:1025],2,function(x){mean(x,na.rm = T)})) # faire un pooling des tuiles sur leurs moyennes respectives 
}




myTable_mean <- as.data.frame(myTable_mean)
myTable_mean$patient <- patient

#myTable_mean$patient<-str_sub(myTable_mean$patient,1,12)
write.csv(myTable_mean,"/work/shared/ptbc/CNN_Pancreas_V2/results/CP/macenko_prog_rep/data_related_survival/without_color_aug.csv", row.names = FALSE)


# On associe la clinique de bsc
#clinique <- read.csv("/work/shared/ptbc/CNN_Pancreas_V2/database/TCGA/1-s2.0-S0092867418302290-mmc1.csv")
#clinique <- clinique[,c("HES_x","DCD","OS","PROG","PFS")]
#colnames(clinique)[1] <- "patient"
#clinique$PFS <- as.numeric(as.vector(clinique$PFS))#/30.41
#clinique$OS <- as.numeric(as.vector(clinique$OS))#/30.41
#clinique$PROG <- as.numeric(as.vector(clinique$PROG))
#clinique$DCD <- as.numeric(as.vector(clinique$DCD))

# On associe la cliniquede tcga

  
clinique <- read.csv("/work/shared/ptbc/CNN_Pancreas_V2/database/canceropol/cp_clinical_data.csv")
clinique <- clinique[,c("patient","os","osev","dfs","dfsev")]
colnames(clinique)[1] <- "patient"
myTable_mean <- merge(clinique,myTable_mean)
write.csv(myTable_mean,"/work/shared/ptbc/CNN_Pancreas_V2/results/CP/macenko_prog_rep/data_related_survival/without_color_aug.csv")




clinique$PFI.time <- as.numeric(as.vector(clinique$PFI.time))/30.41
clinique$PFI <- as.numeric(as.vector(clinique$PFI))

#myTable_mean<-read.csv("/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/data_related_survival/TCGA_mean_macenko.csv")
myTable_mean <- merge(clinique,myTable_mean)
write.csv(myTable_mean,"/work/shared/ptbc/CNN_Pancreas_V2/results/TCGA/macenko_prog_rep/data_related_survival/mean_densenet_v6_2.csv", row.names = FALSE)


fit.null <- survfit(Surv(PFI.time, PFI) ~ 1, data =clinique[which(clinique$type =="PAAD"),])
ggsurvplot(fit.null, pval=FALSE,risk.table = TRUE, surv.median.line = 'h',xlab ='PFI _ Months', ggtheme = theme_grey())


######################################













