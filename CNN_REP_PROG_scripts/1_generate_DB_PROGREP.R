# ---------------------------------------------
# Programme: 
# Auteur EB
# Date: FEVRIER 2022
# Generation des bases de donnees train, validation et test
# pour le modele PROG et REP sur le besoncon
# Les tuiles selectionnees sont tumorales et avec une proba >0.95
# Si un patient a plusieurs lames une seule est choisie au hasard
# et sur cette lame on prend au maximum 400 tuiles pour limiter les
# temps de calcul du CNN
# --------------------------------------------------------
library("survival")
library("ggplot2")
library("survminer")
library('gtsummary')
library('FactoMineR')
library(maxstat)
require(dplyr)

# On charge les donnees cliniques
clinical <- read.csv("/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/clinical_data_HES_BSC.csv")
clinical$PFS <- as.numeric(as.vector(clinical$PFS))
clinical$PROG <- as.numeric(as.vector(clinical$PROG))
clinical$OS <- as.numeric(as.vector(clinical$OS))
clinical$DCD <- as.numeric(as.vector(clinical$DCD))

# REMOVE DFS NAN
#clinical<-clinical[-which(is.na(clinical[,'PFS'])),]

median(clinical$PFS)

clinical <- clinical[which(clinical$PFS>0),]
which(is.na(clinical$PROG))
fit.null <- survfit(Surv(PFS, PROG) ~ 1, data = clinical)
surv_median(fit.null)
# KAPLAN MEIR  / 
ggsurvplot(fit.null, pval=TRUE,risk.table = TRUE, surv.median.line = 'h',xlab ='PFS _ Months', ggtheme = theme_grey())

length(which(clinical$PFS<6&clinical$PROG==1)) # 53 patients progresseurs
length(which(clinical$PFS>36))  # 33 patients consideres comme repondeurs

# POur chaque patient on lui attribue sa classe : REP ou PROG + sa BD : train , val ou test
PROG_bcr <- as.vector(clinical$HES_x[which(clinical$PFS<6&clinical$PROG==1)])
REP_bcr <- as.vector(clinical$HES_x[which(clinical$PFS>36)])

train_PROG <- sample(PROG_bcr,3/5*length(PROG_bcr)) # Train contient 3/5 des donnees 
temp <- PROG_bcr[-which(PROG_bcr %in% train_PROG)]
val_PROG <- sample(temp,1/2*length(temp)) # Val contient 1/5 des donnees 
test_PROG <- temp[-which(temp %in% val_PROG)] # Test contient 1/5 des donnees 

train_REP <- sample(REP_bcr,3/5*length(REP_bcr)) # Train contient 3/5 des donnees 
temp <- REP_bcr[-which(REP_bcr %in% train_REP)]
val_REP <- sample(temp,1/2*length(temp)) # Val contient 1/5 des donnees 
test_REP <- temp[-which(temp %in% val_REP)] # Test contient 1/5 des donnees 

patients <- as.data.frame(c(train_PROG,val_PROG,test_PROG, train_REP, val_REP, test_REP))
patients$folder <- c(rep("Train",length(train_PROG)),rep("Val",length(val_PROG)),rep("Test",length(test_PROG)),
                   rep("Train",length(train_REP)),rep("Val",length(val_REP)),rep("Test",length(test_REP)))
patients$subfolder <- c(rep("PROG",length(train_PROG)+length(val_PROG)+length(test_PROG)),rep("REP",length(train_REP)+length(val_REP)+length(test_REP)))

colnames(patients) <- c("patient",'folder','subfolder')
patients$patient <- as.vector(patients$patient)

# enlever la lame 150 car n'est pas bien numerisé
patients=patients[-c(83),]

#write.csv(patients,"/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/prog_rep_patients.csv", row.names = FALSE)

prediction_file <- list.files("/work/shared/ptbc/CNN_Pancreas_V2/Donnees/project_kernel03_scan21000500/Results/SCNN/TumorVSNormal/")

# On distribue les tuiles des patients suivant la BD et la classe attribuees
for (patient in patients$patient)
{
  folder <- patients$folder[which(patients$patient %in% patient)]
  subfolder <- patients$subfolder[which(patients$patient %in% patient)]
  print(paste0("Je traite le patient: ",patient," qui sera dans le dossier ",folder,"/",subfolder))
  myFiles <- prediction_file[grep(patient,prediction_file)]
  # Si le patient a plus dune lame on en selectionne quune seule au hasard (a moduler)
  if(length(myFiles) > 1)
  {
    myFiles <- sample(myFiles,1)
  }
  for (myfile in myFiles)
  {
    pred <- read.csv(paste0("/work/shared/ptbc/CNN_Pancreas_V2/Donnees/project_kernel03_scan21000500/Results/SCNN/TumorVSNormal/",myfile), header = F)
    pred$pred_class <- ifelse(pred$V2>0.5,"Normal","Tumor")
    # On ne prendra que les tuiles qui sont predites a 95% tumeur pour ne travailler que sur la tumeur 
    tumor_tiles <- as.vector(pred$V1[which(pred$pred_class == "Tumor"&pred$V3>0.95)])
    a <- unlist(strsplit(as.vector(tumor_tiles), "[/]"))
    tumor_tiles <- a[seq(2,length(a), by = 2)]
    #slide_name <- gsub("myPrediction_","",myfile)
    slide_name <- gsub(".csv","",myfile)
    input_folder <- paste0('/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/N_Tselect_macenko/',slide_name)
    output_folder <- paste0('/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/BSC_DB_survie/',folder,'/',subfolder)
    # On limite enormement le nombre de tuiles pour chaque lame pour gagner du temps de calcul pour le CNN a moduler en fonction des resultats
    if(length(tumor_tiles) > 400) 
    {
      tumor_tiles <- sample(tumor_tiles,400)
    }
    for(tile in tumor_tiles)
    {
      file.copy(paste0(input_folder,"/",tile),
                paste0(output_folder,"/",tile))
    }
    print(paste0("Il aura ", length(tumor_tiles), " tuiles."))
  }
}


table(patients$folder, patients$subfolder)
































