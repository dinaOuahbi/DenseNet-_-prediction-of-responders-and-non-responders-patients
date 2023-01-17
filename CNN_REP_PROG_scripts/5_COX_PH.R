

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
library(survivalROC)


# On recharge une des trois bases
set.seed(7)
mean_densenet_v6 <- '/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/data_related_survival/myTable_mean_macenko_2.csv'
densenet_color_aug <- '/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/color_augmentation/data_related_survival/bsc_mean_macenko_color_aug_2.csv'
#densenet_v6_julieDatagen <- '/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/color_augmentation/data_related_survival/bsc_mean_macenko_color_aug_comb_datagen_2.csv'

myTable <- read.csv(densenet_julieDatagen)
myTable$patient <- as.vector(myTable$patient)
myTable$PFS <- ifelse(myTable$PFS>36,36,myTable$PFS)
myTable$PROG <- ifelse(myTable$PFS>36,0,myTable$PROG)
#SPLIT
train_patient <- sample(myTable$patient,3/4*length(myTable$patient))
myTable_train <- myTable[which(myTable$patient %in% train_patient),]
myTable_test <- myTable[-which(myTable$patient %in% train_patient),]



# ACP 
first_cnn_var_ind <- grep("^V2$", colnames(myTable_train))
last_cnn_var_ind <- grep("^V1025$", colnames(myTable_train))

res.pca <- PCA(myTable_train[first_cnn_var_ind:last_cnn_var_ind], graph = FALSE)
hc <- HCPC(res.pca, nb.clust=2, graph = FALSE)
fviz_cluster(hc,repel = TRUE,show.clust.cent = TRUE,ggtheme = theme_minimal(),palette = c("#072eee", "#ff0000"))
fviz_pca_var(res.pca, col.var="contrib")+
  scale_color_gradient2(low="black", mid="gray",high="yellow", midpoint=0.1, space ="Lab") +
  theme_minimal()

myTable_train$GROUPE <- hc$data.clust$clust


#limiter la survie a 2ans et demi
#myTable_train$PROG <- ifelse(myTable_train$PFS>30,0,myTable_train$PROG)
#myTable_train$PFS <- ifelse(myTable_train$PFI.time>30,30,myTable_train$dfs)

km <- survfit(Surv(PFS,PROG)~GROUPE,data=myTable_train)
ggsurvplot(km, legend = "bottom",
           legend.title = "HC _ All_CNN_Variables",
           palette = c("#072eee", "#ff0000"),pval = TRUE,risk.table = TRUE, 
           risk.table.y.text.col = TRUE,surv.median.line = "h")

# SLECT VAR WITH HIGHT VARIANCE
rm_var <- names(which(apply(myTable_train[,first_cnn_var_ind:last_cnn_var_ind],2,function(x){var(x,na.rm = T)})<0.0000001))
myTable_train_new <- myTable_train[,-which(colnames(myTable_train) %in% rm_var)]

fit.null <- survfit(Surv(PFS,PROG)~GROUPE,data=myTable_train)
surv_median(fit.null)
ggsurvplot(fit.null, legend = "bottom",
           legend.title = "unsupervised model _ All_CNN_Variables",
           palette = c("#072eee", "#ff0000"),pval = TRUE,risk.table = TRUE, 
           risk.table.y.text.col = TRUE,surv.median.line = "h")


# modele de cox univarié
Cox_univ <- Cox_univ_dicho <-  c()
pfs_index_col <- grep("^PFS$", colnames(myTable_train))
prog_index_col <- grep("^PROG$", colnames(myTable_train))

first_cnn_var_ind <- grep("^V4$", colnames(myTable_train_new))
last_cnn_var_ind <- grep("^V1025$", colnames(myTable_train_new))

for (i in first_cnn_var_ind:last_cnn_var_ind){
  print(i)
  temp <- myTable_train_new[,c(prog_index_col,pfs_index_col,i)]
  a <- summary(coxph(Surv(temp$PFS,temp$PROG)~temp[,3],data=temp))
  Cox_univ <- rbind(Cox_univ,cbind(a$conf.int[,1],a$conf.int[,3],a$conf.int[,4],a$coefficients[,5])) # HR / LOW / HIGHT / PVAL
  c <- maxstat.test(Surv(temp$PFS,temp$PROG)~temp[,3],data=temp, smethod="LogRank",minprop = 0.3, maxprop = 0.7,pmethod="condMC", B = 9999)$estimate
  temp$DICHO <- relevel(as.factor(ifelse(temp[,3]>c,"high","low") ), ref='low')
  a <- summary(coxph(Surv(temp$PFS,temp$PROG)~DICHO,data=temp))
  Cox_univ_dicho <- rbind(Cox_univ_dicho,cbind(a$conf.int[,1],a$conf.int[,3],a$conf.int[,4],a$coefficients[,5]))

}

Cox_univ <- as.data.frame(Cox_univ)
Cox_univ_dicho <- as.data.frame(Cox_univ_dicho)
rownames(Cox_univ) <- rownames(Cox_univ_dicho) <- colnames(myTable_train_new)[first_cnn_var_ind:last_cnn_var_ind]
colnames(Cox_univ) <- colnames(Cox_univ_dicho) <- c("HR","lower","upper","Cox_univ")

select <- rownames(Cox_univ)[which(Cox_univ$Cox_univ<0.01)]

# ACP
res.pca <- PCA(myTable_train_new[,select], graph = FALSE) 
hc <- HCPC(res.pca, nb.clust=2, graph=FALSE) 
fviz_cluster(hc,repel = TRUE,palette = c("#072eee", "#ff0000"),show.clust.cent = TRUE,ggtheme = theme_minimal()) #LES CLUSTERS
fviz_pca_var(res.pca, col.var="contrib")+
  scale_color_gradient2(low="black", mid="green",high="red", midpoint=10, space ="Lab") +
  theme_minimal()
# KM
myTable_train_new$GROUPE <- hc$data.clust$clust
km <- survfit(Surv(PFS,PROG)~GROUPE,data=myTable_train_new)
ggsurvplot(km, legend = "bottom",
           legend.title = "Variable significatives en univarié",
           palette = c("#072eee", "#ff0000"),pval = TRUE,risk.table = TRUE, 
           risk.table.y.text.col = TRUE,surv.median.line = "hv")



myTable_train_new$PROG <- as.numeric(as.vector(myTable_train_new$PROG))
myTable_train_new$PFS <- as.numeric(as.vector(myTable_train_new$PFS))

# LASSO
x <- as.matrix(data.frame(myTable_train[,select]))
y <- Surv(myTable_train$PFS,myTable_train$PROG)
myResult <- as.data.frame(colnames(x))
myResult$app <- rep(0,length(select))


for (i in 1:100){
  print(i)
  cv.fit <- cv.glmnet(x,y,alpha = 1,family = "cox")
  #plot(cv.fit)
  a <- coef(cv.fit, s = "lambda.min")
  myResult$app[which(coef(cv.fit, s = "lambda.min")!=0)] <- myResult$app[which(coef(cv.fit, s = "lambda.min")!=0)]+1
}

var_select <- as.vector(myResult$`colnames(x)`[which(myResult$app>50)])
res.pca <- PCA(myTable_train[,var_select], graph = FALSE)
fviz_pca_var(res.pca, col.var="contrib")+
  scale_color_gradient2(low="white", mid="blue",high="red", midpoint=7, space ="Lab") +
  theme_minimal()
hc <- HCPC(res.pca, nb.clust=2, graph = FALSE)
fviz_cluster(hc,repel = TRUE,palette = c("#00AFBB", "#E7B800"),show.clust.cent = TRUE,ggtheme = theme_minimal()) #LES CLUSTERS

myTable_train$GROUPE <- hc$data.clust$clust
km <- survfit(Surv(PFS,PROG)~GROUPE,data=myTable_train)
ggsurv <- ggsurvplot(km, legend = "bottom",
                     legend.title = "LASSO COX : variable à coeff de regression > 50 (sur 100 iterations)",
                     palette = c("#00AFBB", "#E7B800"),pval = TRUE,risk.table = TRUE, 
                     risk.table.y.text.col = TRUE,linetype = "strata",surv.median.line = "hv")



# MULTIVARIÉ
modele_cont <- coxph(Surv(myTable_train$PFS,myTable_train$PROG)~ V179+V301+V304+V581+V628+V865+V929+V974, data= myTable_train)
myTable_train$pred <- modele_cont$linear.predictors

# 3 seuils
#c <- median(myTable_train$pred, na.rm=T)
#c <- quantile(myTable_train$pred, probs=1/3)
c <- maxstat.test(Surv(PFS,PROG)~pred,data=myTable_train, smethod="LogRank",minprop = 0.3, maxprop = 0.7,pmethod="condMC", B = 9999)$estimate

myTable_train$pred_dicho <- relevel(as.factor(ifelse(myTable_train$pred>c,"high","low")), ref='low')
km <- survfit(Surv(PFS,PROG)~pred_dicho,data=myTable_train)
a <- summary(coxph(Surv(PFS,PROG)~pred_dicho,data=myTable_train))

ggsurv <- ggsurvplot(km, legend = "bottom",
                     legend.title = "Train",
                     palette = c("#072eee", "#ff0000"),
                     pval = TRUE,
                     risk.table = TRUE, 
                     risk.table.y.text.col = TRUE)
  
ggsurv$plot <- ggsurv$plot+
  ggplot2::annotate("text",
                    x = 25, y = 0.5, # x and y coordinates of the text
                    label = paste0('HR : ',
                                   round(a$coefficients[2],3),
                                   '[',
                                   round(a$conf.int[,3],3),
                                   ',',
                                   round(a$conf.int[,4],3),
                                   '] P : ',
                                   round(a$coefficients[5], 5)),
                    size = 5)
ggsurv
a

# SUR LE TEST
myTable_test$pred <- predict(modele_cont, newdata = myTable_test)
#seuil
#c <- median(myTable_test$pred, na.rm=T)
#c <- quantile(myTable_test$pred, probs=1/3)
#c <- maxstat.test(Surv(PFS,PROG)~pred,data=myTable_test, smethod="LogRank",minprop = 0.3, maxprop = 0.7,pmethod="condMC", B = 9999)$estimate

myTable_test$pred_dicho <- relevel(as.factor(ifelse(myTable_test$pred>c,"high","low")), ref='low')
km <- survfit(Surv(PFS,PROG)~pred_dicho,data=myTable_test)
a <- summary(coxph(Surv(PFS,PROG)~pred_dicho,data=myTable_test))
ggsurv <- ggsurvplot(km, legend = "bottom",
                     legend.title = "Test",
                     palette = c("#072eee", "#ff0000"),
                     pval = TRUE,risk.table = TRUE, 
                     risk.table.y.text.col = TRUE)

ggsurv$plot <- ggsurv$plot+
  ggplot2::annotate("text",
                    x = 25, y = 0.5, # x and y coordinates of the text
                    label = paste0('HR : ',
                                   round(a$coefficients[2],3),
                                   '[',
                                   round(a$conf.int[,3],3),
                                   ',',
                                   round(a$conf.int[,4],3),
                                   '] P : ',
                                   round(a$coefficients[5],3)),
                    size = 5)
ggsurv
a

################# AJOUT DES VARIABLES CLINIQUES 
clinique <- read.csv("/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/clinical_data_HES_BSC_v2.csv")
clinique <- clinique[,c("HES_x","Ctadj","Histo","PN8","pN")]
colnames(clinique)[1] <- "patient"
clin_train <- merge(clinique,myTable_train)
clin_test <- merge(clinique,myTable_test)



modele_cont <- coxph(Surv(clin_train$PFS,clin_train$PROG)~ V35+V60+V420+V518+V557+V564+V595+V640+V661+V701+V765+Ctadj+Histo+PN8+pN,data= clin_train)
clin_train=clin_train[complete.cases(clin_train[,c('V35','V60','V420','V518','V557','V564','V595','V640','V661','V701','V765','Ctadj','Histo','PN8','pN')]), ] 
clin_train$pred <- modele_cont$linear.predictors

c <- median(clin_train$pred, na.rm=T)
c <- quantile(clin_train$pred, probs=1/3)
c <- maxstat.test(Surv(PFS,PROG)~pred,data=clin_train, smethod="LogRank",minprop = 0.3, maxprop = 0.7,pmethod="condMC", B = 9999)$estimate

clin_test$pred <- predict(modele_cont, newdata = clin_test)
clin_test$pred_dicho <- relevel(as.factor(ifelse(clin_test$pred>c,"high","low")), ref='low')
km <- survfit(Surv(PFS,PROG)~pred_dicho,data=clin_test)
a <- summary(coxph(Surv(PFS,PROG)~pred_dicho,data=clin_test))
ggsurvplot(km, legend = "bottom",legend.title = "CNN+clinical data",palette = c("#00AFBB", "#E7B800"),pval = TRUE,risk.table = TRUE,risk.table.y.text.col = TRUE,linetype = "strata")
a

#========================> ANOVA (Train puis test)
cnn_model <- coxph(Surv(clin_test$PFS,clin_test$PROG)~ V35+V60+V420+V518+V557+V564+V595+V640+V661+V701+V765, data= clin_test)
clin_model <- coxph(Surv(clin_test$PFS,clin_test$PROG)~ Ctadj+Histo+PN8+pN, data= clin_test)
merge_model <- coxph(Surv(clin_test$PFS,clin_test$PROG)~ V35+V60+V420+V518+V557+V564+V595+V640+V661+V701+V765+Ctadj+Histo+PN8+pN,data= clin_test)

anova(merge_model,cnn_model)
anova(merge_model,clin_model)

cutoff = 24
cnn_roc= survivalROC(Stime=clin_test$PFS,status=clin_test$PROG,marker = cnn_model$linear.predictors,predict.time = cutoff,span = 0.25*NROW(clin_test)^(-0.20))
clin_roc= survivalROC(Stime=clin_test$PFS,status=clin_test$PROG,marker = clin_model$linear.predictors,predict.time = cutoff,span = 0.25*NROW(clin_test)^(-0.20))
merge_roc= survivalROC(Stime=clin_test$PFS,status=clin_test$PROG,marker = merge_model$linear.predictors,predict.time = cutoff,span = 0.25*NROW(clin_test)^(-0.20))


AUC <- data.frame(name=c('COX ph CNN','COX ph Clinical', 'COX ph Merge') ,value=c(cnn_roc$AUC, clin_roc$AUC, merge_roc$AUC))
ggplot(AUC, aes(x=name, y=value)) + 
  geom_bar(stat = "identity", width=0.2)



#========================> ROC
plot(cnn_roc$FP, cnn_roc$TP,type="l", lty = 5, col='red');par(new=TRUE)
plot(clin_roc$FP, clin_roc$TP,type="l",lty = 5, col='blue',ann=FALSE);par(new=TRUE)
plot(merge_roc$FP, merge_roc$TP,type="l",lty = 5, col='green',ann=FALSE);par(new=TRUE)
legend(0.55,0.2,legend=c('COX ph CNN','COX ph Clinical', 'COX ph Merge'),col=c("red", "blue",'green'), lty=1:2, cex=0.8)
abline(0,1, type='l', lty=2)


####################################### INITIALIZE MY COX PH MODEL

modele_cont <- coxph(Surv(myTable_train$PFS,myTable_train$PROG)~ V179+V301+V304+V581+V628+V865+V929+V974, data= myTable_train)
myTable_train$pred <- modele_cont$linear.predictors
c <- maxstat.test(Surv(PFS,PROG)~pred,data=myTable_train, smethod="LogRank",minprop = 0.3, maxprop = 0.7,pmethod="condMC", B = 9999)$estimate
c <- quantile(myTable_train$pred, probs=1/3)
c <- median(myTable_train$pred, na.rm=T)

myTable_train$pred_dicho <- relevel(as.factor(ifelse(myTable_train$pred>c,"high","low")), ref='low')
km <- survfit(Surv(PFS,PROG)~pred_dicho,data=myTable_train)
a <- summary(coxph(Surv(PFS,PROG)~pred_dicho,data=myTable_train))

ggsurv <- ggsurvplot(km, legend = "bottom",
                     legend.title = "Train",
                     palette = c("#072eee", "#ff0000"),
                     pval = TRUE,
                     risk.table = TRUE, 
                     risk.table.y.text.col = TRUE)

ggsurv$plot <- ggsurv$plot+
  ggplot2::annotate("text",
                    x = 25, y = 0.5, # x and y coordinates of the text
                    label = paste0('HR : ',
                                   round(a$coefficients[2],3),
                                   '[',
                                   round(a$conf.int[,3],3),
                                   ',',
                                   round(a$conf.int[,4],3),
                                   '] P : ',
                                   round(a$coefficients[5], 5)),
                    size = 5)
ggsurv
a


###################################### test notre model cox ph CNN sur les données du TCGA
tcga <- read.csv("/work/shared/ptbc/CNN_Pancreas_V2/results/TCGA/macenko_prog_rep/data_related_survival/mean_densenet_v6.csv")
tcga=tcga[complete.cases(tcga[,c('V64','V107','V172','V363','V548','V562','V606','V617','V643','V831','V857','V909','V953')]), ]

tcga <- read.csv("/work/shared/ptbc/CNN_Pancreas_V2/results/TCGA/macenko_prog_rep/data_related_survival/densenet_julieDatagen.csv")
tcga=tcga[complete.cases(tcga[,c('V179','V301','V304','V581','V628','V865','V929','V974')]), ]

tcga <- read.csv("/work/shared/ptbc/CNN_Pancreas_V2/results/TCGA/macenko_prog_rep/data_related_survival/densenet_v6_julieDatagen.csv")
tcga=tcga[complete.cases(tcga[,c('V533','V566','V638','V805','V834','V901','V964','V970')]), ]

## SOME PROCESSING
tcga$PFI.time <- tcga$PFI.time/30.42
tcga$PFI <- ifelse(tcga$PFI.time>36,0,tcga$PFI)
tcga$PFI.time <- ifelse(tcga$PFI.time>36,36,tcga$PFI.time)
tcga$pred <- predict(modele_cont, newdata = tcga)

c <- maxstat.test(Surv(PFI.time,PFI)~pred,data=tcga, smethod="LogRank",minprop = 0.3, maxprop = 0.7,pmethod="condMC", B = 9999)$estimate
c <- quantile(tcga$pred, probs=1/3)
c <- median(tcga$pred, na.rm=T)

tcga$pred_dicho <- relevel(as.factor(ifelse(tcga$pred>c,"high","low")), ref='low')
km <- survfit(Surv(PFI.time,PFI)~pred_dicho,data=tcga)
a <- summary(coxph(Surv(PFI.time,PFI)~pred_dicho,data=tcga))
ggsurv <- ggsurvplot(km, legend = "bottom",
                     legend.title = "TCGA MACENKO",
                     palette = c("#072eee", "#ff0000"),
                     pval = TRUE,
                     risk.table = TRUE, 
                     risk.table.y.text.col = TRUE)

ggsurv$plot <- ggsurv$plot+
  ggplot2::annotate("text",
                    x = 25, y = 0.5, # x and y coordinates of the text
                    label = paste0('HR : ',
                                   round(a$coefficients[2],3),
                                   '[',
                                   round(a$conf.int[,3],3),
                                   ',',
                                   round(a$conf.int[,4],3),
                                   '] P : ',
                                   round(a$coefficients[5], 5)),
                    size = 5)
ggsurv


########################## sur canceropoel
cp <- read.csv("/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/data_related_survival/myTable_mean_cp_2.csv")
cp=cp[complete.cases(cp[,c('V35','V60','V420','V518','V557','V564','V595','V640','V661','V701','V765')]), ]
cp$dfsev <- as.numeric(as.vector(cp$dfsev))
cp$pred <- predict(modele_cont, newdata = cp)

c <- maxstat.test(Surv(dfs,dfsev)~pred,data=cp, smethod="LogRank",minprop = 0.3, maxprop = 0.7,pmethod="condMC", B = 9999)$estimate
c <- quantile(cp$pred, probs=1/3)
c <- median(cp$pred, na.rm=T)

cp$pred_dicho <- relevel(as.factor(ifelse(cp$pred>c,"high","low")), ref='low')
km <- survfit(Surv(dfs,dfsev)~pred_dicho,data=cp)
a <- summary(coxph(Surv(dfs,dfsev)~pred_dicho,data=cp))
ggsurvplot(km, legend = "bottom",legend.title = "Linear predictor",palette = c("#00AFBB", "#E7B800"),pval = TRUE,risk.table = TRUE,risk.table.y.text.col = TRUE,linetype = "strata")
a



############################### dist des predicteur lineaires
hist(myTable_train$pred, col=rgb(0, 0, 1, 0.5), breaks=30, xlim=c(-2,2));par(new=TRUE)
hist(myTable_test$pred,col=rgb(0, 1, 0, 0.5),breaks=30,xlim=c(-2,2),yaxt='n', ann=FALSE);par(new=TRUE)
#hist(cp$pred,col='red',breaks=20,xlim=c(-12,5),yaxt='n', ann=FALSE);par(new=TRUE)

legend(1,5,legend=c('Train','Test'),col=c('blue','green'), lty=1:2, cex=0.8)

abline(v=median, col='red', lwd=3, lty='dashed')
abline(v=maxstat, col='blue', lwd=3, lty='dashed')
abline(v=quantile, col='orange', lwd=3, lty='dashed')


legend(1,7,legend=c('Median','Maxstat', 'Quantile'),col=c('red','blue', 'orange'), lty=1:2, cex=0.8)

############### BOXPLOT
box = c()
#df1 <- cp[,c('V35','V60','V420','V518','V557','V564','V595','V640','V661','V701','V765')]
#df1$data <- 'CP'

df2 <- tcga[,c('V64','V107','V172','V363','V548','V562','V606','V617','V643','V831','V857','V909','V953')]
df2$data <- 'TCGA'

df3 <- myTable[,c('V64','V107','V172','V363','V548','V562','V606','V617','V643','V831','V857','V909','V953')]
df3$data <- 'BSC'


box <- rbind(box,df1, df2, df3)
#Box plots basiques
#p <- ggplot(box, aes(x=data, y=V35, color=data)) + 
#  geom_boxplot(notch=TRUE,outlier.colour="red", outlier.shape=8,outlier.size=4)
#p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
#  scale_color_brewer(palette="Dark2")+
#  geom_jitter(shape=16, position=position_jitter(0.2))+
#  theme(legend.position="top")


pair <- pairwise.wilcox.test(box$V64, box$data, p.adjust.method='BH')
pair

for (var in c('V64','V107','V172','V363','V548','V562','V606','V617','V643','V831','V857','V909','V953')){
  print(paste('----------------------',var))
  box$temp <- box[,var]
  pair <- pairwise.wilcox.test(box$temp, box$data, p.adjust.method='BH')
  print(pair$p.value)
}






library(reshape2)
box.m <- melt(box,id.vars='data', measure.vars=c('V64','V107','V172','V363','V548','V562','V606','V617','V643','V831','V857','V909','V953'))
p <- ggplot(box.m) +
  geom_boxplot(aes(x=data, y=value, color=variable))
  
p


write.csv(box, '/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/data_related_survival/macenko_survival_signifCnnVar.csv')





















