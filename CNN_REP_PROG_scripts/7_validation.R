
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

# validation v6 color aug (entrainer sur macenko) sur les données de canceropol
set.seed(7)
densenet_color_aug <- '/work/shared/ptbc/CNN_Pancreas_V2/PROG_REP_results/CNN_REP_PROG_bscdata/color_augmentation/data_related_survival/bsc_mean_macenko_color_aug_2.csv'
mean_densenet_v6 <- '/work/shared/ptbc/CNN_Pancreas_V2/PROG_REP_results/CNN_REP_PROG_bscdata/data_related_survival/myTable_mean_macenko_2.csv'

myTable <- read.csv(densenet_color_aug)
myTable$patient <- as.vector(myTable$patient)
myTable$PFS <- ifelse(myTable$PFS>36,36,myTable$PFS)
myTable$PROG <- ifelse(myTable$PFS>36,0,myTable$PROG)
#SPLIT
train_patient <- sample(myTable$patient,3/4*length(myTable$patient))
myTable_train <- myTable[which(myTable$patient %in% train_patient),]
myTable_test <- myTable[-which(myTable$patient %in% train_patient),]

# VARS
#vars <- c('V64','V107','V172','V363','V548','V562','V606','V617','V643','V831','V857','V909','V953') #mean_densenet_v6
vars <- c('V179','V301','V304','V581','V628','V865','V929','V974')
# V179+V301+V304+V581+V628+V865+V929+V974
# V64+V107+V172+V363+V548+V562+V606+V617+V643+V831+V857+V909+V953
# TRAIN
modele_cont <- coxph(Surv(myTable_train$PFS,myTable_train$PROG)~ V179+V301+V304+V581+V628+V865+V929+V974, data= myTable_train)
myTable_train$pred <- modele_cont$linear.predictors
c <- median(myTable_train$pred, na.rm=T)
c <- maxstat.test(Surv(PFS,PROG)~pred,data=myTable_train, smethod="LogRank",minprop = 0.3, maxprop = 0.7,pmethod="condMC", B = 9999)$estimate
c <- quantile(myTable_train$pred, probs=1/3)

myTable_train$pred_dicho <- relevel(as.factor(ifelse(myTable_train$pred>c,"high","low")), ref='low')
km <- survfit(Surv(PFS,PROG)~pred_dicho,data=myTable_train)
a <- summary(coxph(Surv(PFS,PROG)~pred_dicho,data=myTable_train))

#TEST
myTable_test$pred <- predict(modele_cont, newdata = myTable_test)
myTable_test$pred_dicho <- relevel(as.factor(ifelse(myTable_test$pred>c,"high","low")), ref='low')
km <- survfit(Surv(PFS,PROG)~pred_dicho,data=myTable_test)
a <- summary(coxph(Surv(PFS,PROG)~pred_dicho,data=myTable_test))

#TCGA
tcga <- read.csv("/work/shared/ptbc/CNN_Pancreas_V2/results/TCGA/macenko_prog_rep/data_related_survival/densenet_julieDatagen.csv")
#tcga <- read.csv("/work/shared/ptbc/CNN_Pancreas_V2/results/TCGA/macenko_prog_rep/data_related_survival/mean_densenet_v6.csv")
tcga=tcga[complete.cases(tcga[,vars]), ]
tcga$PFI.time <- tcga$PFI.time/30.42
tcga$PFI <- ifelse(tcga$PFI.time>36,0,tcga$PFI)
tcga$PFI.time <- ifelse(tcga$PFI.time>36,36,tcga$PFI.time)
tcga$pred <- predict(modele_cont, newdata = tcga)

#c <- maxstat.test(Surv(PFI.time,PFI)~pred,data=tcga, smethod="LogRank",minprop = 0.3, maxprop = 0.7,pmethod="condMC", B = 9999)$estimate
#c <- quantile(tcga$pred, probs=1/3)
#c <- median(tcga$pred, na.rm=T)

tcga$pred_dicho <- relevel(as.factor(ifelse(tcga$pred>c,"high","low")), ref='low')
km <- survfit(Surv(PFI.time,PFI)~pred_dicho,data=tcga)
a <- summary(coxph(Surv(PFI.time,PFI)~pred_dicho,data=tcga))

#CANCEROPOL
cp <- read.csv('/work/shared/ptbc/CNN_Pancreas_V2/results/CP/macenko_prog_rep/data_related_survival/with_color_aug.csv')
#cp <- read.csv('/work/shared/ptbc/CNN_Pancreas_V2/results/CP/macenko_prog_rep/data_related_survival/without_color_aug.csv')
cp=cp[complete.cases(cp[,vars]), ]
cp$dfsev <- ifelse(cp$dfsev>36,0,cp$dfsev)
cp$dfs <- ifelse(cp$dfs>36,36,cp$dfs)
cp$pred <- predict(modele_cont, newdata = cp)

#c <- maxstat.test(Surv(dfs,dfsev)~pred,data=cp, smethod="LogRank",minprop = 0.3, maxprop = 0.7,pmethod="condMC", B = 9999)$estimate
#c <- quantile(cp$pred, probs=1/3)
#c <- median(cp$pred, na.rm=T)

cp$pred_dicho <- relevel(as.factor(ifelse(cp$pred>c,"high","low")), ref='low')
km <- survfit(Surv(dfs,dfsev)~pred_dicho,data=cp)
a <- summary(coxph(Surv(dfs,dfsev)~pred_dicho,data=cp))

####### PLOTING
ggsurv <- ggsurvplot(km, legend = "bottom",
                     legend.title = "",
                     palette = c("#072eee", "#ff0000"),
                     pval = TRUE,
                     risk.table = TRUE, 
                     risk.table.y.text.col = TRUE)

ggsurv$plot <- ggsurv$plot+
  ggplot2::annotate("text",
                    x = 26, y = 0.75, # x and y coordinates of the text
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


####### DISTRIBUTION DES PREDICTEUR LINEAIRE
breaks<- 30
hist(myTable_test$pred, col=rgb(0.87, 0.91, 0.070, 0.8),breaks=breaks, yaxt='n', ann=FALSE, axes = FALSE);par(new=TRUE)
hist(tcga$pred,col=rgb(0.12, 0.21, 0.6, 0.7),breaks=breaks,yaxt='n', ann=FALSE, axes = FALSE);par(new=TRUE)
hist(cp$pred,col=rgb(0.12, 0.6, 0.12, 0.6),breaks=breaks,yaxt='n', ann=FALSE, axes = FALSE);par(new=TRUE)
hist(myTable_test$pred,col=rgb(0.91, 0.21, 0.070, 0.5),breaks=breaks,yaxt='n', ann=FALSE, axes = FALSE);par(new=TRUE)
hist(myTable_train$pred, col=rgb(0.6, 0.12, 0.58, 0.9), xlim=c(-5,5),breaks=breaks);par(new=TRUE)
legend(2,10,legend=c('BSC_TRAIN','BSC_TEST','TCGA','CP','BSC_TEST'),col=c('pink','yellow','blue','green','red'), lwd = 2, cex=0.8)


############ tester les variables entres les cohortes
box = c()

df1 <- myTable_train[,vars]
df1$data <- 'BSC_TRAIN'

df2 <- myTable_test[,vars]
df2$data <- 'BSC_TEST'

df3 <- tcga[,vars]
df3$data <- 'TCGA'

df4 <- cp[,vars]
df4$data <- 'CANCEROPOL'

box <- rbind(box,df1, df2, df3, df4)



for (var in vars){
  print(paste('----------------------',var))
  box$temp <- box[,var]
  #pair <- pairwise.wilcox.test(box$temp, box$data, p.adjust.method='BH')
  pair <- kruskal.test(box$temp, box$data, p.adjust.method='BH')
  print(pair$p.value)
}

library(reshape2)
box.m <- melt(box,id.vars='data', measure.vars=vars)
p <- ggplot(box.m) +
  geom_boxplot(aes(x=data, y=value, color=variable)) +
  theme_bw()

p


























































