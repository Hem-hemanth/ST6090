library(survival)
library(plyr)
library(randomForestSRC)
library(dplyr)
library(glmnet)
library(caret)
library(ggplot2)
library(survminer)
library(plotrix)
library(gbm)

#data preprocessing for clinical variables

#importing dataset
clinical_dataset <- read.csv('D:\\ST6090\\Data_project.csv',
                      header = T)
#filtering dataset
Filtered_data <- clinical_dataset %>% filter(!is.na(BCR_FreeTime) & 
                                   !is.na(BCR_Event))
#set sample id t row names
rownames(Filtered_data) <- Filtered_data[, 1]
Filtered_data <- Filtered_data[, -1]

#import rna
genomic_data <- read.table("D:\\ST6090\\MSKCC_PCa_mRNA_data.txt", header = T)
genomic_data <- t(genomic_data)
genomic_data <- genomic_data[-1,]
colnames(genomic_data) <- genomic_data[1, ]
genomic_data <- genomic_data[-1, ]
dim(genomic_data)

na.omit(genomic_data)
#combine clinical with rna
merged_data <- merge(Filtered_data, genomic_data, by = 0)
rownames(merged_data ) <- merged_data [, 1]
merged_data <- merged_data [, -1]

dim(merged_data)
#mrna variables are off character format, convert into numeric

merged_data[,39:26485] <- lapply(merged_data[,39:26485], function(x) {
  if(is.character(x)) as.numeric(x)
})

temp_data <- merged_data

#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#........................................DATA PREPROCESSING............................................#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#

merged_data <- subset(merged_data,
                    select = -c(2,20,21,28:35))
## Change BCR_Event, SMS, ECE and SVI to numeric
merged_data$BCR_Event <- ifelse(merged_data$BCR_Event == "NO", 0, 1)
merged_data$SMS.x <- ifelse(merged_data$SMS.x == "Negative", 0, 1)
merged_data$ECE <- ifelse(merged_data$ECE == "NONE", 0, 1)
merged_data$SVI <- ifelse(merged_data$SVI == "Negative", 0, 1)

## Change Staging to factors
merged_data$BxGG1 <- as.factor(merged_data$BxGG1)
merged_data$BxGG2 <- as.factor(merged_data$BxGG2)
merged_data$BxGGS <- as.factor(merged_data$BxGGS)
merged_data$ClinT_Stage <- as.factor(merged_data$ClinT_Stage)
merged_data$PathStage <- as.factor(merged_data$PathStage)
merged_data$PathGG1 <- as.factor(merged_data$PathGG1)
merged_data$PathGG2 <- as.factor(merged_data$PathGG2)

###Removed observation with PAthGGS = Tx Effects
merged_data <- merged_data[!(row.names(merged_data) %in% 'PCA0046'), ]
merged_data$PathGGS <- as.factor(merged_data$PathGGS)

## Change LNI to factor
merged_data$LNI <- as.factor(merged_data$LNI)
## Create a therapy column
merged_data$therapy <- ifelse(!is.na(merged_data$NeoAdjRadTx)|
                              !is.na(merged_data$ChemoTx)|
                              !is.na(merged_data$HormTx)|
                              !is.na(merged_data$RadTxType),1,0)
## Remove original therapy columns
merged_data <- subset(merged_data, select = -c(10:13))

#dim now
dim(merged_data)
#140 26467

## Remove rows with NA values
final_data <- na.omit(merged_data)

dim(final_data)

final_data$Race <- mapvalues(final_data$Race, from = c("Black Non Hispanic", "White Non Hispanic","Black Hispanic", "White Hispanic", "Asian","Unknown"),
                           to =c("2Black", "1White", "2Black", "1White", "3Asian","4Unknown"))
final_data$Race <- as.factor(final_data$Race)

## Change Type to Categorical
final_data$Type <- as.factor(final_data$Type)

## Drop Observation with PSA 506 (error)
final_data <- final_data[!(row.names(final_data) %in% 
                         c('PCA0045','PCA0187')), ]
#dim
dim(final_data)
#134 226467

##Change levels of Biopsy Gleason score lowest == 3, and sum == 6 
final_data$BxGG2 <- mapvalues(final_data$BxGG2, from = c(2, 3, 4, 5),
                            to = c(3, 3, 4, 5))
final_data$BxGGS <- mapvalues(final_data$BxGGS, from = c(5, 6, 7, 8, 9),
                            to = c(6, 6, 7, 8, 9))
final_data$newBGG1 <- mapvalues(final_data$BxGG1,
                              from = c(3, 4, 5), 
                              to = c("3", "4+", "4+"))
final_data$newPGG1 <- mapvalues(final_data$PathGG1,
                              from = c(3, 4, 5), 
                              to = c("3", "4+", "4+"))
final_data$newClinStage <- mapvalues(final_data$ClinT_Stage,
                                   from = c("T1C", "T2", "T2A",
                                            "T2B", "T2C", "T3",
                                            "T3A", "T3B", "T3C"), 
                                   to = c("T1", "T2", "T2",
                                          "T2", "T2", "T3",
                                          "T3", "T3", "T3"))
final_data$newPathStage <- mapvalues(final_data$PathStage,
                                   from = c("T2A","T2B", "T2C",
                                            "T3A", "T3B", "T3C",
                                            "T4"), 
                                   to = c("T2", "T2", "T2",
                                          "T3", "T3", "T3",
                                          "T4"))





final_data$SMS.x <- as.factor(final_data$SMS.x)
final_data$ECE <- as.factor(final_data$ECE)
final_data$SVI <- as.factor(final_data$SVI)
final_data$Copy.number.Cluster <- as.numeric(final_data$Copy.number.Cluster)
final_data$ERG.fusion.aCGH <- as.factor(final_data$ERG.fusion.aCGH)
final_data$ERG.fusion.gex <- as.factor(final_data$ERG.fusion.gex)

#log rank
Y = Surv(final_data$BCR_FreeTime,final_data$BCR_Event==1)

pvallist <- numeric(21)

l1 <- survdiff(Y~final_data$Type)
pvallist[1] <- 2e-15
l2 <- survdiff(Y~final_data$Race)
pvallist[2] <- 0.6
l3 <- survdiff(Y~final_data$PreDxBxPSA)
pvallist[3] <- 2e-16
l4 <- survdiff(Y~final_data$DxAge)
pvallist[4] <- 2e-16
l5 <- survdiff(Y~final_data$newBGG1)
pvallist[5] <- 2e-08
l6 <- survdiff(Y~final_data$BxGG2)
pvallist[6] <- 0.001
l7 <- survdiff(Y~final_data$BxGGS)
pvallist[7] <- 5e-05
l8 <- survdiff(Y~final_data$PreTxPSA)
pvallist[8] <- 2e-16
l9 <- survdiff(Y~final_data$newClinStage)
pvallist[9] <- 0.002
l10 <- survdiff(Y~final_data$RP_Type)
pvallist[10] <- 0.4
l11 <- survdiff(Y~final_data$SMS.x)
l11
pvallist[11] <- 0.05
l12 <- survdiff(Y~final_data$ECE)
pvallist[12] <- 0.08
l13 <- survdiff(Y~final_data$SVI)
l13
pvallist[13] <- 1e-08
l14<- survdiff(Y~final_data$LNI)
l14
pvallist[14] <- 2e-16
l15 <- survdiff(Y~final_data$newPathStage)
pvallist[15] <- 2e-06
l16 <- survdiff(Y~final_data$newPGG1)
pvallist[16] <- 1e-10
l17 <- survdiff(Y~final_data$PathGG2)
pvallist[17] <- 1e-12
l18 <- survdiff(Y~final_data$PathGGS)
pvallist[18] <- 2e-16
l19 <- survdiff(Y~final_data$Copy.number.Cluster)
pvallist[19] <- 0.007
l20 <- survdiff(Y~final_data$ERG.fusion.aCGH)
pvallist[20] <- 0.2
l21 <- survdiff(Y~final_data$ERG.fusion.gex)
pvallist[20] <- 0.6
l22 <- survdiff(Y~final_data$therapy)
pvallist[21] <- 2e-16

final_data$newBGG1 <- mapvalues(final_data$BxGG1,
                              from = c(3, 4, 5), 
                              to = c("3", "4+", "4+"))
final_data$newPGG1 <- mapvalues(final_data$PathGG1,
                              from = c(3, 4, 5), 
                              to = c("3", "4+", "4+"))
final_data$newClinStage <- mapvalues(final_data$ClinT_Stage,
                                   from = c("T1C", "T2", "T2A",
                                            "T2B", "T2C", "T3",
                                            "T3A", "T3B", "T3C"), 
                                   to = c("T1", "T2", "T2",
                                          "T2", "T2", "T3",
                                          "T3", "T3", "T3"))
final_data$newPathStage <- mapvalues(final_data$PathStage,
                                   from = c("T2A","T2B", "T2C",
                                            "T3A", "T3B", "T3C",
                                            "T4"), 
                                   to = c("T2", "T2", "T2",
                                          "T3", "T3", "T3"))
cbind(pvallist)

cbind(p.adjust(pvallist[1:21],method = 'fdr'))



#univariate cox test
coxpval = numeric(21)
Y <- Surv(final_data$BCR_FreeTime,final_data$BCR_Event)


l1 <- coxph(Y~final_data$Type)
coxpval[1] <- 1.37e-05
l2 <- coxph(Y~final_data$Race)
coxpval[2] <- 0.5025
l3 <- coxph(Y~final_data$PreDxBxPSA)
coxpval[3] <- 0.01427
l4 <- coxph(Y~final_data$DxAge)
coxpval[4] <- 0.3365
l5 <- coxph(Y~final_data$newBGG1)
coxpval[5] <- 3.422e-06
l6 <- coxph(Y~final_data$BxGG2)
coxpval[6] <- 0.02812
l7 <- coxph(Y~final_data$BxGGS)
coxpval[7] <- 0.001068
l8 <- coxph(Y~final_data$PreTxPSA)
coxpval[8] <- 0.0003364
l9 <- coxph(Y~final_data$newClinStage)
coxpval[9] <- 0.02864
l10 <- coxph(Y~final_data$RP_Type)
coxpval[10] <- 0.3571
l11 <- coxph(Y~final_data$SMS.x)
l11
coxpval[11] <- 0.06238
l12 <- coxph(Y~final_data$ECE)
coxpval[12] <- 0.06459
l13 <- coxph(Y~final_data$SVI)
l13
coxpval[13] <- 9.538e-06
l14<- coxph(Y~final_data$LNI)
l14
coxpval[14] <- 2.542e-07
l15 <- coxph(Y~final_data$newPathStage)
coxpval[15] <- 1.866e-05
l16 <- coxph(Y~final_data$newPGG1)
coxpval[16] <- 9.539e-09
l17 <- coxph(Y~final_data$PathGG2)
coxpval[17] <- 1.048e-05
l18 <- coxph(Y~final_data$PathGGS)
coxpval[18] <- 1.691e-09
l19 <- coxph(Y~final_data$Copy.number.Cluster)
coxpval[19] <- 0.4795
l20 <- coxph(Y~final_data$ERG.fusion.aCGH)
coxpval[20] <- 0.2201
l21 <- coxph(Y~final_data$ERG.fusion.gex)
coxpval[20] <- 0.6464
l22 <- coxph(Y~final_data$therapy)
coxpval[21] <- 1.027e-15
cbind(coxpval,p.adjust(coxpval,method = 'fdr'))

final_data <- subset(final_data,
                   select = -c(BxGG1, PathGG1, ClinT_Stage, 
                               PathStage))
#dim
dim(final_data)
#134 26467

## Remove therapy
final_data <- subset(final_data, select = -therapy)

final_data$Copy.number.Cluster<-as.factor(final_data$Copy.number.Cluster)
final_data$ERG.fusion.aCGH<-as.factor(final_data$ERG.fusion.aCGH)
final_data$ERG.fusion.gex <- as.factor(final_data$ERG.fusion.gex)
final_data$RP_Type=as.factor(final_data$RP_Type)

dim(final_data)


#data analysis done.

#correlation matrix

data <- final_data
library(dplyr)
cutoff <- data[,c(15, 16)]
cutoff <- cutoff %>% filter((BCR_FreeTime >= 24 & BCR_Event == 0) 
                            |
                              (BCR_FreeTime <= 24 & BCR_Event == 
                                 1))
feature_selection_data <- merge(cutoff, data, by = 0)
rownames(feature_selection_data) <- feature_selection_data[, 1]
feature_selection_data <- feature_selection_data[, -c(1:3)]

#Correlation between BCR and mRNA variables
#from 26463 is encoded variables
correlation_data <- cor(feature_selection_data[20:26466], feature_selection_data[15])
correlation_data <- data.frame(correlation_data)
colnames(correlation_data) <- "Cor"

##Select significant values 
correlation_data <- data.frame(mRNA = row.names(correlation_data), correlation_data)
correlation_data <- correlation_data %>% 
  arrange(desc(abs(Cor)))
#threshold value
correlation_data <- subset(correlation_data, correlation_data$Cor <= -0.45| correlation_data$Cor >= 
                    0.45)
mRNA_col <- c(correlation_data$mRNA)

#sub_data has the gene attributes that satisfy minimum correlation
#constraint with BCR_freetime
sub_data <- data[, names(data) %in% mRNA_col]
#new_data has all the clinical variables with satisfied mrna variables
new_data <- cbind(data[c(1:19, 26467:26470)], sub_data)
f_data <- cbind(data[c(1:19, 26467:26470)], sub_data)

new_data$Race=as.numeric(new_data$Race)
new_data$RP_Type=as.numeric(new_data$RP_Type)


temp <- new_data
temp$Copy.number.Cluster <- NULL
new_data <- temp
temp$ERG.fusion.aCGH <- NULL

new_data$ERG.fusion.aCGH <- NULL

new_data$Type <- as.numeric(new_data$Type)
new_data$BxGGS <- as.numeric(new_data$BxGGS)
new_data$BxGG2 <- as.numeric(new_data$BxGG2)
new_data$LNI <- as.numeric(new_data$LNI)
new_data$PathGG2 <- as.numeric(new_data$PathGG2)
new_data$PathGGS <- as.numeric(new_data$PathGGS)
new_data$Copy.number.Cluster <- as.numeric(new_data$Copy.number.Cluster)
new_data$SMS.x <- as.numeric(new_data$SMS.x)
new_data$newBGG1 <- as.numeric(new_data$newBGG1)
new_data$newPGG1 <- as.numeric(new_data$newPGG1)
new_data$newClinStage <- as.numeric(new_data$newClinStage)
new_data$newPathStage <- as.numeric(new_data$newPathStage)
new_data$ERG.fusion.gex <- as.numeric(new_data$ERG.fusion.gex)
new_data$ECE <- as.numeric(new_data$ECE)
new_data$SVI <- as.numeric(new_data$SVI)
new_data$Copy.number.Cluster= NULL
#---------------------------------------------------------------------------------------------------------

Y=Surv(f_data$BCR_FreeTime, f_data$BCR_Event == 1)

#KM plot

#overall
kmfit1 = survfit(Y~1,data = f_data)


ggsurvplot(kmfit1,xlab = "BCR time in months",title="Overall BCR Curve",ylab = "BCR 
probabilities",conf.int = TRUE)

#Type
kmfit2 = survfit(Y~f_data$Type,data = f_data)
ggsurvplot(kmfit2,xlab = "BCR time in months",title="BCR Curve by Type",ylab = "BCR 
probabilities",pval=TRUE)

#Race

kmfit3 = survfit(Y ~ f_data$Race,data = f_data)
ggsurvplot(kmfit3,xlab = "BCR time in months",title="BCR Curve by Race",ylab = "BCR 
probabilities",pval=TRUE)

#BGG1

kmfit4 = survfit(Y ~ f_data$newBGG1,data = f_data)
ggsurvplot(kmfit4,xlab = "BCR time in months",title="BCR Curve by Primary Gleason Score",ylab = "BCR 
probabilities",pval=TRUE)

#BGG2
kmfit5 = survfit(Y ~ f_data$BxGG2,data = f_data)
ggsurvplot(kmfit5,xlab = "BCR time in months",title="BCR Curve by Secondary Gleason Score",ylab = "BCR 
probabilities",pval=TRUE)

# BGGS
kmfit6 = survfit(Y ~ f_data$BxGGS,data = f_data)
ggsurvplot(kmfit6,xlab = "BCR time in months",title="BCR Curve by Combined Gleason Score",ylab = "BCR 
probabilities",pval=TRUE)

# Clinical Stage
kmfit7 = survfit(Y ~ f_data$newClinStage,data = f_data)
ggsurvplot(kmfit7,xlab = "BCR time in months",title="BCR Curve by Clinical Stage",ylab = "BCR 
probabilities",pval=TRUE)


# SMS
kmfit8 = survfit(Y ~ f_data$SMS.x,data = f_data)
ggsurvplot(kmfit8,xlab = "BCR time in months",title="BCR Curve by SMS",ylab = "BCR 
probabilities",pval=TRUE)

# ECE
kmfit9 = survfit(Y ~ f_data$ECE,data = f_data)
ggsurvplot(kmfit9,xlab = "BCR time in months",title="BCR Curve by ECE",ylab = "BCR 
probabilities",pval=TRUE)

# SVI
kmfit10 = survfit(Y ~ f_data$SVI,data = f_data)
ggsurvplot(kmfit10,xlab = "BCR time in months",title="BCR Curve by SVI",ylab = "BCR 
probabilities",pval=TRUE)

## LNI
kmfit11 = survfit(Y ~ f_data$LNI,data = f_data)
ggsurvplot(kmfit11,xlab = "BCR time in months",title="BCR Curve by LNI",ylab = "BCR 
probabilities",pval=TRUE)

## Path Stage
kmfit12 = survfit(Y ~ f_data$newPathStage,data = f_data)
ggsurvplot(kmfit12,xlab = "BCR time in months",title="BCR Curve by Pathalogical Stage",ylab = "BCR 
probabilities",pval=TRUE)



##PGG1
kmfit13 = survfit(Y ~ f_data$newPGG1,data = f_data)
ggsurvplot(kmfit13,xlab = "BCR time in months",title="BCR Curve by Primary Pathalogical Gleason Score",ylab = "BCR 
probabilities",pval=TRUE)

## PGG2
kmfit14 = survfit(Y ~ f_data$PathGG2,data = f_data)
ggsurvplot(kmfit14,xlab = "BCR time in months",title="BCR Curve by Secondary Pathalogical Gleason Score",ylab = "BCR 
probabilities",pval=TRUE)

## PGGS
kmfit15 = survfit(Y ~ f_data$PathGGS,data = f_data)
ggsurvplot(kmfit15,xlab = "BCR time in months",title="BCR Curve by Combined Pathalogical Gleason Score",ylab = "BCR 
probabilities",pval=TRUE)


kmfit19 = survfit(Y~f_data$ERG.fusion.gex,data = f_data)
ggsurvplot(kmfit19,xlab = "BCR time in months",title="BCR Curve by ERG fusion gex",ylab = "BCR 
probabilities",pval=TRUE)


#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#........................................GNERIC MODEL.................................................#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#

#stepwise model data
step_full <- new_data
step_pre <- cbind(new_data[c(1:6, 15:18, 20)], 
                  sub_data)

#lasso ridge model
penalized_full <- new_data


penalized_pre <- cbind(penalized_full[c(1:6, 15:18, 20)], 
                       sub_data)

#rfsrc
rf_full <- temp
rf_pre <- cbind(temp[c(1:6, 15:18, 20)], 
                sub_data)

#gbm model

gbm_full <- penalized_full
gbm_pre <- penalized_pre


#folds
K=5
gbm_counter = 1
counter = 1

#-------------------------------------------------------------------------------------------
#stepmodel variables

#feature matrix post
subset_cols_step_post = ncol(step_full)-2
subset_rows_setp_post = nrow(step_full)
feat_matrx_step_post = matrix(0,nrow=50,ncol=subset_cols_step_post)
subset_nms_step_post = colnames(subset(step_full, select = -c(BCR_FreeTime, BCR_Event)))
colnames(feat_matrx_step_post) = subset_nms_step_post

#fetaure matrix pre
subset_cols_step_pre = ncol(step_pre)-2
subset_rows_step_pre = nrow(step_pre)
feat_matrx_step_pre = matrix(0,nrow=50,ncol=subset_cols_step_pre)
subset_nms_step_pre = colnames(subset(step_pre, select = -c(BCR_FreeTime, BCR_Event)))
colnames(feat_matrx_step_pre) = subset_nms_step_pre

#pred matrix post
pred_matrx_step_post = matrix(0,nrow=50,ncol=134)


#pred matrix pre
pred_matrx_step_pre = matrix(0,nrow=50,ncol=134)

#list post
full_step_model = matrix(nrow = 10,ncol=5)
full_cv_conc = matrix(nrow = 10,ncol=5)
full_validation_conc = matrix(nrow = 10,ncol=5)
full_validation_se = matrix(nrow = 10,ncol=5)

#list pre
pre_op_step_model <- matrix(nrow = 10,ncol=5)
pre_op_cv_conc <- matrix(nrow = 10,ncol=5)
pre_op_validation_conc <- matrix(nrow = 10,ncol=5)
pre_op_validation_se <- matrix(nrow = 10,ncol=5)

#-------------------------------------------------------------------------------------------
#lasso variables

#pred matrix post
pred_matrx_lasso_test_post = matrix(0,nrow=50,ncol=134)


#pred matrix pre
pred_matrx_lasso_test_pre = matrix(0,nrow=50,ncol=134)

#lasso cindex
lasso.cin.test.post = lasso.cin.test.pre = lasso.cin.train.post = lasso.cin.train.pre =numeric(50)

#lasso coef
subset_nms_lasso_post = colnames(subset(penalized_full, select = -c(BCR_FreeTime, BCR_Event)))
subset_nms_lasso_pre = colnames(subset(penalized_pre, select = -c(BCR_FreeTime, BCR_Event)))


lasso.coef.post = matrix(0,nrow=50,ncol=length(subset_nms_lasso_post))
lasso.coef.pre = matrix(0,nrow=50,ncol=length(subset_nms_lasso_pre))

colnames(lasso.coef.post) = subset_nms_lasso_post
colnames(lasso.coef.pre) = subset_nms_lasso_pre

#--------------------------------------------------------------------------------------------
#ridge variables

#pred matrix post
pred_matrx_ridge_test_post = matrix(0,nrow=50,ncol=134)


#pred matrix pre
pred_matrx_ridge_test_pre = matrix(0,nrow=50,ncol=134)

#ridge cindex
ridge.cin.test.post = ridge.cin.test.pre = ridge.cin.train.post = ridge.cin.train.pre = numeric(50)

#ridge coef
subset_nms_ridge_post = colnames(subset(penalized_full, select = -c(BCR_FreeTime, BCR_Event)))
subset_nms_ridge_pre = colnames(subset(penalized_pre, select = -c(BCR_FreeTime, BCR_Event)))


ridge.coef.post = matrix(0,nrow=50,ncol=length(subset_nms_ridge_post))
ridge.coef.pre = matrix(0,nrow=50,ncol=length(subset_nms_ridge_pre))

colnames(ridge.coef.post) = subset_nms_ridge_post
colnames(ridge.coef.pre) = subset_nms_ridge_pre




#------------------------------------------------------------------------------------------------
#rfsrc variables

#features

rfsrc_features_post <- matrix(0,nrow=50,ncol=subset_cols_step_post)
rfsrc_features_pre <- matrix(0,nrow=50,ncol=subset_cols_step_pre)

rsf_model_full = rsf_pred_full = rsf_model_pre = rsf_pred_pre =  vector("list",50)
rsf_conc_train_full = rsf_conc_test_full =rsf_conc_train_pre = rsf_conc_test_pre = vector("list",50)

subset_col_names_rf_full = colnames(subset(rf_full, select = -c(BCR_FreeTime,BCR_Event)))
subset_col_names_rf_pre = colnames(subset(rf_pre, select = -c(BCR_FreeTime,BCR_Event)))


colnames(rfsrc_features_post) = subset_nms_step_post
colnames(rfsrc_features_pre) = subset_nms_step_pre
#-----------------------------------------------------------------------------------------------
#gbm

#pred matrix post
pred_matrx_gbm_test_post = matrix(0,nrow=50,ncol=134)


#pred matrix pre
pred_matrx_gbm_test_pre = matrix(0,nrow=50,ncol=134)

gbm.coef.post = matrix(0, nrow=50, ncol= subset_cols_step_post)

gbm.coef.pre = matrix(0, nrow=50, ncol= subset_cols_step_pre)

colnames(gbm.coef.post) = subset_nms_step_post

colnames(gbm.coef.pre) = subset_nms_step_pre


gbm_summary_post = gbm.test.conc.post = gbm.train.conc.post =vector("list",50) 
gbm_summary_pre = gbm.test.conc.pre = gbm.train.conc.pre =vector("list",50) 

#---------------------------------------------------------------------------------------------------
set.seed(6090)

step_full[1:20]

rfe_train_conc = vector("integer",100)
rfe_test_conc = vector("integer",100)
variables_picked = list()
for(i in 1:10){
  
  
  #shuffle 
  data_shuffle <- new_data[sample(1:nrow(new_data)),]
  cvIndex <- createFolds(factor(data_shuffle$BCR_Event), K, 
                         returnTrain = T)
  
  
  for(j in 1:length(cvIndex)){
    #train.index = which(cvIndex!=j)
    #test.index = which(cvIndex==j)
    
    #----------------------------------------------------------------------------
    
    
    
    #stepwise model
    train_step_full <- step_full[cvIndex[[j]],]
    train_step_pre <- step_pre[cvIndex[[j]],]
    
    test_step_full <- step_full[-cvIndex[[j]],]
    test_step_pre <- step_pre[-cvIndex[[j]],]
    
    #start model
    start_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ 1, data = 
                         train_step_full)
    #stepwise for post-rp(pre-rp+post-rp)
    full_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ ., data = 
                        train_step_full)
    #stepwise for pre-rp
    pre_op_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ ., data = 
                          train_step_pre)
    
    
    full_fit_step = step(start_cox, direction = "both", scope = 
                           full_cox$formula,
                         steps = 10)
    full_form<-full_fit_step$formula[3]
    full_form1<- formula(paste("Surv(BCR_FreeTime,BCR_Event)~",full_form))
    #full_pro_haz[i,j] <- cox.zph(full_fit_step)
    
    
    #stepwise seletion for pre-rp
    pre_op_fit_step = step(start_cox, direction = "both", scope = 
                             pre_op_cox$formula,
                           steps = 10)
    pre_op_form<-pre_op_fit_step$formula[3]
    pre_op_form1<- formula(paste("Surv(BCR_FreeTime,BCR_Event)~",pre_op_form))
    
    full_coefficients <- names(full_fit_step$coefficients)
    pre_coefficients <- names(pre_op_fit_step$coefficients)
    
    for(m in 1:length(full_coefficients)){
      feat_matrx_step_post[counter,full_coefficients[m]] = 1
    }
    
    
    for(m in 1:length(pre_coefficients)){
      feat_matrx_step_pre[counter,pre_coefficients[m]] = 1
    }
    
    #cox model post
    cox_model_post=coxph(full_form1,data=train_step_full)
    cox_conc_post=predict(cox_model_post,newdata = test_step_full,type="lp")
    
    #cox model pre
    cox_model_pre=coxph(pre_op_form1,data=train_step_pre)
    cox_conc_pre=predict(cox_model_pre,newdata = test_step_pre,type="lp")
    
    
    #for km curve
    for(r in 1:length(cox_conc_post)){
      some_var = cox_conc_post[r]
      indx = which(rownames(step_full)==names(some_var))
      pred_matrx_step_post[counter,indx] = as.numeric(some_var)[1]
    }
    
    for(r in 1:length(cox_conc_pre)){
      some_var = cox_conc_pre[r]
      indx = which(rownames(step_pre)==names(some_var))
      pred_matrx_step_pre[counter,indx] = as.numeric(some_var)[1]
    }
    
    
    
    full_validation_post <- concordance(cox_model_post, newdata = test_step_full)
    full_train_validation<-concordance(cox_model_post, newdata = train_step_full)
    full_cv_conc[i,j] <- full_train_validation$concordance
    full_validation_conc[i,j] <- full_validation_post$concordance
    full_step_model[i,j] <- as.character(full_form)
    full_validation_se[i,j] <- sqrt(full_validation_post$var)
    
    
    pre_op_validation <- concordance(cox_model_pre, newdata = test_step_pre)
    pre_op_train_validation <- concordance(cox_model_pre, newdata = train_step_pre)
    pre_op_cv_conc[i,j]<-pre_op_train_validation$concordance
    pre_op_validation_conc[i,j]<-pre_op_validation$concordance
    pre_op_step_model[i,j]<-as.character(pre_op_form)
    pre_op_validation_se[i,j]<-sqrt(pre_op_validation$var)
    
    #--------------------------------------------------------------------------------------------------
    #lasso & ridge
    # train_penal_full <- penalized_full[cvIndex[[j]],]
    # train_penal_pre <- penalized_pre[cvIndex[[j]],]
    # 
    # test_penal_full <- penalized_full[-cvIndex[[j]],]
    # test_penal_pre <- penalized_pre[-cvIndex[[j]],]
    
    xtrain_post=subset(train_step_full,select=-c(15,16))
    xtrain_pre=subset(train_step_pre,select=-c(7,8))
    ytrain=subset(train_step_full,select=c(15,16))
    
    #Creating the Survival Object -TRAIN DATA
    new_surv_obj=Surv(ytrain$BCR_FreeTime,ytrain$BCR_Event)
    
    #Creating Model matrix for glmnet input -TRAIN DATA
    x_model_post= model.matrix(new_surv_obj~.+0, data=xtrain_post)
    x_model_pre = model.matrix(new_surv_obj~.+0, data=xtrain_pre)
    
    #test data
    x.val = subset(test_step_full, select = -c(BCR_FreeTime, BCR_Event))
    x.val.pre = subset(test_step_pre,select = -c(BCR_FreeTime, BCR_Event))
    y.valid = subset(test_step_full, select = c(BCR_FreeTime, BCR_Event))
    
    #Creating Survival Object -TEST DATA
    surv.val = Surv(y.valid$BCR_FreeTime,y.valid$BCR_Event)
    
    surv.train=Surv(ytrain$BCR_FreeTime,ytrain$BCR_Event)
    
    #Creating Model Matrix for glmnet-TEST data
    x.val.m = model.matrix(surv.val ~.+0,data=x.val)
    x.val.m2 <- model.matrix(surv.val~.+0,data=x.val.pre)
    
    #lasso
    cr.cv.lasso.post = cv.glmnet(x_model_post, new_surv_obj, family="cox", alpha=1)
    cr.lasso.post = glmnet(x_model_post, new_surv_obj, family="cox", lambda=cr.cv.lasso.post$lambda.min, alpha=1)
    pred.lasso.post = predict(cr.lasso.post,newx = x.val.m,type="link")[,1]
    pred.lasso.post.train = predict(cr.lasso.post, newx = x_model_post,type="link")[,1]
    
    for(r in 1:length(pred.lasso.post)){
      some_var = pred.lasso.post[r]
      indx = which(rownames(penalized_full)==names(some_var))
      pred_matrx_lasso_test_post[counter,indx] = as.numeric(some_var)[1]
    }
    
    cr.cv.lasso.pre = cv.glmnet(x_model_pre, new_surv_obj, family="cox", alpha=1)
    cr.lasso.pre = glmnet(x_model_pre, new_surv_obj, family="cox", lambda=cr.cv.lasso.pre$lambda.min, alpha=1)
    pred.lasso.pre = predict(cr.lasso.pre,newx = x.val.m2,type="link")[,1]
    pred.lasso.pre.train = predict(cr.lasso.pre,newx = x_model_pre,type="link")[,1]
    for(r in 1:length(pred.lasso.pre)){
      some_var = pred.lasso.pre[r]
      indx = which(rownames(penalized_pre)==names(some_var))
      pred_matrx_lasso_test_pre[counter,indx] = as.numeric(some_var)[1]
    }
    
    lasso.cin.test.post[counter] = Cindex(pred.lasso.post, y=surv.val)
    lasso.cin.test.pre[counter] =  Cindex(pred.lasso.pre, y=surv.val)
    
    lasso.cin.train.post[counter] = Cindex(pred.lasso.post.train, y=surv.train)
    lasso.cin.train.pre[counter] =  Cindex(pred.lasso.pre.train, y=surv.train)
    
    lasso.coef.post[counter,] = as.numeric(coef(cr.lasso.post)!=0)
    lasso.coef.pre[counter,] = as.numeric(coef(cr.lasso.pre)!=0)
    
    #--------------------------------------------------------------------------------------------
    #ridge
    cr.cv.ridge.post = cv.glmnet(x_model_post, new_surv_obj, family="cox", alpha=0)
    cr.ridge.post = glmnet(x_model_post, new_surv_obj, family="cox", lambda=cr.cv.ridge.post$lambda.min, alpha=0)
    pred.ridge.post = predict(cr.ridge.post,newx = x.val.m,type="link")[,1]
    pred.ridge.post.train = predict(cr.ridge.post, newx = x_model_post,type="link")[,1]
    
    for(r in 1:length(pred.ridge.post)){
      some_var = pred.ridge.post[r]
      indx = which(rownames(penalized_full)==names(some_var))
      pred_matrx_ridge_test_post[counter,indx] = as.numeric(some_var)[1]
    }
    
    cr.cv.ridge.pre = cv.glmnet(x_model_pre, new_surv_obj, family="cox", alpha=0)
    cr.ridge.pre = glmnet(x_model_pre, new_surv_obj, family="cox", lambda=cr.cv.ridge.pre$lambda.min, alpha=0)
    pred.ridge.pre = predict(cr.ridge.pre,newx = x.val.m2,type="link")[,1]
    pred.ridge.pre.train = predict(cr.ridge.pre,newx = x_model_pre,type="link")[,1]
    for(r in 1:length(pred.ridge.pre)){
      some_var = pred.ridge.pre[r]
      indx = which(rownames(penalized_pre)==names(some_var))
      pred_matrx_ridge_test_pre[counter,indx] = as.numeric(some_var)[1]
    }
    
    ridge.cin.test.post[counter] = Cindex(pred.ridge.post, y=surv.val)
    ridge.cin.test.pre[counter] =  Cindex(pred.ridge.pre, y=surv.val)
    
    ridge.cin.train.post[counter] = Cindex(pred.ridge.post.train, y=surv.train)
    ridge.cin.train.pre[counter] =  Cindex(pred.ridge.pre.train, y=surv.train)
    
    #calulation
    post_coef_sum = sum(abs(coef(cr.ridge.post)))
    post_list = numeric(length  = length(coef(cr.ridge.post)))
    for(cf in 1:length(post_list)){
      post_list[cf]=as.numeric(abs(coef(cr.ridge.post)[cf])/post_coef_sum)
    }
    dec_list_post<-sort.list(post_list,decreasing = TRUE)
    threshold_post <- dec_list_post[10]
    ridge.coef.post[counter,]<- as.numeric(abs(coef(cr.ridge.post))>=abs(coef(cr.ridge.post)[threshold_post]))
    
    
    pre_coef_sum = sum(abs(coef(cr.ridge.pre)))
    pre_list = numeric(length  = length(coef(cr.ridge.pre)))
    for(cf in 1:length(pre_list)){
      pre_list[cf]=as.numeric(abs(coef(cr.ridge.pre)[cf])/pre_coef_sum)
    }
    dec_list_pre<-sort.list(pre_list,decreasing = TRUE)
    threshold_pre <- dec_list_pre[10]
    ridge.coef.pre[counter,]<- as.numeric(abs(coef(cr.ridge.pre))>=abs(coef(cr.ridge.pre)[threshold_pre]))
    #--------------------------------------------------------------------------------------------------------------
    
    #gbm
    gbm_model_post = gbm(surv.train~., data=xtrain_post, distribution="coxph")
    gbm_model_pre = gbm(surv.train~., data=xtrain_pre, distribution="coxph")
    
    
    gbm_summary_post[[gbm_counter]] = summary(gbm_model_post)
    gbm_summary_pre[[gbm_counter]] = summary(gbm_model_pre)
    
    gbm_pred_post= predict(gbm_model_post, newdata=x.val, type="link")
    gbm_pred_pre= predict(gbm_model_pre, newdata=x.val.pre, type="link")
    
    gbm_pred_post_train= predict(gbm_model_post, newdata=xtrain_post, type="link")
    gbm_pred_pre_train= predict(gbm_model_pre, newdata=xtrain_pre, type="link")
    
    
    post_index =  rownames(x.val)
    pre_index = rownames(x.val.pre)
    
    
    
    
    
    
    
    
    for(r in 1:length(post_index)){
      some_var = gbm_pred_post[r]
      name = post_index[r]
      indx = which(rownames(gbm_full)==name)
      pred_matrx_gbm_test_post[gbm_counter,indx] = as.numeric(some_var)[1]
    }
    
    for(r in 1:length(gbm_pred_pre)){
      some_var = gbm_pred_pre[r]
      name = pre_index[r]
      indx = which(rownames(gbm_pre)==name)
      pred_matrx_gbm_test_pre[counter,indx] = as.numeric(some_var)[1]
    }
    
    gbm.test.conc.post[gbm_counter] = Cindex(gbm_pred_post, y=surv.val)
    gbm.test.conc.pre[gbm_counter] = Cindex(gbm_pred_pre, y=surv.val)
    
    gbm.train.conc.post[gbm_counter] = Cindex(gbm_pred_post_train, y=surv.train)
    gbm.train.conc.pre[gbm_counter] = Cindex(gbm_pred_pre_train, y=surv.train)
    #---------------------------------------------------------------------------------------------------
    #rfsrc
    
    #post rp
    rsf_full = rfsrc(Surv(BCR_FreeTime,BCR_Event) ~ .,data = train_step_full, importance = TRUE)
    pred_full = predict.rfsrc(rsf_full,newdata = x.val,type="lp")
    
    
    #prerp
    rsf_pre = rfsrc(Surv(BCR_FreeTime,BCR_Event) ~ .,data = train_step_pre, importance = TRUE)
    pred_pre = predict.rfsrc(rsf_pre,newdata = x.val.pre,type="lp")
    
    
    rsf_model_full[[counter]] = rsf_full
    rsf_pred_full[[counter]] = pred_full
    rsf_conc_test_full[counter]= Cindex(pred_full$predicted,y=surv.val)
    rsf_conc_train_full[counter] =Cindex(rsf_full$predicted.oob,y=surv.train)
    
    rsf_model_pre[[counter]] = rsf_pre
    rsf_pred_pre[[counter]] = pred_pre
    rsf_conc_test_pre[counter]= Cindex(pred_pre$predicted,y=surv.val)
    rsf_conc_train_pre[counter] =Cindex(rsf_pre$predicted.oob,y=surv.train)
    
    rsf_feats_post = head(subset_col_names_rf_full[order(rsf_full$importance, decreasing=TRUE)],10)
    for(m in 1:length(rsf_feats_post)){
      rfsrc_features_post[counter,rsf_feats_post[m]]=1
    }
    
    rsf_feats_pre = head(subset_col_names_rf_pre[order(rsf_pre$importance, decreasing=TRUE)],10)
    for(m in 1:length(rsf_feats_pre)){
      rfsrc_features_pre[counter,rsf_feats_pre[m]]=1
    }
    
  
    counter = counter + 1
    gbm_counter = gbm_counter + 1
  }
  
}

split.kms <- function(zos,os,oevn,NQ=100,zqmin=.05,zqmax=.95){
  qmin = quantile(zos,zqmin,na.rm=TRUE)
  qmax = quantile(zos,zqmax,na.rm=TRUE)
  p0q = numeric(NQ)
  med0s = seq(qmin,qmax,length=NQ)
  izs = matrix(NA,ncol=NQ,nrow=length(zos))
  for(iq in 1:NQ){
    IZ = as.integer(zos<med0s[iq])
    p0q[iq] = 1-pchisq(survdiff(Surv(os,oevn)~factor(IZ))$chisq,1)
    izs[,iq] = IZ
  }
  best.med0 = med0s[which.min(p0q)]
  IZ = as.integer(zos<best.med0)
  return(list(iz=IZ,izs=izs,best.t=best.med0,
              ps=p0q,best.p=p0q[which.min(p0q)]))
}

#result

#stepwise post -rp
median(train_set)
std.error(train_set)
median(unlist(full_validation_conc))

train_set = numeric(50)
for(i in 1:50){
  train_set[i]=full_cv_conc[i]
}


test_set = numeric(50)
for(i in 1:50){
  test_set[i]=full_validation_conc[i]
}
std.error(test_set)

#lasso
median(lasso.cin.train.post)
std.error(lasso.cin.train.post)

median(lasso.cin.test.post)
std.error(lasso.cin.test.post)



#ridge
median(ridge.cin.train.post)
std.error(ridge.cin.train.post)

median(ridge.cin.test.post)
std.error(ridge.cin.test.post)

#gbm
median(unlist(gbm.train.conc.post))
std.error(unlist(gbm.train.conc.post))

median(unlist(gbm.test.conc.post))
std.error(unlist(gbm.test.conc.post))

#rfsrc
median(unlist(rsf_conc_train_full))
std.error(unlist(rsf_conc_train_full))

median(unlist(rsf_conc_test_full))
std.error(unlist(rsf_conc_test_full))


#post rp features

which(colSums(feat_matrx_step_post)>=15)


which(colSums(lasso.coef.post)>=15)

which(colSums(ridge.coef.post)>=15)


which(colSums(rfsrc_features_post)>=15)

which(colSums(gbm))

for(m in 1:50 ){
  for(n in 1:36){
    gbm_index_var=which(subset_nms_step_post==gbm_summary_post[[m]][n,]$var)
    gbm.coef.post[m,gbm_index_var]=gbm_summary_post[[m]][n,]$rel.inf
  }
}

sort(colSums(gbm.coef.post))



boxplot(train_set,lasso.cin.train.post,ridge.cin.train.post,unlist(gbm.train.conc.post),unlist(rsf_conc_train_full),names=c("Stepwise","Lasso","Ridge", "GBM","RSF"),
        col=c(0,2,7,4,5),main="Train Concordances")
boxplot(test_set,lasso.cin.test.post,ridge.cin.test.post,unlist(gbm.test.conc.post),unlist(rsf_conc_test_full),names=c("Stepwise","Lasso","Ridge", "GBM","RSF"),
        col=c(0,2,7,4,5),main="Test Concordances")


boxplot(train_set,test_set,lasso.cin.train.post,lasso.cin.test.post,ridge.cin.train.post,ridge.cin.test.post,unlist(gbm.train.conc.post),unlist(gbm.test.conc.post),unlist(rsf_conc_train_full),unlist(rsf_conc_test_full),names=c("Step Train","Step Test","Lasso Train","Lasso Test","Ridge Train","Ridge Test","GBM Train","GBM Test","RSF Train","RSF Test"),col=c(0,0,2,2,7,7,4,4,5,5))
?boxplot

#kmcurve

stp_mean = as.array(apply(pred_matrx_step_post,2,mean))
stp_kmo = split.kms(stp_mean, os=step_full$BCR_FreeTime, oevn=step_full$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(step_full$BCR_FreeTime, step_full$BCR_Event)~stp_temp_var, data=step_full)
ggsurvplot(stp_km.split,title="BCR Free Survival for Stepwise Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"),)
?ggsurvplot


lasso_mean = as.array(apply(pred_matrx_lasso_test_post,2,mean))
stp_kmo = split.kms(lasso_mean, os=step_full$BCR_FreeTime, oevn=step_full$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(step_full$BCR_FreeTime, step_full$BCR_Event)~stp_temp_var, data=step_full)
ggsurvplot(stp_km.split,title="BCR Free Survival for Lasso Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"),surv.median.line = 'hv')



hridge_mean = as.array(apply(pred_matrx_ridge_test_post,2,mean))
stp_kmo = split.kms(ridge_mean, os=step_full$BCR_FreeTime, oevn=step_full$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(step_full$BCR_FreeTime, step_full$BCR_Event)~stp_temp_var, data=step_full)
ggsurvplot(stp_km.split,title="BCR Free Survival for Ridge Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))


gbm_mean = as.array(apply(pred_matrx_gbm_test_post,2,mean))
stp_kmo = split.kms(gbm_mean, os=step_full$BCR_FreeTime, oevn=step_full$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(step_full$BCR_FreeTime, step_full$BCR_Event)~stp_temp_var, data=step_full)
ggsurvplot(stp_km.split,title="BCR Free Survival for GBM Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))


#kmcurve pre

stp_mean = as.array(apply(pred_matrx_step_pre,2,mean))
stp_kmo = split.kms(stp_mean, os=step_pre$BCR_FreeTime, oevn=step_pre$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(step_pre$BCR_FreeTime, step_pre$BCR_Event)~stp_temp_var, data=step_pre)
ggsurvplot(stp_km.split,title="BCR Free Survival for Stepwise Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))



lasso_mean = as.array(apply(pred_matrx_lasso_test_pre,2,mean))
stp_kmo = split.kms(lasso_mean, os=step_pre$BCR_FreeTime, oevn=step_pre$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(step_pre$BCR_FreeTime, step_pre$BCR_Event)~stp_temp_var, data=step_pre)
ggsurvplot(stp_km.split,title="BCR Free Survival for Lasso Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))



ridge_mean = as.array(apply(pred_matrx_ridge_test_pre,2,mean))
stp_kmo = split.kms(ridge_mean, os=step_pre$BCR_FreeTime, oevn=step_pre$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(step_pre$BCR_FreeTime, step_pre$BCR_Event)~stp_temp_var, data=step_pre)
ggsurvplot(stp_km.split,title="BCR Free Survival for Ridge Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))


gbm_mean = as.array(apply(pred_matrx_gbm_test_pre,2,mean))
stp_kmo = split.kms(gbm_mean, os=step_pre$BCR_FreeTime, oevn=step_pre$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(step_pre$BCR_FreeTime, step_pre$BCR_Event)~stp_temp_var, data=step_pre)
ggsurvplot(stp_km.split,title="BCR Free Survival for GBM Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))


#pre rp


#stepwise pre -rp
median(train_set)
std.error(train_set)
median(test_set)

train_set = numeric(50)
for(i in 1:50){
  train_set[i]=pre_op_cv_conc[i]
}


test_set = numeric(50)
for(i in 1:50){
  test_set[i]=pre_op_validation_conc[i]
}
std.error(test_set)

#lasso
median(lasso.cin.train.pre)
std.error(lasso.cin.train.pre)

median(lasso.cin.test.pre)
std.error(lasso.cin.test.pre)



#ridge
median(ridge.cin.train.pre)
std.error(ridge.cin.train.pre)

median(ridge.cin.test.pre)
std.error(ridge.cin.test.pre)

#gbm
median(unlist(gbm.train.conc.pre))
std.error(unlist(gbm.train.conc.pre))

median(unlist(gbm.test.conc.pre))
std.error(unlist(gbm.test.conc.pre))

#rfsrc
median(unlist(rsf_conc_train_pre))
std.error(unlist(rsf_conc_train_pre))

median(unlist(rsf_conc_test_pre))
std.error(unlist(rsf_conc_test_pre))


boxplot(train_set,lasso.cin.train.pre,ridge.cin.train.pre,unlist(gbm.train.conc.pre),unlist(rsf_conc_train_pre),names=c("Stepwise","Lasso","Ridge", "GBM","RSF"),
        col=c(0,2,7,4,5),main="Train Concordances")


boxplot(test_set,lasso.cin.test.pre,ridge.cin.test.pre,unlist(gbm.test.conc.pre),unlist(rsf_conc_test_pre),names=c("Stepwise","Lasso","Ridge", "GBM","RSF"),
        col=c(0,2,7,4,5),main="Test Concordances")


boxplot(train_set,test_set,lasso.cin.train.pre,lasso.cin.test.pre,ridge.cin.train.pre,ridge.cin.test.pre,unlist(gbm.train.conc.pre),unlist(gbm.test.conc.pre),unlist(rsf_conc_train_full),unlist(rsf_conc_test_full),names=c("Step Train","Step Test","Lasso Train","Lasso Test","Ridge Train","Ridge Test","GBM Train","GBM Test","RSF Train","RSF Test"),col=c(0,0,2,2,7,7,4,4,5,5))


#pre rp features

which(colSums(feat_matrx_step_pre)>=15)


which(colSums(lasso.coef.pre)>=15)

which(colSums(ridge.coef.pre)>=15)


which(colSums(rfsrc_features_pre)>=15)

which(colSums(gbm))

for(m in 1:50 ){
  for(n in 1:26){
    gbm_index_var=which(subset_nms_step_pre==gbm_summary_pre[[m]][n,]$var)
    gbm.coef.pre[m,gbm_index_var]=gbm_summary_pre[[m]][n,]$rel.inf
  }
}

sort(colSums(gbm.coef.pre))

#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#........................................OPTIMAL MODEL.................................................#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
K=5
gbm_counter = 1
counter = 1

#-------------------------------------------------------------------------------------------
#stepmodel variables




post_data_cox <- subset(post_data_cox,select = c("PathGGS","newPGG1","newPathStage","SVI","LNI","PreTxPSA","EMP1","SRD5A2","KLF4","TIPARP","BCR_FreeTime","BCR_Event" ))
pre_data_cox <- subset(pre_data_cox,select = c("PreDxBxPSA" ,"Type","newBGG1","newClinStage","EMP1","SRD5A2","KLF4","SGK1","FERMT2","BCR_FreeTime","BCR_Event"))

#rfsrc
#posr

post_data_rf <- subset(post_data_rf,select = c("PathGGS","newPGG1","newPathStage","SVI","LNI","PreTxPSA","EMP1","SRD5A2","KLF4","TIPARP","BCR_FreeTime","BCR_Event"))
pre_data_rf <- subset(pre_data_rf,select = c("PreDxBxPSA" ,"Type","newBGG1","newClinStage","EMP1","SRD5A2","KLF4","SGK1","FERMT2","BCR_FreeTime","BCR_Event"))

penalized_full <- post_data_cox
penalized_pre <- pre_data_cox
#feature matrix post
subset_cols_step_post = ncol(post_data_cox)-2
subset_rows_setp_post = nrow(post_data_cox)
feat_matrx_step_post = matrix(0,nrow=50,ncol=subset_cols_step_post)
subset_nms_step_post = colnames(subset(post_data_cox, select = -c(BCR_FreeTime, BCR_Event)))
colnames(feat_matrx_step_post) = subset_nms_step_post

#fetaure matrix pre
subset_cols_pre_data_cox = ncol(pre_data_cox)-2
subset_rows_pre_data_cox = nrow(pre_data_cox)
feat_matrx_pre_data_cox = matrix(0,nrow=50,ncol=subset_cols_pre_data_cox)
subset_nms_pre_data_cox = colnames(subset(pre_data_cox, select = -c(BCR_FreeTime, BCR_Event)))
colnames(feat_matrx_pre_data_cox) = subset_nms_pre_data_cox

#pred matrix post
pred_matrx_step_post = matrix(0,nrow=50,ncol=134)


#pred matrix pre
pred_matrx_pre_data_cox = matrix(0,nrow=50,ncol=134)

#list post
full_step_model = matrix(nrow = 10,ncol=5)
full_cv_conc = matrix(nrow = 10,ncol=5)
full_validation_conc = matrix(nrow = 10,ncol=5)
full_validation_se = matrix(nrow = 10,ncol=5)

#list pre
pre_op_step_model <- matrix(nrow = 10,ncol=5)
pre_op_cv_conc <- matrix(nrow = 10,ncol=5)
pre_op_validation_conc <- matrix(nrow = 10,ncol=5)
pre_op_validation_se <- matrix(nrow = 10,ncol=5)

#-------------------------------------------------------------------------------------------
#lasso variables

#pred matrix post
pred_matrx_lasso_test_post = matrix(0,nrow=50,ncol=134)


#pred matrix pre
pred_matrx_lasso_test_pre = matrix(0,nrow=50,ncol=134)

#lasso cindex
lasso.cin.test.post = lasso.cin.test.pre = lasso.cin.train.post = lasso.cin.train.pre =numeric(50)

#lasso coef
subset_nms_lasso_post = colnames(subset(penalized_full, select = -c(BCR_FreeTime, BCR_Event)))
subset_nms_lasso_pre = colnames(subset(penalized_pre, select = -c(BCR_FreeTime, BCR_Event)))


lasso.coef.post = matrix(0,nrow=50,ncol=length(subset_nms_lasso_post))
lasso.coef.pre = matrix(0,nrow=50,ncol=length(subset_nms_lasso_pre))

colnames(lasso.coef.post) = subset_nms_lasso_post
colnames(lasso.coef.pre) = subset_nms_lasso_pre

#--------------------------------------------------------------------------------------------
#ridge variables

#pred matrix post
pred_matrx_ridge_test_post = matrix(0,nrow=50,ncol=134)


#pred matrix pre
pred_matrx_ridge_test_pre = matrix(0,nrow=50,ncol=134)

#ridge cindex
ridge.cin.test.post = ridge.cin.test.pre = ridge.cin.train.post = ridge.cin.train.pre = numeric(50)

#ridge coef
subset_nms_ridge_post = colnames(subset(penalized_full, select = -c(BCR_FreeTime, BCR_Event)))
subset_nms_ridge_pre = colnames(subset(penalized_pre, select = -c(BCR_FreeTime, BCR_Event)))


ridge.coef.post = matrix(0,nrow=50,ncol=length(subset_nms_ridge_post))
ridge.coef.pre = matrix(0,nrow=50,ncol=length(subset_nms_ridge_pre))

colnames(ridge.coef.post) = subset_nms_ridge_post
colnames(ridge.coef.pre) = subset_nms_ridge_pre




#------------------------------------------------------------------------------------------------
#rfsrc variables

#features

rfsrc_features_post <- matrix(0,nrow=50,ncol=subset_cols_step_post)
rfsrc_features_pre <- matrix(0,nrow=50,ncol=subset_cols_pre_data_cox)

rsf_model_full = rsf_pred_full = rsf_model_pre = rsf_pred_pre =  vector("list",50)
rsf_conc_train_full = rsf_conc_test_full =rsf_conc_train_pre = rsf_conc_test_pre = vector("list",50)

subset_col_names_post_data_rf = colnames(subset(post_data_rf, select = -c(BCR_FreeTime,BCR_Event)))
subset_col_names_pre_data_rf = colnames(subset(pre_data_rf, select = -c(BCR_FreeTime,BCR_Event)))


colnames(rfsrc_features_post) = subset_nms_step_post
colnames(rfsrc_features_pre) = subset_nms_pre_data_cox
#-----------------------------------------------------------------------------------------------
#gbm

#pred matrix post
pred_matrx_gbm_test_post = matrix(0,nrow=50,ncol=134)


#pred matrix pre
pred_matrx_gbm_test_pre = matrix(0,nrow=50,ncol=134)

gbm.coef.post = matrix(0, nrow=50, ncol= subset_cols_step_post)

gbm.coef.pre = matrix(0, nrow=50, ncol= subset_cols_pre_data_cox)

colnames(gbm.coef.post) = subset_nms_step_post

colnames(gbm.coef.pre) = subset_nms_pre_data_cox


gbm_summary_post = gbm.test.conc.post = gbm.train.conc.post =vector("list",50) 
gbm_summary_pre = gbm.test.conc.pre = gbm.train.conc.pre =vector("list",50) 

#---------------------------------------------------------------------------------------------------
set.seed(6090)





for(i in 1:10){
  #shuffle 
  data_shuffle <- new_data[sample(1:nrow(new_data)),]
  cvIndex <- createFolds(factor(data_shuffle$BCR_Event), K, 
                         returnTrain = T)
  
  
  for(j in 1:length(cvIndex)){
    #stepwise model
    train_post_data_cox <- post_data_cox[cvIndex[[j]],]
    train_pre_data_cox <- pre_data_cox[cvIndex[[j]],]
    
    test_post_data_cox <- post_data_cox[-cvIndex[[j]],]
    test_pre_data_cox <- pre_data_cox[-cvIndex[[j]],]
    
    
    #post model
    full_mod = coxph(Surv(BCR_FreeTime,BCR_Event)~.,data = train_post_data_cox)
    full_validation_post <- concordance(full_mod, newdata = test_post_data_cox)
    full_train_validation<-concordance(full_mod, newdata = train_post_data_cox)
    full_cv_conc[i,j] <- full_train_validation$concordance
    full_validation_conc[i,j] <- full_validation_post$concordance
    
    cox_conc_post=predict(full_mod,newdata = test_post_data_cox,type="lp")

    #pre model
    pre_mod = coxph(Surv(BCR_FreeTime,BCR_Event)~.,data = train_pre_data_cox)
    pre_op_validation <- concordance(pre_mod, newdata = test_pre_data_cox)
    pre_op_train_validation <- concordance(pre_mod, newdata = train_pre_data_cox)
    pre_op_cv_conc[i,j]<-pre_op_train_validation$concordance
    pre_op_validation_conc[i,j]<-pre_op_validation$concordance
    
    cox_conc_pre=predict(pre_mod,newdata = test_pre_data_cox,type="lp")
    #for km curve
    for(r in 1:length(cox_conc_post)){
      some_var = cox_conc_post[r]
      indx = which(rownames(step_full)==names(some_var))
      pred_matrx_step_post[counter,indx] = as.numeric(some_var)[1]
    }
    
    for(r in 1:length(cox_conc_pre)){
      some_var = cox_conc_pre[r]
      indx = which(rownames(pre_data_cox)==names(some_var))
      pred_matrx_pre_data_cox[counter,indx] = as.numeric(some_var)[1]
    }
    

    
    xtrain_post=subset(train_post_data_cox,select=-c(11,12))
    xtrain_pre=subset(train_pre_data_cox,select=-c(10,11))
    ytrain=subset(train_post_data_cox,select=c(11,12))
    
    #Creating the Survival Object -TRAIN DATA
    new_surv_obj=Surv(ytrain$BCR_FreeTime,ytrain$BCR_Event)
    
    #Creating Model matrix for glmnet input -TRAIN DATA
    x_model_post= model.matrix(new_surv_obj~.+0, data=xtrain_post)
    x_model_pre = model.matrix(new_surv_obj~.+0, data=xtrain_pre)
    
    #test data
    x.val = subset(test_post_data_cox, select = -c(BCR_FreeTime, BCR_Event))
    x.val.pre = subset(test_pre_data_cox,select = -c(BCR_FreeTime, BCR_Event))
    y.valid = subset(test_post_data_cox, select = c(BCR_FreeTime, BCR_Event))
    
    #Creating Survival Object -TEST DATA
    surv.val = Surv(y.valid$BCR_FreeTime,y.valid$BCR_Event)
    
    surv.train=Surv(ytrain$BCR_FreeTime,ytrain$BCR_Event)
    
    #Creating Model Matrix for glmnet-TEST data
    x.val.m = model.matrix(surv.val ~.+0,data=x.val)
    x.val.m2 <- model.matrix(surv.val~.+0,data=x.val.pre)
    
    #lasso
    cr.cv.lasso.post = cv.glmnet(x_model_post, new_surv_obj, family="cox", alpha=1)
    cr.lasso.post = glmnet(x_model_post, new_surv_obj, family="cox", lambda=cr.cv.lasso.post$lambda.min, alpha=1)
    pred.lasso.post = predict(cr.lasso.post,newx = x.val.m,type="link")[,1]
    pred.lasso.post.train = predict(cr.lasso.post, newx = x_model_post,type="link")[,1]
    
    for(r in 1:length(pred.lasso.post)){
      some_var = pred.lasso.post[r]
      indx = which(rownames(penalized_full)==names(some_var))
      pred_matrx_lasso_test_post[counter,indx] = as.numeric(some_var)[1]
    }
    
    cr.cv.lasso.pre = cv.glmnet(x_model_pre, new_surv_obj, family="cox", alpha=1)
    cr.lasso.pre = glmnet(x_model_pre, new_surv_obj, family="cox", lambda=cr.cv.lasso.pre$lambda.min, alpha=1)
    pred.lasso.pre = predict(cr.lasso.pre,newx = x.val.m2,type="link")[,1]
    pred.lasso.pre.train = predict(cr.lasso.pre,newx = x_model_pre,type="link")[,1]
    for(r in 1:length(pred.lasso.pre)){
      some_var = pred.lasso.pre[r]
      indx = which(rownames(penalized_pre)==names(some_var))
      pred_matrx_lasso_test_pre[counter,indx] = as.numeric(some_var)[1]
    }
    
    lasso.cin.test.post[counter] = Cindex(pred.lasso.post, y=surv.val)
    lasso.cin.test.pre[counter] =  Cindex(pred.lasso.pre, y=surv.val)
    
    lasso.cin.train.post[counter] = Cindex(pred.lasso.post.train, y=surv.train)
    lasso.cin.train.pre[counter] =  Cindex(pred.lasso.pre.train, y=surv.train)

    #--------------------------------------------------------------------------------------------
    #ridge
    cr.cv.ridge.post = cv.glmnet(x_model_post, new_surv_obj, family="cox", alpha=0)
    cr.ridge.post = glmnet(x_model_post, new_surv_obj, family="cox", lambda=cr.cv.ridge.post$lambda.min, alpha=0)
    pred.ridge.post = predict(cr.ridge.post,newx = x.val.m,type="link")[,1]
    pred.ridge.post.train = predict(cr.ridge.post, newx = x_model_post,type="link")[,1]
    
    for(r in 1:length(pred.ridge.post)){
      some_var = pred.ridge.post[r]
      indx = which(rownames(penalized_full)==names(some_var))
      pred_matrx_ridge_test_post[counter,indx] = as.numeric(some_var)[1]
    }
    
    cr.cv.ridge.pre = cv.glmnet(x_model_pre, new_surv_obj, family="cox", alpha=0)
    cr.ridge.pre = glmnet(x_model_pre, new_surv_obj, family="cox", lambda=cr.cv.ridge.pre$lambda.min, alpha=0)
    pred.ridge.pre = predict(cr.ridge.pre,newx = x.val.m2,type="link")[,1]
    pred.ridge.pre.train = predict(cr.ridge.pre,newx = x_model_pre,type="link")[,1]
    for(r in 1:length(pred.ridge.pre)){
      some_var = pred.ridge.pre[r]
      indx = which(rownames(penalized_pre)==names(some_var))
      pred_matrx_ridge_test_pre[counter,indx] = as.numeric(some_var)[1]
    }
    
    ridge.cin.test.post[counter] = Cindex(pred.ridge.post, y=surv.val)
    ridge.cin.test.pre[counter] =  Cindex(pred.ridge.pre, y=surv.val)
    
    ridge.cin.train.post[counter] = Cindex(pred.ridge.post.train, y=surv.train)
    ridge.cin.train.pre[counter] =  Cindex(pred.ridge.pre.train, y=surv.train)
    

    #gbm
    gbm_model_post = gbm(surv.train~., data=xtrain_post, distribution="coxph")
    gbm_model_pre = gbm(surv.train~., data=xtrain_pre, distribution="coxph")
    
    
    gbm_summary_post[[gbm_counter]] = summary(gbm_model_post)
    gbm_summary_pre[[gbm_counter]] = summary(gbm_model_pre)
    
    gbm_pred_post= predict(gbm_model_post, newdata=x.val, type="link")
    gbm_pred_pre= predict(gbm_model_pre, newdata=x.val.pre, type="link")
    
    gbm_pred_post_train= predict(gbm_model_post, newdata=xtrain_post, type="link")
    gbm_pred_pre_train= predict(gbm_model_pre, newdata=xtrain_pre, type="link")
    
    
    post_index =  rownames(x.val)
    pre_index = rownames(x.val.pre)
    
    
    
    for(r in 1:length(post_index)){
      some_var = gbm_pred_post[r]
      name = post_index[r]
      indx = which(rownames(post_data_cox)==name)
      pred_matrx_gbm_test_post[gbm_counter,indx] = as.numeric(some_var)[1]
    }
    
    for(r in 1:length(gbm_pred_pre)){
      some_var = gbm_pred_pre[r]
      name = pre_index[r]
      indx = which(rownames(pre_data_cox)==name)
      pred_matrx_gbm_test_pre[counter,indx] = as.numeric(some_var)[1]
    }
    
    gbm.test.conc.post[gbm_counter] = Cindex(gbm_pred_post, y=surv.val)
    gbm.test.conc.pre[gbm_counter] = Cindex(gbm_pred_pre, y=surv.val)
    
    gbm.train.conc.post[gbm_counter] = Cindex(gbm_pred_post_train, y=surv.train)
    gbm.train.conc.pre[gbm_counter] = Cindex(gbm_pred_pre_train, y=surv.train)
    #---------------------------------------------------------------------------------------------------
    #rfsrc
    
    #post rp
    rsf_full = rfsrc(Surv(BCR_FreeTime,BCR_Event) ~ .,data = train_post_data_cox, importance = TRUE)
    pred_full = predict.rfsrc(rsf_full,newdata = x.val,type="lp")
    
    
    #prerp
    rsf_pre = rfsrc(Surv(BCR_FreeTime,BCR_Event) ~ .,data = train_pre_data_cox, importance = TRUE)
    pred_pre = predict.rfsrc(rsf_pre,newdata = x.val.pre,type="lp")
    
    
    rsf_model_full[[counter]] = rsf_full
    rsf_pred_full[[counter]] = pred_full
    rsf_conc_test_full[counter]= Cindex(pred_full$predicted,y=surv.val)
    rsf_conc_train_full[counter] =Cindex(rsf_full$predicted.oob,y=surv.train)
    
    rsf_model_pre[[counter]] = rsf_pre
    rsf_pred_pre[[counter]] = pred_pre
    rsf_conc_test_pre[counter]= Cindex(pred_pre$predicted,y=surv.val)
    rsf_conc_train_pre[counter] =Cindex(rsf_pre$predicted.oob,y=surv.train)
    

    counter = counter + 1
    gbm_counter = gbm_counter + 1
    
    
    
  }
  
  
}

#result

#stepwise post -rp
median(train_set)
std.error(train_set)
median(unlist(full_validation_conc))

train_set = numeric(50)
for(i in 1:50){
  train_set[i]=full_cv_conc[i]
}


test_set = numeric(50)
for(i in 1:50){
  test_set[i]=full_validation_conc[i]
}
std.error(test_set)

#lasso
median(lasso.cin.train.post)
std.error(lasso.cin.train.post)

median(lasso.cin.test.post)
std.error(lasso.cin.test.post)



#ridge
median(ridge.cin.train.post)
std.error(ridge.cin.train.post)

median(ridge.cin.test.post)
std.error(ridge.cin.test.post)

#gbm
median(unlist(gbm.train.conc.post))
std.error(unlist(gbm.train.conc.post))

median(unlist(gbm.test.conc.post))
std.error(unlist(gbm.test.conc.post))

#rfsrc
median(unlist(rsf_conc_train_full))
std.error(unlist(rsf_conc_train_full))

median(unlist(rsf_conc_test_full))
std.error(unlist(rsf_conc_test_full))

boxplot(train_set,lasso.cin.train.post,ridge.cin.train.post,unlist(gbm.train.conc.post),unlist(rsf_conc_train_full),names=c("Stepwise","Lasso","Ridge", "GBM","RSF"),
        col=c(0,2,7,4,5),main="Train Concordances")
boxplot(test_set,lasso.cin.test.post,ridge.cin.test.post,unlist(gbm.test.conc.post),unlist(rsf_conc_test_full),names=c("Stepwise","Lasso","Ridge", "GBM","RSF"),
        col=c(0,2,7,4,5),main="Test Concordances")
boxplot(train_set,test_set,lasso.cin.train.post,lasso.cin.test.post,ridge.cin.train.post,ridge.cin.test.post,unlist(gbm.train.conc.post),unlist(gbm.test.conc.post),unlist(rsf_conc_train_full),unlist(rsf_conc_test_full),names=c("Step Train","Step Test","Lasso Train","Lasso Test","Ridge Train","Ridge Test","GBM Train","GBM Test","RSF Train","RSF Test"),col=c(0,0,2,2,7,7,4,4,5,5))
?boxplot

#kmcurve

stp_mean = as.array(apply(pred_matrx_step_post,2,mean))
stp_kmo = split.kms(stp_mean, os=post_data_cox$BCR_FreeTime, oevn=post_data_cox$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(post_data_cox$BCR_FreeTime, post_data_cox$BCR_Event)~stp_temp_var, data=post_data_cox)
ggsurvplot(stp_km.split,title="BCR Free Survival for Stepwise Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"),)
?ggsurvplot


lasso_mean = as.array(apply(pred_matrx_lasso_test_post,2,mean))
stp_kmo = split.kms(lasso_mean, os=post_data_cox$BCR_FreeTime, oevn=post_data_cox$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(post_data_cox$BCR_FreeTime, post_data_cox$BCR_Event)~stp_temp_var, data=post_data_cox)
ggsurvplot(stp_km.split,title="BCR Free Survival for Lasso Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))



hridge_mean = as.array(apply(pred_matrx_ridge_test_post,2,mean))
stp_kmo = split.kms(ridge_mean, os=post_data_cox$BCR_FreeTime, oevn=post_data_cox$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(post_data_cox$BCR_FreeTime, post_data_cox$BCR_Event)~stp_temp_var, data=post_data_cox)
ggsurvplot(stp_km.split,title="BCR Free Survival for Ridge Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))


gbm_mean = as.array(apply(pred_matrx_gbm_test_post,2,mean))
stp_kmo = split.kms(gbm_mean, os=post_data_cox$BCR_FreeTime, oevn=post_data_cox$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(post_data_cox$BCR_FreeTime, post_data_cox$BCR_Event)~stp_temp_var, data=post_data_cox)
ggsurvplot(stp_km.split,title="BCR Free Survival for GBM Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))


#pre

#pre rp


#stepwise pre -rp
median(train_set)
std.error(train_set)
median(test_set)

train_set = numeric(50)
for(i in 1:50){
  train_set[i]=pre_op_cv_conc[i]
}


test_set = numeric(50)
for(i in 1:50){
  test_set[i]=pre_op_validation_conc[i]
}
std.error(test_set)

#lasso
median(lasso.cin.train.pre)
std.error(lasso.cin.train.pre)

median(lasso.cin.test.pre)
std.error(lasso.cin.test.pre)



#ridge
median(ridge.cin.train.pre)
std.error(ridge.cin.train.pre)

median(ridge.cin.test.pre)
std.error(ridge.cin.test.pre)

#gbm
median(unlist(gbm.train.conc.pre))
std.error(unlist(gbm.train.conc.pre))

median(unlist(gbm.test.conc.pre))
std.error(unlist(gbm.test.conc.pre))

#rfsrc
median(unlist(rsf_conc_train_pre))
std.error(unlist(rsf_conc_train_pre))

median(unlist(rsf_conc_test_pre))
std.error(unlist(rsf_conc_test_pre))


boxplot(train_set,lasso.cin.train.pre,ridge.cin.train.pre,unlist(gbm.train.conc.pre),unlist(rsf_conc_train_pre),names=c("Stepwise","Lasso","Ridge", "GBM","RSF"),
        col=c(0,2,7,4,5),main="Train Concordances")


boxplot(test_set,lasso.cin.test.pre,ridge.cin.test.pre,unlist(gbm.test.conc.pre),unlist(rsf_conc_test_pre),names=c("Stepwise","Lasso","Ridge", "GBM","RSF"),
        col=c(0,2,7,4,5),main="Test Concordances")


boxplot(train_set,test_set,lasso.cin.train.pre,lasso.cin.test.pre,ridge.cin.train.pre,ridge.cin.test.pre,unlist(gbm.train.conc.pre),unlist(gbm.test.conc.pre),unlist(rsf_conc_train_full),unlist(rsf_conc_test_full),names=c("Step Train","Step Test","Lasso Train","Lasso Test","Ridge Train","Ridge Test","GBM Train","GBM Test","RSF Train","RSF Test"),col=c(0,0,2,2,7,7,4,4,5,5))


#kmcurve pre

stp_mean = as.array(apply(pred_matrx_pre_data_cox,2,mean))
stp_kmo = split.kms(stp_mean, os=pre_data_cox$BCR_FreeTime, oevn=pre_data_cox$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(pre_data_cox$BCR_FreeTime, pre_data_cox$BCR_Event)~stp_temp_var, data=pre_data_cox)
ggsurvplot(stp_km.split,title="BCR Free Survival for Stepwise Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))



lasso_mean = as.array(apply(pred_matrx_lasso_test_pre,2,mean))
stp_kmo = split.kms(lasso_mean, os=pre_data_cox$BCR_FreeTime, oevn=pre_data_cox$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(pre_data_cox$BCR_FreeTime, pre_data_cox$BCR_Event)~stp_temp_var, data=pre_data_cox)
ggsurvplot(stp_km.split,title="BCR Free Survival for Lasso Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))



ridge_mean = as.array(apply(pred_matrx_ridge_test_pre,2,mean))
stp_kmo = split.kms(ridge_mean, os=pre_data_cox$BCR_FreeTime, oevn=pre_data_cox$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(pre_data_cox$BCR_FreeTime, pre_data_cox$BCR_Event)~stp_temp_var, data=pre_data_cox)
ggsurvplot(stp_km.split,title="BCR Free Survival for Ridge Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))


gbm_mean = as.array(apply(pred_matrx_gbm_test_pre,2,mean))
stp_kmo = split.kms(gbm_mean, os=pre_data_cox$BCR_FreeTime, oevn=pre_data_cox$BCR_Event)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(pre_data_cox$BCR_FreeTime, pre_data_cox$BCR_Event)~stp_temp_var, data=pre_data_cox)
ggsurvplot(stp_km.split,title="BCR Free Survival for GBM Cox model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))

