rm(list=ls())
library(tidyverse)
library(gridExtra)
library(glmnet)
library(boot)
library(HDDAData)

source('auc_cv_glm.r')
source('get_data.r')
# 1. Get the data
all_data <- get_data(seed=42,training_pct = .7)
trainX <- all_data$tr_x
trainY <- all_data$tr_y
testX <- all_data$test_x
testY <- all_data$test_y

# 2. Get the candidate models using crossvalidation
# 2a model_pcr1 PCR for the first k principal components
### new dataframe with reject_status and principal components
pca_x <- prcomp(trainX)
Z <- pca_x$x
fit_data <- data.frame(Y=trainY, Z)

### 10-fold Cross-validation on this one particular model, using AUC
### 2 step-procedure to detect k: big steps (10) first
set.seed(42)
trials <- seq(1,nrow(fit_data),10)
folds <- 10
cvpcr_run1 <- auc_cv_glm(trials, fit_data,folds)
nr_pc_opt <- cvpcr_run1 %>% filter(auc==max(auc)) %>% .$nr_of_pc
ggplot(cvpcr_run1)+
  geom_line(mapping=aes(x=nr_of_pc,y=auc),col='blue')+
  geom_vline(mapping=aes(xintercept = nr_pc_opt),col='firebrick')+
  labs(title ='cross validation pcr: steps of 10')+
  theme_bw()
## largest auc for nr_of_pc <= 20

### zoom in - same seed
set.seed(42)
trials <- seq(1,20)
folds <- 10
cvpcr_run2 <- auc_cv_glm(trials, fit_data,folds)
nr_pc_opt <- cvpcr_run2 %>% filter(auc==max(auc)) %>% .$nr_of_pc
ggplot(auc_run2)+
  geom_line(mapping=aes(x=nr_of_pc,y=auc),col='blue')+
  geom_vline(mapping=aes(xintercept = nr_pc_opt),col='firebrick')+
  labs(title ='cross validation pcr: steps of 1')+
  theme_bw()

# nr_pcr_opt = 12, but was very local
# 4 pc seems enough (auc .88)
nr_pc_chosen <- 4
model_pcr1 <- glm(Y ~ ., data = fit_data[1:(nr_pc_chosen+1)],family = 'binomial')
auc_pcr1 <- cvpcr_run2$auc[nr_pc_chosen]

# 3. PCR from glmnet
## elastic net cross validation on principal components
## alpha > makes number principal components sparse
## retains principal components relevant for Y, instead of 1:k
set.seed(42)
pca_x <- svd(trainX)
scores_x <- pca_x$u %*% diag(pca_x$d)
cverr <- cv.glmnet(x=scores_x,y=trainY, family="binomial",
                   alpha=.5,  nfolds= 10, type.measure = "auc")
lambda_cv <- cverr$lambda.1se
auc_pcr2 <- max(cverr$cvm)
model_pcr2 <- glmnet(x=scores_x,y=trainY, family="binomial",alpha=.5,lambda= lambda_cv)

coef_pcr2 <- coef(model_pcr2)
id_non_null <- which(abs(coef_pcr2)>0)
pcr2_relevant <-  rownames(coef_pcr2)[id_non_null[-1]]

# relevant principal components, V1, V2, V3, V10, auc = .87

# 4. ridge
set.seed(42)
folds <- 10
grid <- seq(1, 1000, by = 10)  # 1 to 1000 with steps of 10
ridge_mod_grid <- cv.glmnet(as.matrix(trainX), as.matrix(trainY), 
                            alpha = 0, lambda = grid,nfolds=folds, cost='auc', 
                            family="binomial")
lambda_cv <- ridge_mod_grid$lambda.1se
id_lambda_cv <- which(ridge_mod_grid$lambda == lambda_cv)
auc_ridge <- ridge_mod_grid$cvm[id_lambda_cv]

model_ridge <- glmnet(as.matrix(trainX), as.matrix(trainY), alpha = 0, 
                      lambda = lambda_cv, family="binomial")


# lasso
set.seed(42)
lasso_mod_grid <- cv.glmnet(as.matrix(trainX), as.matrix(trainY),
                            alpha = 1, nfolds=folds, cost='auc', family="binomial")
lambda_cv <- lasso_mod_grid$lambda.1se
id_lambda_cv <- which(lasso_mod_grid$lambda == lambda_cv)
auc_lasso <- lasso_mod_grid$cvm[id_lambda_cv]

model_lasso = glmnet(as.matrix(trainX), as.matrix(trainY),
                     alpha = 1,lambda=lambda_cv, family="binomial")

auc_models <- c(auc_pcr1,auc_pcr2,auc_ridge,auc_lasso)
names(auc_models) <- c('pcr1','pcr2','ridge','lasso')
# pcr and lasso perform best on trainingset,
# lasso is easier to interpret

# evaluate models on testset
# for pcr-models we use the z-scores for the observations
V <- pca_x$v
scores_test <- as.matrix(testX) %*% V
df_scores_test <- as.data.frame(scores_test)
names(df_scores_test) <- gsub('V','PC',names(df_scores_test))
eval_pcr1 <- predict(model_pcr1, newdata= df_scores_test)
eval_pcr2 <- predict(model_pcr2, newx = scores_test)
eval_ridge <- predict(model_ridge,newx = as.matrix(testX))
eval_lasso <- predict(model_lasso,newx = as.matrix(testX))
y_test <- factor(testY, levels=c(0,1),  labels=c('ok','reject'))
obs <- names(eval_pcr1)
for (mods in c('pcr1','pcr2','ridge','lasso')){
  yhat <- as.numeric(eval(as.name(paste0('eval_',mods))))
  data_val <- data.frame(Obs= obs, Yhat = yhat, Y = y_test)
  assign(paste0('plot_',mods), 
         ggplot(data=data_val, groups = Y)+
           geom_density(mapping=aes(x=Yhat, color=Y,fill=Y), alpha=.3)+
           labs(x= paste0('model ',mods))+
           theme_bw())
}
# densities
grid.arrange(plot_pcr1,plot_pcr2,
             plot_ridge,plot_lasso, 
             ncol=2)

# i prefer lasso here
