# Project assignment High dimensional data analysis
# 18 nov 2021

library(HDDAData)
data("Einecke2010Kidney")

# Data dimensions 
dim(Einecke2010Kidney)
# 250 10001

# Showing first 6 rows and 10 columns
head(Einecke2010Kidney[, 1:10])

# Extract gene expression data as matrix X
X <- as.matrix(Einecke2010Kidney[, -1])
dim(X)
# 250 10000
str(X)

# Extract Reject_Status column as vector
reject_status <- Einecke2010Kidney$Reject_Status
names(reject_status) <- rownames(Einecke2010Kidney)
table(reject_status) # number of 0's (accepts) and 1's (rejects)
#> reject_status
#>   0   1 
#> 183  67

######################
## data exploration  #
######################

library(tidyverse)
# range of the means
X %>% colMeans %>%range
# so not yet centered

# looking for missing values
which(is.na(X))
# no missing values

# singular value decomposition
svdX <- svd(X)
k <- 2
Zk <- svdX$u[,1:k] %*% diag(svdX$d[1:k])
colnames(Zk) <- paste0("Z",1:k)
Vk <- svdX$v[,1:k]

Zk %>%
  as.data.frame %>%
  mutate(reject = reject_status %>% as.factor) %>%
  ggplot(aes(x= Z1, y = Z2, color = reject)) +
  geom_point(size = 3)
# no seperation to be seen
# using only the first 2 PCs does not seem to separate the rejection and accepted cases clearly.

grid.arrange(
  Vk %>%
    as.data.frame %>%
    mutate(geneID = 1:nrow(Vk)) %>%
    ggplot(aes(x = geneID, y = V1)) +
    geom_point(pch=21) +
    geom_hline(yintercept = c(-2,0,2)*sd(Vk[,1]), col = "red") ,
  Vk %>%
    as.data.frame %>%
    mutate(geneID = 1:nrow(Vk)) %>%
    ggplot(aes(x = geneID, y = V2)) +
    geom_point(pch=21) +
    geom_hline(yintercept = c(-2,0,2)*sd(Vk[,2]), col = "red"),
  ncol=2)
# almost impossible to interpret genes
# no mega large outliers

# other pca things
Z <- svdX$u %*% diag(svdX$d) # Calculate the scores
V <- svdX$v                 # Calculate the loadings

# Plotting parameters
par(pch = 19, mfrow = c(1, 2))
plot(svdX$d, type = "b", ylab = "Singular values", xlab = "PCs")

# Percentage variance explained for each PC
var_explained <- svdX$d^2 / sum(svdX$d^2)
plot(var_explained,
     type = "b", ylab = "Percent variance explained", xlab = "PCs",
     col = 2
)
# a subset of the pcs already explains quite a bit of variance

# Plot histograms of the loadings of the first and second PCs. 
par(mfrow = c(2, 1))
# First
hist(V[, 1], breaks = 50, xlab = "PC 1 loadings", main = "")
# Add vertical line at 95% quantile
abline(v = c(quantile(V[,1],0.05), quantile(V[, 1], 0.95)), col = "red", lwd = 2)

# Second
hist(V[, 2], breaks = 50, xlab = "PC 2 loadings", main = "")
abline(v = c(
  quantile(V[, 2], 0.05),
  quantile(V[, 2], 0.95)
), col = "red", lwd = 2)
# remember that the PC loadings reflect the contributions of each feature (in this case: gene) to the PC. 
# From these histograms it should be clear that only a minor fraction of the genes are really driving
# these first 2 PCs, 

# you could do sparce pca
# but this already uses ridge regression?

library(MASS)
# Perform LDA
reject_lda <- lda(x = X, grouping = reject_status)
V1 <- reject_lda$scaling
Z1 <- X %*% V1
par(mfrow = c(1, 1))
boxplot(Z1 ~ reject_status, ylab = expression("Z"[1]),
        main = "Separation of rejected and accepted cases by LDA")


###################
# controlling FDR #
###################

data <- as.matrix(Einecke2010Kidney[, -1])
group <- Einecke2010Kidney$Reject_Status

# Use `apply` to loop over columns of gene data and perform t-tests, extracting
# p-values, test stastistics, and degrees of freedom
ttest_results <- t(apply(data, 2, function(x) {
  t_test <- t.test(x ~ group)
  p_val <- t_test$p.value
  stat <- t_test$statistic
  df <- t_test$parameter
  ## Return values in named vector
  c(stat, "p_val" = p_val, df)
}))

# Take a look at results
head(ttest_results)
genes <- rownames(ttest_results)
ttest_results <- as.data.frame(ttest_results)

ggplot(data=ttest_results) + 
  geom_point(mapping=aes(x=seq_along(genes),y=t),size=.3)+
  geom_hline(yintercept=qnorm(.975),col='red')+
  geom_hline(yintercept=-qnorm(.975),col='red')+
  theme_bw()

# number of discoveries when alpha = 0.05
p_vals <- ttest_results[, "p_val"]
alpha <- 0.05
sum(p_vals < alpha) # 3661 (alpha = 0.10) , 2883 (alpha = 0.05)
table(p_vals < alpha)

# FDR method to adjust the original p-values
fdr <- p.adjust(p_vals, method = "BH")
fdr <- p.adjust(p_vals, method = 'fdr')

plot(
  p_vals[order(p_vals)], fdr[order(p_vals)],
  pch = 19, cex = 0.6, xlab = "p-value", ylab = "FDR-adjusted p-value", col = 4
)
abline(a = 0, b = 1)

fdr_level <- 0.10
sum(fdr < fdr_level)
table(fdr < fdr_level)


fdr_genes <- filter(ttest_results, fdr < fdr_level)
rownames(fdr_genes)

## local fdr
library(locfdr)
z_val <- qnorm(ttest_results$p_val)
localfdr <- locfdr(z_val, nulltype=0)
localfdr<- locfdr(z_val,plot=2)
localfdr<- locfdr(z_val,plot=3) 

localfdr$Efdr  
# Efdr = .2598
# local fdr for a typical non-null feature is expected to be 26%

lfdr_level <- localfdr$Efdr[1]
localfdr$fdr

ttest_results <- ttest_results %>%
  mutate(
    lfdr = localfdr$fdr,
    zfdr = (lfdr < lfdr_level) * z_val)
head(ttest_results)

ttest_results %>%
  filter(lfdr < lfdr_level) %>%
  summarize(nr_of_genes=n(), mean_fdr=mean(lfdr))


###############################################################################

## regression ##

Xscale <- scale(X, center = TRUE, scale = TRUE)

set.seed(1)
# Sample 80 random IDs from the rows of X (120 total)
trainID <- sample(nrow(Xscale), 0.7*nrow(Xscale))

# Training data
trainX <- Xscale[trainID, ]
trainY <- reject_status[trainID]
table(trainY)

# Test data
testX <- Xscale[-trainID, ]
testY <- reject_status[-trainID]
table(testY)

train_data <- data.frame("reject" = trainY, trainX)
test_data <- data.frame("Treject" = testY, testX)

library(glmnet)

##################################
# principal component regression #
##################################

# Calculate PCA and extract scores
pca_X <- prcomp(trainX)
Z <- pca_X$x
head(Z)

# Total number of available PCs
n_PC <- ncol(Z)
n_PC

# cv.glm() requires the response and predictors in one data.frame, so we need
# to combine them back together
fit_data <- data.frame(trainY, Z)
head(fit_data)

# Example of PC Log. Reg. with all PCs
full_model <- glm(trainY ~ ., data = fit_data, family = "binomial")
summary(full_model)

library(boot)
library(pROC)

# 4-fold Cross-validation on this one particular model, using AUC and K = 4
set.seed(2)
full_model_cv <- cv.glm(
  data = fit_data,  glmfit = full_model,
  cost = auc, K = 4  # note: specify the auc function (from pROC) without`()`!
)

full_model_cv$delta[1] # This is the AUC for this particular model estimated by AUC

# Now we'll wrap this code in a for-loop and repeat for each number of PCs
cv_auc <- vector("numeric", length = n_PC)
set.seed(12) # seed for reproducibility
for (i in seq_len(n_PC)) {
  # Prepare fit_data; subset number of PCs to i
  fit_data <- data.frame(trainY, Z[, 1:i, drop = FALSE])  #use drop = FALSE to avoid problems when subsetting single column
  pcr_mod <- suppressWarnings(
    glm(trainY ~ ., data = fit_data, family = "binomial")
  )
  
  # Do 4-fold CV while suppressing Warnings and Messages
  cv <- suppressWarnings(
    suppressMessages(
      cv.glm(fit_data, pcr_mod, cost = pROC::auc, K = 4)
    )
  )
  cv_auc[i] <- cv$delta[1]
}
names(cv_auc) <- seq_along(cv_auc)
cv_auc

# Finding the optimal nr. of PCs corresponds to finding the max. AUC
optim_nPC <- names(which.max(cv_auc))
optim_nPC

plot(names(cv_auc), cv_auc, xlab = "n PCs", ylab = "AUC", type = "l")
abline(v = optim_nPC, col = "red")

cv_auc[optim_nPC]

# pcregression with chosen number of components
library(pls)

# X is already scaled and centered, so that's not needed.
pcr_model <- pcr(trainY ~ trainX, ncomp = optim_nPC, family = "binomial")
summary(pcr_model)

#set.seed(123)
#K <- 10
# The 'Y ~ .' notation means: fit Y by every other variable in the data
#pcr_cv2 <- pcr(trainY ~ ., data = train_data, validation = "CV", segments = K)
#summary(pcr_cv2)
#plot(pcr_cv2, plottype = "validation")
#optimal_ncomp <- selectNcomp(pcr_cv, method = "onesigma", plot = TRUE)

####################
# Lasso regression #
####################

# Lasso
mLasso <- glmnet(x = trainX, y = trainY, alpha = 1, family  = "binomial")  # alpha = 1 -> Lasso
plot(mLasso, xvar = "lambda", xlim = c(-6,-1.5))

# cross validation to find optimal lambda
mCvLasso <- cv.glmnet(x = trainX, y = trainY, alpha = 1, family  = "binomial", type.measure = "auc")
plot(mCvLasso)

# Lasso with optimal lambda
mLassoOpt <- glmnet(x = trainX, y = trainY, alpha = 1, family  = "binomial", lambda = mCvLasso$lambda.min)
summary(coef(mLassoOpt))
mCvLasso$lambda.min
# with optimal lambda (largest auc) the output shows the non-zero estimated regression coefficients
qplot(summary(coef(mLassoOpt))[-1,1],
      summary(coef(mLassoOpt))[-1,3]) + xlab("gene") + ylab("beta-hat") + geom_hline(yintercept = 0, color = "red")


mLassoOpt2 <- glmnet(x = trainX, y = trainY, alpha = 1, family  = "binomial", lambda = mCvLasso$lambda.1se)
summary(coef(mLassoOpt2))

# plot roc curve for Lasso
#BiocManager::install("plotROC")
library(plotROC)

dfLassoOpt <- data.frame(
  pi = predict(mCvLasso,
               newx = testX,
               s = mCvLasso$lambda.min,
               type = "response") %>% c(.),
  known.truth = testY)

roc <-
  dfLassoOpt  %>%
  ggplot(aes(d = known.truth, m = pi)) +
  geom_roc(n.cuts = 0) +
  xlab("1-specificity (FPR)") +
  ylab("sensitivity (TPR)")
roc

calc_auc(roc)

# lasso on test data
lasso_preds <- predict(mCvLasso, s = mCvLasso$lambda.min, newx = as.matrix(testX))

# Inputs:
# * obs: vector of test observations
# * pred: vector of model predictions
# * cutoff_values: vector of prediction thresholds c
calculate_misclass_error <- function(obs, pred, cutoff_values) {
  stopifnot(length(obs) == length(pred))
  misclass_errors <- rep(NA, length(cutoff_values))
  for (i in seq_along(cutoff_values)) {
    cutoff <- cutoff_values[i]
    ypred <- as.numeric(pred > cutoff) # translates TRUE/FALSE to 1/0
    misclass_errors[i] <- mean(ypred != obs) # proportion of misclassifications
  }
  data.frame(
    "cutoff" = cutoff_values,
    "misclass" = misclass_errors
  )
}
result <- calculate_misclass_error(testY, lasso_preds, seq(0.1,0.9,by=0.01))
id <- which.min(result[,2])
cutoff <- result[id, 1]
newY <- ifelse(lasso_preds > cutoff, 1, 0)  # everything bigger than cutoff -> predict 1
#newY

# calculate sensitivity and specificity
sens <- sum(newY == testY & testY==1)/ sum(testY == 1)
spec <- sum(newY == testY & testY==0)/ sum(testY == 0)
sens 
spec


#####################
# Ridge regression ##
#####################

# Rigde regression
mRidge <- glmnet(x = trainX, y = trainY, alpha = 0, family  = "binomial")
plot(mRidge, xvar = "lambda")

# Rigde regression with cross-validation
mCvRidge <- cv.glmnet(x = trainX, y = trainY, alpha = 0, family  = "binomial", type.measure = "auc")  
# alpha = 0 -> ridge
plot(mCvRidge)

# Ridge regression with optimal lambda
mRidgeOpt <- glmnet(x = trainX, y = trainY, alpha = 0, family  = "binomial", lambda = mCvRidge$lambda.min)
summary(coef(mRidgeOpt))
mCvRidge$lambda.min
# with optimal lambda (largest auc) the output shows the non-zero estimated regression coefficients
qplot(summary(coef(mRidgeOpt))[-1,1],
      summary(coef(mRidgeOpt))[-1,3]) + xlab("gene") + ylab("beta-hat") + geom_hline(yintercept = 0, color = "red")

# roc for ridge regression
dfRidgeOpt <- data.frame(
  pi = predict(mCvRidge,
               newx = testX,
               s = mCvRidge$lambda.min,
               type = "response") %>% c(.),
  known.truth = testY)

roc <-
  dfRidgeOpt  %>%
  ggplot(aes(d = known.truth, m = pi)) +
  geom_roc(n.cuts = 0) +
  xlab("1-specificity (FPR)") +
  ylab("sensitivity (TPR)")
roc

calc_auc(roc)

# ridge on test data
ridge_preds <- predict(mCvRidge, s = mCvRidge$lambda.min, newx = as.matrix(testX))

result <- calculate_misclass_error(testY, ridge_preds, seq(0.1,0.9,by=0.01))
id <- which.min(result[,2])
cutoff <- result[id, 1]
cutoff
newY <- ifelse(ridge_preds > cutoff, 1, 0)  # everything bigger than cutoff -> predict 1

# calculate sensitivity and specificity
sens <- sum(newY == testY & testY==1)/ sum(testY == 1)
spec <- sum(newY == testY & testY==0)/ sum(testY == 0)
sens 
spec

