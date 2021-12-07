#packages
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("statOmics/HDDAData")

library(tidyverse)
library(ggbiplot)
install.packages(c("glmnet", "pls", "NormalBetaPrime", "boot"))
library(NormalBetaPrime)

library(glmnet)
library(pls)
library(MASS)
library(boot)
library(gridExtra)
library(pROC)

###########_0_DATA_PREPARATION_#############

library(HDDAData)
data("Einecke2010Kidney")

X <- scale(Einecke2010Kidney[, -1])
Y <- as.factor(Einecke2010Kidney[, 1])
table(Y)

n <- nrow(X)
H <- diag(n)-1/n*matrix(1,ncol=n,nrow=n)
X <- H%*%X



###########_1_DATA_EXPLORATION_#############

###_1.1_Plot on first 2 variables_###

#First we plot the observations to  the two first genes to see if there is a noticeable difference of
#transplant acceptation or rejection. This is not the case, as the plot illustrates no specific trend
Einecke2010Kidney$rejection_fac <- as.factor(ifelse(Einecke2010Kidney$Reject_Status == 0, "Accepted", "Rejected"))
ggplot(Einecke2010Kidney) +
  geom_point(aes(X211352_s_at, X239275_at, col = rejection_fac)) +
  ggtitle("Visualizing data")


###_1.2_PCA_###

svd_X <- svd(X)
Z <- svd_X$u %*% diag(svd_X$d) # Calculate the scores
V <- svd_X$v                 # Calculate the loadings
pca_x <- prcomp(X, center = FALSE, scale. = FALSE)

par(pch = 19, mfrow = c(1, 2))
plot(svd_X$d, type = "b", ylab = "Singular values", xlab = "PCs")

var_explained <- svd_X$d^2 / sum(svd_X$d^2)
plot(var_explained,
     type = "b", ylab = "Percent variance explained", xlab = "PCs",
     col = 2
)

#plot data on first 2 PC's. Does not separate clearly between accepted and rejected cases. 
#The rejected are more present in the bottom left corner, the accepted more in top right corner.
par(mfrow = c(1, 1))
cols <- c("0" = "red", "1" = "blue")
plot(Z[, 1], Z[, 2],
     col = cols[Y],
     xlab = "PC1", ylab = "PC2", pch = 19
)
legend("bottomright", c("Accepted", "Rejected"),
       col = c("red", "blue"),
       pch = 19, title = "Kidney"
)

#We see that only a small amount of genes are driving the PC's as most of the genes has loadings near 0
#Therefore, we think sparse PCA is worthwhile
par(mfrow = c(1, 1))
hist(V[, 1], breaks = 50, xlab = "PC 1 loadings", main = "")
abline(v = c(
  quantile(V[, 1], 0.05),
  quantile(V[, 1], 0.95)
), col = "red", lwd = 2)

hist(V[, 2], breaks = 50, xlab = "PC 2 loadings", main = "")
abline(v = c(
  quantile(V[, 2], 0.05),
  quantile(V[, 2], 0.95)
), col = "red", lwd = 2)

## other way to look at this
par(mfrow=c(1,1))
TOL <- 10^(-5) 
zero_load <- apply(svd_X$v,2,function(column){
  sum(column< TOL)})
hist(zero_load, main = "Histogram for loadings", 
     xlab = "number of irrelevant genes: value < 1e-5")

###_1.3_Sparse_PCA_###

# For PC1
set.seed(45)
fit_loadings1 <- cv.glmnet(X, Z[, 1],
                           alpha = 0.5, nfolds = 5)
plot(fit_loadings1, main = "PC1")

# For PC2
set.seed(45)
fit_loadings2 <- cv.glmnet(X, Z[, 2], alpha = 0.5, nfolds = 5)
plot(fit_loadings2, main = "PC2")

par(mfrow = c(1, 1))
plot(fit_loadings1$glmnet.fit, main = "PC1", xvar = "lambda")
abline(v = log(fit_loadings1$lambda.min), lty = 3)
abline(v = log(fit_loadings1$lambda.1se), lty = 3)
plot(fit_loadings2$glmnet.fit, main = "PC2", xvar = "lambda")
abline(v = log(fit_loadings2$lambda.min), lty = 3)
abline(v = log(fit_loadings2$lambda.1se), lty = 3)

# We see that for PC1 201 to 217 genes are the most important and 198 to 205 for PC2.
fit_loadings1
fit_loadings2


sparse_loadings1 <- as.vector(coef(fit_loadings1, s = fit_loadings1$lambda.1se))
sparse_loadings2 <- as.vector(coef(fit_loadings2, s = fit_loadings2$lambda.1se))

## How many non-zero loadings do we have (excluding the intercept)?

(non_zero1 <- sum(abs(sparse_loadings1[-1]) > 0))
(non_zero2 <- sum(abs(sparse_loadings2[-1]) > 0))

SPC1 <- X %*% sparse_loadings1[-1] # without the intercept
SPC2 <- X %*% sparse_loadings2[-1] # without the intercept


#What we see is that 201 genes for PC1 and 198 genes for PC2 are useful.
par(mfrow = c(1, 2))
plot(Z[, 1], Z[, 2],
     col = cols[Y], xlab = "PC1", ylab = "PC2", pch = 16,
     main = "All genes \nfor PC1 and PC2"
)
plot(SPC1, SPC2,
     col = cols[Y], xlab = "SPC1", ylab = "SPC2", pch = 16,
     main = paste(non_zero1, "genes for SPC1 \n and", non_zero2, "genes for SPC2")
)


###_1.4_LDA_###

#PCA only considers the genes in the decomposition and not the class memebership (Accepted (0) or rejected (1))
#Therefore, LDA can be used to get a better insight
library(HDDAData)
data("Einecke2010Kidney")

X <- scale(Einecke2010Kidney[, -1])
Y <- as.factor(Einecke2010Kidney[, 1])

kidney.lda <- MASS::lda(x = X, grouping = Y)

Vlda <- kidney.lda$scaling
colnames(Vlda) <- paste0("V",1:ncol(Vlda))
Zlda <- X%*%Vlda
colnames(Zlda) <- paste0("Z",1:ncol(Zlda))

par(mfrow = c(1, 1))
boxplot(Zlda ~ Y, col = cols, ylab = expression("Z"[1]),
        main = "Separation of accepted and rejected samples by LDA")

#Sparse LDA
set.seed(45)
lda_loadings <- cv.glmnet(X, Zlda, alpha = 0.5, nfolds = 5)
plot(lda_loadings)

sparse_lda_loadings <- as.vector(
  coef(lda_loadings, s = lda_loadings$lambda.1se)
)

SLDA <- X %*% sparse_lda_loadings[-1]

# number of non-zero loadings
n_nonzero <- sum(sparse_lda_loadings != 0)
boxplot(SLDA ~ Y,
        col = cols, ylab = "SLDA",
        main = sprintf("Subset of %d genes", n_nonzero)
)



###########_2_Hypothesis_Testing_#############

gene_data <- as.matrix(Einecke2010Kidney[, -1])
group <- Einecke2010Kidney[, 1]

ttest_results <- t(apply(gene_data, 2, function(x) {
  t_test <- t.test(x ~ group)
  p_val <- t_test$p.value
  stat <- t_test$statistic
  df <- t_test$parameter
  ## Return values in named vector
  c(stat, "p_val" = p_val, df)
}))

head(ttest_results)

p_vals <- ttest_results[, "p_val"]
hist(
  p_vals,
  breaks = seq(0, 1, by = 0.05), main = "", xlab = "p-value",
  ylim = c(0, 4000)
)

#2883 significant genes
alpha <- 0.05
sum(p_vals < alpha)

fdr <- p.adjust(p_vals, method = "BH")
fdr
plot(
  p_vals[order(p_vals)], fdr[order(p_vals)],
  pch = 19, cex = 0.6, xlab = "p-value", ylab = "FDR-adjusted p-value", col = 4
)
abline(a = 0, b = 1)

#2227 significant genes
sum(fdr < 0.10)

#The 20 most significant genes
head(p_vals[fdr < 0.10], 20)



###########_3_Rejection_Prediction_#############

X <- scale(Einecke2010Kidney[, -1], center = TRUE, scale = TRUE)
Y <- as.factor(Einecke2010Kidney[, 1])

set.seed(1)
trainID <- sample(nrow(X), 0.7*nrow(X))
trainX <- X[trainID, ]
trainY <- Y[trainID]
testX <- X[-trainID, ]
testY <- Y[-trainID]
train_data <- data.frame(trainY, trainX)
test_data <- data.frame(testY, testX)


###_Number_of_PCs_###

## Calculate PCA and extract scores
pca_X <- prcomp(trainX)
Z <- pca_X$x
head(Z)

n_PC <- ncol(Z)

fit_data <- data.frame(trainY, Z)
head(fit_data)

full_model <- glm(trainY ~ ., data = fit_data, family = "binomial")
summary(full_model)

full_model_cv <- cv.glm(
  data = fit_data,  glmfit = full_model,
  cost = auc, K = 4  # note: specify the auc function (from pROC) without`()`!
)

full_model_cv$delta[1]

cv_auc <- vector("numeric", length = n_PC)
set.seed(12) # seed for reproducibility
for (i in seq_len(n_PC)) {
  ## Prepare fit_data; subset number of PCs to i
  fit_data <- data.frame(trainY, Z[, 1:i, drop = FALSE])  # use drop = FALSE to avoid problems when subsetting single column
  pcr_mod <- suppressWarnings(
    glm(trainY ~ ., data = fit_data, family = "binomial")
  )
  
  ## Do 4-fold CV while suppressing Warnings and Messages
  cv <- suppressWarnings(
    suppressMessages(
      cv.glm(fit_data, pcr_mod, cost = pROC::auc, K = 4)
    )
  )
  cv_auc[i] <- cv$delta[1]
}
names(cv_auc) <- seq_along(cv_auc)
cv_auc  
optim_nPC <- names(which.max(cv_auc))
optim_nPC
cv_auc[optim_nPC]



plot(names(cv_auc), cv_auc, xlab = "n PCs", ylab = "AUC", type = "l")
abline(v = optim_nPC, col = "red")
#optimal nr of PC's = 17 for me
#auc = 87.17%
#decision on what k should be

l <- 17
pca <- prcomp(trainX)
Vk <- pca$rotation[, 1:l] # the loadings matrix
Zk <- pca$x[, 1:l]

pcr_model1 <- glm(trainY ~ Zk, family = "binomial")
summary(pcr_model1)

################################

###_Ridge_###

mRidge <- glmnet(
  x = trainX,
  y = trainY,
  alpha = 0,
  family="binomial")  



set.seed(123)
k <- 15
mCvRidge <- cv.glmnet(
  x = trainX,
  y = trainY,
  alpha = 0, nfolds = k,
  type.measure = "auc",
  family = "binomial")

mCvRidge
#lambda 159.6 / 11.8??
#auc = 83.40%/ 86.33%

################################

###_Lasso_###

mLasso <- glmnet(
  x = trainX,
  y = trainY,
  alpha = 1,
  family="binomial")  # lasso: alpha = 1

plot(mLasso, xvar = "lambda", xlim = c(-6,-1.5))

set.seed(123)
k <- 15
mCvLasso <- cv.glmnet(
  x = trainX,
  y = trainY,
  alpha = 1, nfolds = k,
  type.measure = "auc",
  family = "binomial")

mCvLasso
#lambda is very small??
#auc = 82.74%/ 84.63%

################################

