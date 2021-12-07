remotes::install_github("statOmics/HDDAData")
library(HDDAData)
library(tidyverse)
library(patchwork)

data("Einecke2010Kidney")
# dim
dim(Einecke2010Kidney)
str(Einecke2010Kidney)
n_col <- length(Einecke2010Kidney)
# matrices
Y <- Einecke2010Kidney$Reject_Status
X <- Einecke2010Kidney[,-1]

# soundness
nr_values  <- sapply(X, function(var) {length(var)})
nr_missing <- sapply(X, function(var) {length(which(is.na(var)))})
nr_unique  <- sapply(X, function(var) {n_distinct(var)})
range(nr_values)
range(nr_unique)
range(nr_missing)
# no missing values, no concentrations

# outliers?
# scaling: all variables are intensity parameters
max_per_var <- sapply(Einecke2010Kidney, max)
min_per_var <- sapply(Einecke2010Kidney, min)
mean_per_var <- sapply(Einecke2010Kidney, mean)
tb_summary <- tibble(Name= colnames(Einecke2010Kidney),
                     Nr = seq_along(Einecke2010Kidney),
                     MaxPerVar = max_per_var,
                     MinPerVar = min_per_var,
                     MeanPerVar=mean_per_var)
ggplot(data=tb_summary)+
  geom_point(mapping=aes(x= Nr,y=MinPerVar),col='blue',size=.3)+
  geom_point(mapping=aes(x= Nr,y=MeanPerVar),col='green',size=.3)+
  geom_point(mapping=aes(x= Nr,y=MaxPerVar),col='red',size=.3)+
  theme_bw()
# data looks calibrated, close to centered but not scaled
# differences in variability are not extreme 
# no outliers
# let's scale X, Y is not useful (categorical)
X <- scale(X)

# do principal components match?
# svd
Svd_x <- svd(X)
ggplot()+
  geom_col(mapping= aes(x=seq_along(Svd_x$d),y=Svd_x$d), fill='lightsteelblue')+
  labs(x='right eigenvectors', y = 'singular value', title= 'SVD for X')+
  theme_bw()
tb_singular_values <- tibble(SingularValue= Svd_x$d)
tb_singular_values <- tb_singular_values %>%
  mutate(Nr = seq_along(SingularValue),
         VarianceExplained = SingularValue^2/ sum(SingularValue^2),
         VarianceCumul = cumsum(VarianceExplained))
             
ggplot(data=tb_singular_values)+
  geom_col(mapping= aes(x=Nr,  y=VarianceExplained), fill='lightsteelblue')+
  geom_line(mapping= aes(x=Nr, y=VarianceCumul),col='blue')+
  labs(x='PCA vectors', y = 'Explained Variance', 
       title= 'Dimensions needed to explain X')+
  theme_bw()
# we can explain 63% of X with 50 first PCA's
# we can explain 19% of X with 2 first PCA's
# biplot(prcomp(X)) not very interesting (only 19% and confusing to watch)
# idea of biplot is determining inproduct (e_k, pca_i) for i=1:dim_red, 
# k = nr of variable that interests us



# is rejection associated with pca1 or pca2
Scores <- X %*% Svd_x$v
tb_scores <- tibble(PCA1 = Scores[,1],
                    PCA2 = Scores[,2],
                    Object=rownames(Scores),
                    Rejected = as.factor(Y))
                    
ggplot(tb_scores)+
  geom_point(mapping=aes(x= PCA1,y=PCA2, col=Rejected))+
  labs(title="Scores for 2 PCA's")+
  theme_bw()  
# not really distinguishing


# tracking of the most outspoken genes
# highest ttests
lst_tests <- lapply(colnames(X),function(nm){
       test <- t.test(X[,nm]~Y)
       tibble(gene= nm, tvalue = as.numeric(test$statistic),
              pvalue = as.numeric(test$p.value),
              df = as.numeric(test$parameter))})
tb_tests <- bind_rows(lst_tests)
#highest 2
names_top5 <- tb_tests %>% 
  arrange(pvalue) %>% 
  head(5) %>% 
  .$gene
nm1 <- names_top5[1] 
nm2 <- names_top5[2] 
tb2 <- tibble(gene1 = X[[nm1]],
              gene2 = X[[nm2]],
              y = as.factor(Y)) 
plot1 <- ggplot(data= tb2, groups = y) +
  geom_histogram(mapping=aes(x=gene1,fill = y),
                 position = "dodge",bins=20)+
  theme_bw()
plot2 <- ggplot(data= tb2, groups = y) +
  geom_histogram(mapping=aes(x=gene2,fill = y),
                 position = "dodge", bins=20)+
  theme_bw()
plot1+plot2