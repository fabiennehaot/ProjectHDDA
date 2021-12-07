---
title: "HDDA Project"
author: "L. Barbier, K. De Witte, F. Haot, K. Putseys"
date: "12/7/2021"
output:   
  html_document:
    keep_md: true
---



# Research Question

# Executive Summary

# Technical Report

## Data

```r
data("Einecke2010Kidney")
X_raw <- Einecke2010Kidney[, -1]
Y <- factor(Einecke2010Kidney[, 1], levels = c(0,1),labels = c('accept','reject'))
```

## Exploration of the Data

## Hypothesis Testing 

```r
gene_data <- as.matrix(X_raw)
group <- Y
# non adjusted p-values as named vector
ttest_results <- apply(gene_data, 2, function(x) {
  t_test <- t.test(x ~ group)
  p_val <- t_test$p.value
  stat <- t_test$statistic
  df <- t_test$parameter
  z_val <- case_when(stat < 0 ~ qnorm(p_val/2), TRUE~ qnorm(1-p_val/2))
  tibble(stat= stat, p_val = p_val,z_val = z_val, df = df)})
t_stats <- bind_rows(ttest_results) %>%
  mutate(
    gene = colnames(X_raw),
    relevant = cut(z_val,
                   c(-Inf,qnorm(.025),qnorm(.975),Inf), 
                   labels= c('low','mid','high')))


ggplot(data=as.data.frame(t_stats)) + 
  geom_point(mapping=aes(x=seq_along(stat),y=stat,
                         color = relevant),size=.3)+
  labs(title = "t-statistics", x= "gene",y= "unadjusted t-statistic",
       caption= "There are more observations in the lower tail") +
  theme_bw()
```

![](Report-project_files/figure-html/hyptest_unadjusted-1.png)<!-- -->

```r
table_p_nonadj <- table(t_stats$p_val < 0.05)
no_corr_genes <- t_stats %>% 
  filter(p_val<.05) %>%
  .$gene
```

There are 2883 out of 10.000 genes with p-value of 0.05.


```r
#fdr analysis - current behaviour of p-values: uniform?
t_stats %>%
  ggplot(aes(x = p_val)) +
  geom_histogram(fill = "firebrick",breaks = seq(0,1,.05))+
  labs(title = 'histogram of p-values')+
  theme_bw()
```

![](Report-project_files/figure-html/hyptest_fdr_hist-1.png)<!-- -->

This histogram shows a distribution which is close to a uniform distribution
For the larger p-values, but with more small p-values than expected under 
a uniform distribution. The higher number of genes with small p-values indicates that some genes show a different behaviour regarding kidney rejection



```r
# adjusted p values
fdr_level <- .1
t_stats <- t_stats %>%
  mutate(p_adj = p.adjust(p_val, method="fdr"))

# monotonous transformation:
pp_adj <- t_stats %>%
  ggplot(aes(x=p_val,y=p_adj)) +
  geom_point(size = .3,color='blue') +
  geom_segment(x=0,y=0,xend=1,yend=1) +
  labs(y= "adjusted p-value (BH, 1995)") +
  theme_bw()

# right plot removes observations with high fdr 
# 8090 obs are removed when ylim .05
# 7575 when ylim .1
# plotting seems to take a while
grid.arrange(pp_adj,
             pp_adj + ylim(c(0,fdr_level)),
             ncol=2)
```

```
## Warning: Removed 7773 rows containing missing values (geom_point).
```

![](Report-project_files/figure-html/hyptest_fdr_padj-1.png)<!-- -->

```r
# these are the genes
table_p_adj <- table(t_stats$p_adj < fdr_level)
fdr_genes <- t_stats %>%
  filter(p_adj < fdr_level)%>%
  .$gene
```

After selection still high dimensional. The genes tell us something about the rejection status through the t-tests.


```r
#local fdr

fdr_x <- locfdr(t_stats$z_val,plot=3) 
```

![](Report-project_files/figure-html/localfdr-1.png)<!-- -->
 The probability $Pr[fdr(Z) < \alpha]$ is the probability that a non-null 
 can be detected when the nominal local fdr is set at $\alpha$
 $\alpha$ = O.1 corresponds with proportion 60%.
 Efdr = .215 means that the local fdr for a typical non-null feature is expected to be 21.5%.


 $\delta = -0.605$ indicates that the null distribution is shifted to negative values due to correlations or non-compliance with the null-hypotheses.
 $\sigma = 1.621$: standard deviation becomes larger than 1
The fact that the blue and green line differ on the left side show the signal genes for kidney rejection. We don't expect to find signal genes with high t-values.


```r
lfdr_level <- fdr_x$Efdr[1]
t_stats <- t_stats %>%
  mutate(
    lfdr = fdr_x$fdr,
    zfdr = (lfdr < lfdr_level) * z_val)

t_stats %>%
  filter(lfdr < lfdr_level) %>%
  summarize(nr_of_genes=n(), mean_fdr=mean(lfdr))
```

```
## # A tibble: 1 × 2
##   nr_of_genes mean_fdr
##         <int>    <dbl>
## 1         358   0.0372
```

## Model Selection

## Model Evalutaion

## Conclusions


## Including Plots

