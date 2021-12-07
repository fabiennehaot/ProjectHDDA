# libraries
library(HDDAData)
library(tidyverse)
library(gridExtra)
library(locfdr)

# data
data("Einecke2010Kidney")
Y <- Einecke2010Kidney$Reject_Status
X <- Einecke2010Kidney[,-1]

# tests per gene: multiple testing problem
df_t <- nrow(X)-2
t_vals <- sapply(names(X), function(nm){ t.test(X[[nm]]~Y,alternative="two.sided")}$statistic)
p_vals <- pt(t_vals,df=df_t)
z_vals <- qnorm(p_vals)
t_stats <- tibble(gene = names(t_vals),
                  t_values = as.numeric(t_vals),
                  p_values = as.numeric(p_vals),
                  z_values = as.numeric(z_vals))
ggplot(data=t_stats) + 
  geom_point(mapping=aes(x=seq_along(gene),y=z_values),size=.3)+
  geom_hline(yintercept=qnorm(.975),col='red')+
  geom_hline(yintercept=-qnorm(.975),col='red')+
  theme_bw()

table(t_stats$p_values < 0.05)
no_corr_genes <- t_stats %>% 
  filter(abs(z_values)>qnorm(.975))%>%
  .$gene

#controlling the fdr: short version
fdr_level <- .1
t_stats <- t_stats %>% 
  mutate(p_adj = p.adjust(p_values,'fdr'))
table(t_stats$p_adj < fdr_level)


fdr_genes <- t_stats %>%
  filter(p_adj < fdr_level)%>%
  .$gene

# fdr analysis longer version
# Histogram should show a distribution which is close to a uniform distribution
# for the larger p-values, but with more small p-values than expected under 
# a uniform distribution
fdr_level <- .1
t_stats %>%
  ggplot(aes(x = p_values)) +
  geom_histogram(fill = "firebrick",breaks = seq(0,1,.05))+
  labs(title = 'histogram of p-values')+
  theme_bw()
# too many genes with p-value 1, but there is a signal - uniform split visible


# adjusted p values
t_stats <- t_stats %>%
  mutate(
    p_adj = p.adjust(p_values, method="fdr"),
    z_fdr = (p_adj < fdr_level) * z_values)

pp_adj <- t_stats %>%
  ggplot(aes(p_values,p_adj)) +
  geom_point(size = .3,color='blue') +
  geom_segment(x=0,y=0,xend=1,yend=1) +
  ylab("adjusted p-value (BH, 1995)")+
  theme_bw()

# right plot removes observations with high fdr 
# 8090 obs are removed when ylim .05
# 7575 when ylim .1
# plotting seems to take a while
grid.arrange(pp_adj,
             pp_adj + ylim(c(0,fdr_level)),
             ncol=2)

# BH corrected p-values
table(t_stats$p_adj < fdr_level)

# after selection still high dimensional
# we can perform pcr / ridge / lasso 

#local fdr
fdr_x <- locfdr(z_vals,plot=3) 
# The probability Pf1[fdr(Z)<alpha] is the probability that a non-null 
# can be detected when the nominal local fdr is set at alpha
# alpha = .1 corresponds with proportion 60 pct
# Efdr = .215:
# local fdr for a typical non-null feature is expected to be 21.5%


# delta = -0.605: null distribution is shifted to negative values
# sigma = 1.621: standard deviation becomes larger than 1
# blue and green line differ on the left side = signal genes

fdr_x$Efdr
lfdr_level <- fdr_x$Efdr[1]
t_stats <- t_stats %>%
  mutate(
    lfdr = fdr_x$fdr,
    zfdr = (lfdr < lfdr_level) * z_values)

t_stats %>%
  filter(lfdr < lfdr_level) %>%
  summarize(nr_of_genes=n(), mean_fdr=mean(lfdr))


