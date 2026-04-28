# Authors: 
# Lan Wen, UWaterloo, Canada

# Version: 2026/04/22

library(dplyr)
library(tidyr)
library(MASS)
library(glm2)
library(data.table)


dat = read.csv("/Users/sorenleblanc/Mon disque (soren.harnoisleblanc@gmail.com)/Postdoc Harvard/Workshop SPER 2026/parametric_data.csv", header=T)
str(dat)
n<-dim(dat)[1] #baseline sample size is 10,000

#L is a covariable with values 1 or 0, e.g. sex
table(dat$L)
#0    1 
#4959 5041 

# A is a 3 level exposure with 39% / 45% / 16% distribution
table(dat$A) 
dat$A = as.numeric(dat$A)
#0    1    2 
#3921 4501 1578 
# we will treat it as "fast-food" level 0:<1 time/wk, 1: 1-2 times/wk, 2: >=3times/wk
# this simulated sample tends to consume more than in Project Viva

# C is whether censored yes (1) or no (0)
table(dat$C)
#0    1 
#8306 1694 

# Y is a dichotomous outcome 1/0. 
# Let's use a continuous outcome to represent Homa-ir
# homa-ir is always positive, non-zero, and often right-skewed
set.seed(9999)
dat$Y_homa = rgamma(10000, shape = 3, scale=2)

# set as missing if C=1
dat = dat %>% mutate(Y_homa=if_else(C==1, NA, Y_homa))
hist(dat$Y_homa)
summary(dat$Y_homa)
  
## Here, we will use iterated logistic regression models to model multinomial exposure levels A=0,1,2,. 
## You can also use multinom() in R for parametric models, but iterated logistic regression models are necessary for many existing machine learning algorithms
## We redefine exposure levels to use in logistic regression models:
dat$A0 = ifelse(dat$A==0,1,0); ## A0 = I(A=0)
dat$A1 = ifelse(dat$A==1,1,0); ## A1 = I(A=1)
dat$A2 = ifelse(dat$A==2,1,0); ## A2 = I(A=2)
  
# Fit a model for exposure A
fitA0 = glm(A0 ~ L, data=dat, family=binomial()) ## fit P(A=0|L)
dat$pred_a0 = predict(fitA0, newdata=dat, type="response"); ## predict P(A=0|L)
fitA1 = glm(A1 ~ L, data=dat[!dat$A==0,], family=binomial()) ## fit (A=1|A>0,L) 
dat$Apred1 = predict(fitA1, newdata=dat, type="response"); ## predict P(A=1|A>0,L)
dat$pred_a1 = dat$Apred1*(1-dat$pred_a0); ## predict P(A=1|L) = P(A=1|A>0,L)*P(A>0|L)
dat$pred_a2 = 1-dat$pred_a0-dat$pred_a1 ## predict P(A=2|L) = 1-P(A=0|L)-P(A=1|L)
  
# Fit a model for C
cfit = glm2(C ~ A*L, family=binomial(), data = dat) ; ## censoring model 
dat$pred_c = predict(cfit, newdata = dat, type="response") ## predicted probability of censoring
dat$pred_c = ifelse(dat$C==0, 1-dat$pred_c, dat$pred_c) ## predicted probability of P(C=0|A,L) i.e., not being censored
  
### calculate [I(A=a)/p(A|L)]*[(1-C)/p(C=0|A,L)]
# Treatment weights, static
# I(A==0)/p[A=0|L]
static_0 = I(dat$A==0)/(dat$pred_a0)
# I(A==1)/p[A=1|L]
static_1 =  I(dat$A==1)/(dat$pred_a1) 
static_2 = I(dat$A==2)/(dat$pred_a2) 
## putting the weights (treatment and censoring) together
dat$wt0 = (static_0)*(I(dat$C==0)/dat$pred_c) #[I(A=0)/p(A|L)]*[(1-C)/p(C=0|A,L)]
summary(dat$wt0)
dat$wt1 = (static_1)*(I(dat$C==0)/dat$pred_c) #[I(A=1)/p(A|L)]*[(1-C)/p(C=0|A,L)]
summary(dat$wt1)
dat$wt2 = (static_2)*(I(dat$C==0)/dat$pred_c) #[I(A=2)/p(A|L)]*[(1-C)/p(C=0|A,L)]
summary(dat$wt2)
  
# Treatment weights, pragmatic (g)
# I(A==0)*(1+p[A=1|L]/p[A=0|L])
wtg_0 = I(dat$A==0)*(1+dat$pred_a1/dat$pred_a0)
# I(A==1)*(p[A=2|L]/p[A=1|L])
wtg_1 =  I(dat$A==1)*(dat$pred_a2/dat$pred_a1)
## putting it together
dat$wt_g = (wtg_0 + wtg_1)*(I(dat$C==0)/dat$pred_c) #[p^g(A|L)/p(A|L)]*[(1-C)/p(C=0|A,L)]
summary(dat$wt_g)
  
##################
##################
## fit outcome model in those whose C=0 #remember that Y is continuous here
yfit = lm(Y_homa ~ A*L, data = dat[dat$C==0,])
## predict E(Y|A,L,C=0)
dat$ypred = predict(yfit, newdata = dat, type="response") ## predicted E(Y|A,L,C=0) 
## predict E(Y|A=0,L,C=0) and under A=1, A=2
ydat0 = dat; ydat0$A=0; ## set A values to 0 to predict E(Y|A=0,L,C=0)
ydat1 = dat; ydat1$A=1; ## set A values to 1 to predict E(Y|A=1,L,C=0) 
ydat2 = dat; ydat2$A=2; ## set A values to 2 to predict E(Y|A=2,L,C=0)
dat$ypredintervene0 = predict(yfit, newdata = ydat0, type="response") ## predicted E(Y|A=0,L,C=0) 
dat$ypredintervene1 = predict(yfit, newdata = ydat1, type="response") ## predicted E(Y|A=1,L,C=0) 
dat$ypredintervene2 = predict(yfit, newdata = ydat2, type="response") ## predicted E(Y|A=2,L,C=0) 
  
#under intervention g
ydat = dat; ydat$A=ifelse(ydat$A>0, ydat$A-1, ydat$A); ## set A values according to intervention in testing data for A^g
dat$ypredinterveneg = predict(yfit, newdata = ydat, type="response") ## predicted E(Y|A^g,L,C=0) 
  
  
## calculate aipw
dat$Y_homa = ifelse(dat$C==1, 99999, dat$Y_homa) ## sub in any number so that we can calculate aipw (otherwise will spit out NA)
aipw_0 = mean(dat$wt0*(dat$Y_homa-dat$ypred) + dat$ypredintervene0) ## aipw estimator for E(Y^{a=0})
aipw_1 = mean(dat$wt1*(dat$Y_homa-dat$ypred) + dat$ypredintervene1) ## aipw estimator for E(Y^{a=1})
aipw_2 = mean(dat$wt2*(dat$Y_homa-dat$ypred) + dat$ypredintervene2) ## aipw estimator for E(Y^{a=2})
diff1_0 = aipw_0 - aipw_1 ## difference is the E(Y^{a=0}) - E(Y^{a=1})
diff2_0 = aipw_0 - aipw_2 
  
#under intervention g
aipw_g = mean(dat$wt_g*(dat$Y_homa-dat$ypred) + dat$ypredinterveneg) ## this solves for \psi^g in equation (13)
### in other words, aipw solves for 
### \psi^g = (1/n)*\sum_{i=1}^n [p^g(Ai|Li)/p(Ai|Li)]*[(1-Ci)/p(Ci=0|Ai,Li)(Yi-E(Yi|Ai,Li,Ci=0)) + E(Yi|Ai^g+,Li,Ci=0)
  
#natural course:
aipw_natural = mean((I(dat$C==0)/dat$pred_c)*(dat$Y_homa-dat$ypred) + dat$ypred) ## ## This is AIPW estimator for natural course E(Y) (or E(Y^{c=0}) depending if you eliminate censoring)
diff_g = aipw_g - aipw_natural ## difference is the E(Y^g) - E(Y)
  




