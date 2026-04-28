# Authors: 
# Jessica G Young, HMS, Boston
# Soren Harnois-Leblanc, HMS, Boston

# Version: 2026/04/28


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


##############################
### Parametric g-computation
##############################

### 1. Fit a parametric model for E[Y |C = 0, A, L ], as much as flexible as possible, among c=0, save as object

mod_homa = lm(Y_homa ~ A*L, data=dat)
summary(mod_homa)
nobs(mod_homa) 

### 2. Make a copy of original dataset with 10000 rows for each intervention of interest. Drop values of Y, and replace the values of A_i with the value it would take under intervention (A^g+_i)

###Copy 1: the pragmatic intervention we've conceived.
newdatg = dplyr::select(dat,-Y_homa)
newdatg = newdatg %>% mutate(A=if_else(A>0, A-1, A))
#check if g+ worked:
table(newdatg$A)
table(dat$A)

###Copy 2: the intervention set everyone to exposure level 0.
newdat0 = dplyr::select(dat,-Y_homa)
newdat0 = newdat0 %>% mutate(A=0)
table(newdat0$A)

###Copy 3: the intervention set everyone to exposure level 1.
newdat1 = dplyr::select(dat,-Y_homa)
newdat1 = newdat1 %>% mutate(A=1)
table(newdat1$A)

###Copy 4: the intervention set everyone to exposure level 2.
newdat2 = dplyr::select(dat,-Y_homa)
newdat2 = newdat2 %>% mutate(A=2)
table(newdat2$A)

###Copy 5: the natural course (no intervention), don't change exposure but have a data set that will predict Y for all of the participants.
newdatnc = dplyr::select(dat,-Y_homa)

### 3.predict outcome values with model above, conditional on remaining uncensored, having person i's level of L, but now with exposure set to the value under intervention 

homa_g = predict.lm(mod_homa, newdata=newdatg)
homa_0 = predict.lm(mod_homa, newdata=newdat0)
homa_1 = predict.lm(mod_homa, newdata=newdat1)
homa_2 = predict.lm(mod_homa, newdata=newdat2)
homa_nc = predict.lm(mod_homa, newdata=newdatnc)

### 4. Values are averaged over n, baseline 10,000 sample

mean_homa_g = sum(homa_g)/n 
mean_homa_0 = sum(homa_0)/n 
mean_homa_1 = sum(homa_1)/n 
mean_homa_2 = sum(homa_2)/n 
mean_homa_nc = sum(homa_nc)/n

#################
# Values of interest:
#################

# Pragmatic intervention Y_g - NC
mean_homa_g-mean_homa_nc 

# Strict intervention Y_a=0 - Y_a=2
mean_homa_0-mean_homa_2 

# Strict intervention (so Y_a=0 - Y_a=1)
mean_homa_0-mean_homa_1 



#######################################
###  Estimator inverse probability weighting
#######################################

# here we will need C 
table(dat$C)

#### 1. in the original dataset, create 2 indicator variables taking values 1,0
#A0 = 1 when A = 0 and 0 otherwise, 
#A1 = A1=1 when A=1 and 0 otherwise (people with A=2 will have 0 for both).

table(dat$A)

dat = dat %>% mutate(A0 = if_else(A==0, 1, 0),
                     A1 = if_else(A==1, 1, 0))

#check 
table(dat$A)
table(dat$A0)
table(dat$A1)


#### 2. Fit a model with A0 ~ L, such that predictions from this model give us an estimate of p_0 = P[A=0|L].

mod_A0 = glm(A0 ~ L, data=dat, family="binomial") 
summary(mod_A0)
# predicted probability of A0 = p_0
p0 = predict(mod_A0, type="response")
summary(p0) 

#### 3. Fit a model with A1 ~ L, but restricted to records who have A=1 or A=2 (excluding everyone with A=0). 
# The predictions from this model are an estimate of p_1 = P[A=1|A\neq0,L]

mod_A1 = glm(A1 ~ L,  data=dat[!dat$A==0,], family="binomial") 
summary(mod_A1)
nobs(mod_A1) 

# predicted probability of A1 = p_1
p1 = predict(mod_A1, type="response", newdata =dat)
summary(p1)

#### 4. take 1-p1 = p2 for all n people. 1-p1 is an estimate of p_2 = P[A=2|A\neq0,L]
p2 = 1-p1


#### 5. Calculate the treatment weights tw_i for any level of A under each intervention. 

#For our conceived natural value intervention g
dat = dat %>% mutate(twg=case_when(A=="2" ~ 0,
                                   A=="1" ~ (p2*(1-p0))/(p1*(1-p0)),
                                   A=="0" ~ 1 + ((p1*(1-p0))/p0)))



#For the intervention force everyone's exposure to 0
dat = dat %>% mutate(tw0=case_when(A=="2" ~ 0,
                                   A=="1" ~ 0,
                                   A=="0" ~ 1/p0))


#For the intervention force everyone's exposure to 1
dat = dat %>% mutate(tw1=case_when(A=="2" ~ 0,
                                   A=="1" ~ 1/(p1*(1-p0)),
                                   A=="0" ~ 0))


#For the intervention force everyone's exposure to 2
dat = dat %>% mutate(tw2=case_when(A=="2" ~ 1/(p2*(1-p0)),
                                   A=="1" ~ 0,
                                   A=="0" ~ 0))


#### 6. Calculate the censoring weights cw_i, for outcome homa-ir

mod_cens = glm(C ~  A*L, data=dat, family="binomial")
dat$pcens1 = predict(mod_cens, newdata=dat, type="response")

#calculate the inverse probability of C = 0 (not censored)), and assign 0 for those C=1
dat$cw = ifelse(dat$C==0, 1/(1-dat$pcens1), 0) ## predicted probability of P(C=0|A,L) i.e., not being censored

#distribution in those with Y observed:
summary(dat$cw[!is.na(dat$Y_homa)])
#distribution in those with Y missing:
summary(dat$cw[is.na(dat$Y_homa)])


#### 7. multiply tw_i by cw_i, this is the global weight w_i

dat$wg = dat$twg*dat$cw
dat$w0 = dat$tw0*dat$cw
dat$w1 = dat$tw1*dat$cw
dat$w2 = dat$tw2*dat$cw

summary(dat$wg)
summary(dat$w0)
summary(dat$w1)
summary(dat$w2) #we see that weights for A=2 have higher values

#### 8. multiply value of Y by w_i : this gives the weighted outcome values under intervention g+
# then calculate the average of weighted outcome with n as denominator

#subset data to exclude participants with homa-ir missing 
dat_homa_obs = subset(dat, !is.na(dat$Y_homa))

#product of y by w
weighted_homa_g = dat_homa_obs$Y_homa*dat_homa_obs$wg
weighted_homa_0 = dat_homa_obs$Y_homa*dat_homa_obs$w0
weighted_homa_1 = dat_homa_obs$Y_homa*dat_homa_obs$w1
weighted_homa_2 = dat_homa_obs$Y_homa*dat_homa_obs$w2
weighted_homa_nc = dat_homa_obs$Y_homa*dat_homa_obs$cw

#sum the weighted homa values for the numerator, and divide by n (10,000)
mean_weighted_homa_g = sum(weighted_homa_g)/n 
mean_weighted_homa_g 

mean_weighted_homa_0 = sum(weighted_homa_0)/n 
mean_weighted_homa_0 

mean_weighted_homa_1 = sum(weighted_homa_1)/n 
mean_weighted_homa_1 

mean_weighted_homa_2 = sum(weighted_homa_2)/n 
mean_weighted_homa_2 

mean_weighted_homa_nc = sum(weighted_homa_nc)/n 
mean_weighted_homa_nc


#################
# Values of interest:
#################

# Pragmatic intervention Y_g - NC
Yg_Ync = mean_weighted_homa_g-mean_weighted_homa_nc
Yg_Ync 

# Strict intervention Y_a=0 - Y_a=2
Y2_Y0 = mean_weighted_homa_0-mean_weighted_homa_2 
Y2_Y0 

# Strict intervention Y_a=0 - Y_a=1
Y1_Y0 = mean_weighted_homa_0-mean_weighted_homa_1 
Y1_Y0 

###################################################
# Then use bootstrap for 95%CI - see code on github
####################################################

