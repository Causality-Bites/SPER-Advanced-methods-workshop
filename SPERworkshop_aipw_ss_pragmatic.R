# Authors: Lan Wen, UWaterloo, Canada
# Version: 2025/12/01


  library(MASS);library(ResourceSelection);
  library(dplyr); library(glm2); library(SuperLearner)
  library(data.table); library(nnet);
  
  df = read.csv("/Users/sorenleblanc/Mon disque (soren.harnoisleblanc@gmail.com)/Postdoc Harvard/Workshop SPER 2026/ML_data.csv", header=T)
  ### this dataset contains (L1,L2,A,C,Y), where L1 is binary, L2~Unif[0,1]
  ### A has three levels (0,1,2,3), Y is binary outcome
  ### n=10000
  n=dim(df)[1]
  
  ## Here, we will use iterated logistic regression models to model multinomial exposure levels A=0,1,2,3. 
  ## We redefine exposure levels to use in logistic regression models:
  df$A0 = ifelse(df$A==0,1,0); ## A0 = I(A=0)
  df$A1 = ifelse(df$A==1,1,0); ## A1 = I(A=1)
  df$A2 = ifelse(df$A==2,1,0); ## A2 = I(A=2)
  
  ## define superlearner library
  sl.lib = c("SL.gam","SL.glm", "SL.glm.interaction","SL.nnet") ## add more candidates if desired
  
  ## set seed and split data into two folds 
  set.seed(123456)
  train_ind = sample(unique(df$id),size = n/2)
  
  ### Following is estimation procedure for E(Y^g) and theoretical standard error of the estimate.
  ### Idea: sample splitting and cross-fitting — fit nuisance models on one fold (training), predict on another (testing)
  ### and then alternate folds to obtain approximately unbiased estimates.
  estimate = function(df1,df2) ## df1 is training data, df2 is test data
  {
    # Fit a model for exposure A in training sample
    ## P(A=0|L) first
    l = as.data.frame(cbind("L1" = df1$L1, "L2" = df1$L2)) ## define predictors for P(A=0|L) in training data
    y = df1$A0 ## define outcome for P(A=0|L) in training data
    fitA0 = SuperLearner(Y=y, X=cbind(l), SL.library=sl.lib, family=binomial) ## fit model on training data for P(A=0|L)
    dtmp = as.data.frame(cbind("L1" = df2$L1, "L2" = df2$L2)) ## extract predictor data to predict P(A=0|L) on test set
    df2$pred_a0 = predict.SuperLearner(fitA0, newdata = dtmp, type="response", onlySL = F)$pred; ## predict P(A=0|L) in testing data
    ## P(A=1|A>0,L) next
    tmpdat = df1[df1$A>0,] ## calling rows in which A>0 in training data
    l = as.data.frame(cbind("L1" = tmpdat$L1, "L2" = tmpdat$L2)) ## define predictors for P(A=1|A>0,L) in training data
    y = tmpdat$A1  ## define outcome for P(A=1|A>0,L) in training data
    fitA1 = SuperLearner(Y=y, X=cbind(l), SL.library=sl.lib, family=binomial) ## fit model on training data for P(A=1|A>0,L)
    dtmp = as.data.frame(cbind("L1" = df2$L1, "L2" = df2$L2)) ## extract predictor data to predict P(A=1|A>0,L) on testing data
    df2$Apred1 = predict.SuperLearner(fitA1, newdata = dtmp, type="response", onlySL = F)$pred; ## predict P(A=1|A>0,L) in testing data
    df2$pred_a1 = df2$Apred1*(1-df2$pred_a0); ## predict P(A=1|L) = P(A=1|A>0,L)*P(A>0|L) in testing data
    df2$pred_a2 = 1-df2$pred_a0-df2$pred_a1 ## predict P(A=2|L) = 1-P(A=0|L)-P(A=1|L) in testing data
    
    # Fit a model for C: P(C|A,L)
    la = as.data.frame(cbind("L1" = df1$L1, "L2" = df1$L2, "A" = df1$A))  ## define predictors for censoring model in training data
    y = df1$C # define outcome for censoring model in training data
    cfit = SuperLearner(Y=y, X=cbind(la), SL.library=sl.lib, family=binomial) ## fit censoring model on training data
    dtmp = as.data.frame(cbind("L1" = df2$L1, "L2" = df2$L2, "A" = df2$A)) ## extract predictors to predict probability of censoring on testing data
    df2$pred_c = predict.SuperLearner(cfit, newdata = dtmp, type="response", onlySL = F)$pred; ## predicted probability of censoring in testing data
    df2$pred_c = ifelse(df2$C==0, 1-df2$pred_c, df2$pred_c) ## predicted probability of P(C=0|A,L) i.e., not being censored in testing data
    
    ## predict p^g(A|L)/p(A|L)]*[(1-C)/p(C=0|A,L)] in testing data
    # I(A==0)*(1+p[A=1|L]/p[A=0|L])
    num_0 = I(df2$A==0)*(1+df2$pred_a1/df2$pred_a0)
    # I(A==1)*(p[A=2|L]/p[A=1|L])
    num_1 = I(df2$A==1)*(df2$pred_a2/df2$pred_a1)
    ## putting it together
    df2$num_tot = (num_0 + num_1)*(I(df2$C==0)/df2$pred_c) #[p^g(A|L)/p(A|L)]*[(1-C)/p(C=0|A,L)] in testing data
    
    ## fit outcome model in those whose C=0 in training data: E(Y|C=0,A,L)
    tmpdat = df1[df1$C==0,] ## only take data from those whose C=0 in training data
    la = as.data.frame(cbind("L1" = tmpdat$L1, "L2" = tmpdat$L2, "A" = tmpdat$A)) # define predictors in model for E(Y|C=0,A,L) in training data
    y = tmpdat$Y # define outcome for E(Y|C=0,A,L) in training data
    yfit = SuperLearner(Y=y, X=cbind(la), SL.library=sl.lib, family=binomial) # fit outcome model in training data
    ## predict E(Y|A,L,C=0) in testing data
    dtmp = as.data.frame(cbind("L1" = df2$L1, "L2" = df2$L2, "A" = df2$A))  ## extract predictors for prediction of E(Y|A,L,C=0) in testing data
    df2$ypred= predict.SuperLearner(yfit, newdata = dtmp, type="response", onlySL = F)$pred; ## predicted E(Y|C=0,A,L)
    ## predict E(Y|A^g+,L,C=0) in testing data
    ydat = df2; ydat$A=ifelse(ydat$A>0, ydat$A-1, ydat$A); ## set A values according to intervention in testing data for A^g+
    dtmpint = as.data.frame(cbind("L1" = ydat$L1, "L2" = ydat$L2, "A" = ydat$A)) ## extract predictors for prediction of E(Y|A^g+,L,C=0) in testing data 
    df2$ypredintervene = predict.SuperLearner(yfit, newdata = dtmpint, type="response", onlySL = F)$pred; ## predicted E(Y|A^g+,L,C=0) 
    
    ## calculate aipw in testing
    df2$Y = ifelse(df2$C==1, 99999, df2$Y) ## sub in any number so that we can calculate aipw (otherwise will spit out NA)
    aipw_g = mean(df2$num_tot*(df2$Y-df2$ypred) + df2$ypredintervene) ## this solves for \psi^g in equation (13)
    #########################
    #########################
    #########################
    aipw_natural = mean((I(df2$C==0)/df2$pred_c)*(df2$Y-df2$ypred) + df2$ypred) ## This is AIPW estimator for natural course E(Y) (or E(Y^{c=0}) depending if you eliminate censoring)
    diff = aipw_g - aipw_natural ## difference is the E(Y^g) - E(Y)
    
    ## calculate SE from:
    ## \sum_{i=1}^n { [p^g(Ai|Li)/p(Ai|Li)]*[(1-Ci)/p(Ci=0|Ai,Li)](Yi-E(Yi|Ai,Li,Ci=0)) + E(Yi|Ai^g+,Li,Ci=0) - [(1-Ci)/p(Ci=0|Ai,Li)(Yi-E(Yi|Ai,Li,Ci=0)) + E(Yi|Ai^g+,Li,Ci=0)] - \psi_diff^g }^2
    ## i.e., the sum of the squared values in equation (13) over all individuals i.
    individual = df2$num_tot*(df2$Y-df2$ypred) + df2$ypredintervene - (I(df2$C==0)/df2$pred_c*(df2$Y-df2$ypred) + df2$ypred) - diff
    expected1 = 0;
    for(i in 1:nrow(df2))
    {
      expected1 = expected1 + ((individual[i])^2)
    }
    myvar1 = expected1/(n/2)
    
    return(c(diff,myvar1)) ## returns (estimate, standard error)
  }
  
  ## predict aipw_a
  df1  = df[df$id %in% train_ind,]; df2 = df[!df$id %in% train_ind,]; # df1 is training, df2 is testing
  aipw_a =  estimate(df1,df2)
  df1  = df[!df$id %in% train_ind,]; df2 = df[df$id %in% train_ind,]; # swap training and testing
  aipw_b =  estimate(df1,df2)
  aipw = (aipw_a[1] + aipw_b[1])/2
  se = sqrt(((aipw_a[2] + aipw_b[2])/2)/n)
  aipw # E(Y^g)-E(Y) estimate
  se # corresponding estimated standard error
  
  
