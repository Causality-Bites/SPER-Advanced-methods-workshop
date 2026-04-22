# Authors: 
# Jessica G Young, HMS, Boston
# Soren Harnois-Leblanc, HMS, Boston

# Version: 2026/04/22


#make sure to run the point estimate code to make sure you have the final version of mydat as an object.
mydataorig = dat 
#remove the original mydat
rm(dat)

#use NewID in place of id in the algorithm within each bootstrap sample 
mydataorig$NewID<-as.numeric(as.factor(mydataorig$id))
n<-length(mydataorig$NewID)
K<-10000  #number of replications
matrix.est<-matrix(NA,nrow=K,ncol=8) #8 is the number of estimates you are going after
colnames(matrix.est) = c("Y_a=0", "Y_a=1", "Y_a=2", "Y_a=g", "Y_nc", "StrictA1_A0", "StrictA2_A0", "Pragmatic") #names of the estimates at the end in the output
seed<-456

for (k in 1:K){
  seed<-seed+1
  set.seed(seed)
  bootsamp<-sample(1:n,n,replace=T)
  bootsamp<-data.frame(NewID=bootsamp,bootID=1:length(bootsamp))
  dat<-bootsamp%>% left_join(mydataorig) %>% rename(originalID=NewID) %>% rename(NewID=bootID)
  #length(unique(mydata$NewID))
  
 #HERE COPY AND PASTE YOUR CODE FOR IPW OR G-COMP OR AIPW
  
  outputvec<-rep(NA,each=8)
  outputvec[1]=mean_homa_0
  outputvec[2]=mean_homa_1
  outputvec[3]=mean_homa_2
  outputvec[4]=mean_homa_g
  outputvec[5]=mean_homa_nc
  outputvec[6]=Y1_Y0
  outputvec[7]=Y2_Y0
  outputvec[8]=Yg_Ync
  
  
  matrix.est[k,]<-outputvec
  
}
#calculate variance, CI low, and CI high for each column
VAR=apply(matrix.est,2, FUN=var) 
CIlow=apply(matrix.est,2, FUN=quantile, probs=0.025)
CIhigh=apply(matrix.est,2, FUN=quantile, probs=0.975)
output = rbind(VAR, CIlow, CIhigh)
colnames(output) = colnames(matrix.est)
output  
  
  
  