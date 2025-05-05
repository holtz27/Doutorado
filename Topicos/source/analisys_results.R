### IG
#load("~/topicos/st/sim1/ig/0.1/ig_sa_0.1_1.RData")
#res1=result
#load("~/topicos/st/sim1/ig/0.1/ig_sa_0.1_2.RData")
#res2=result
### Exp
#load("~/topicos/st/sim1/exp/0.1/exp_sa_0.1_1.RData")
#res1=result
#load("~/topicos/st/sim1/exp/0.1/exp_sa_0.1_2.RData")
#res2=result
### PCP
load("~/topicos/st/sim1/pcp/0.1/pcp_sa_0.1_1.RData")
res1=result
load("~/topicos/st/sim1/pcp/0.1/pcp_sa_0.1_2.RData")
res2=result

result = c(res1, res2)
################################################################################
its=length(result)
errors=NULL
for(i in 1:its){
  
  s=result[[i]][[1]]
  if(sum(is.nan(s)+is.infinite(s))!= 0){
    errors = cbind(errors, i) 
  }
  
} 
errors
if(!is.null(errors)){
  Res = result[-errors]
}else{
  Res = result
} 

its=length(Res)
M=result[[1]][[1]]
Reps = matrix(0,nrow=dim(M)[1],ncol=its)
for(i in 1:its){
  Reps[,i] = Res[[i]][[1]][,1]
}

# boxplot
par(mfrow=c(1,5))
boxplot(Reps[1,], main='mu')
boxplot(Reps[2,], main='phi')
boxplot(Reps[3,], main='sh')
boxplot(Reps[4,], main='sa')
boxplot(Reps[5,], main='nu')
par(mfrow=c(1,1))


indx=-1
indx = cbind(indx, which(Reps[1,]>=3))
#indx = cbind(indx, which(Reps[2,] <= 0.8))
#indx = cbind(indx, which(Reps[3,] <= 0.8))
#indx = cbind(indx, which(Reps[4,] <= 0.8))
indx = cbind(indx, which(Reps[5,]>=30))
indx = cbind(indx, which(Reps[5,]<=3))
indx = unique(as.vector(indx))
indx=indx[-1]
indx

Res[[indx[1]]]

#L = result
if(!is.null(indx)) Res=Res[-indx]
its=length(Res)
#M=L[[1]][[1]]
Reps = matrix(0,nrow=dim(M)[1],ncol=its)
for(i in 1:its){
  Reps[,i] = Res[[i]][[1]][,1]
}
# boxplot
par(mfrow=c(1,5))
boxplot(Reps[1,], main='mu')
boxplot(Reps[2,], main='phi')
boxplot(Reps[3,], main='sh')
boxplot(Reps[4,], main='sa')
boxplot(Reps[5,], main='nu')
par(mfrow=c(1,1))

# error percent
its
(length(result)-its)/length(result)*100

#indx2= sample(1:its, 300)
#Reps=Reps[,indx2]
#Res=Res[indx2]
#
theta = theta_vdd
theta[theta == 0] = 1
m=dim(Reps)[1]
n=dim(Reps)[2]
err=min=max=cd=matrix(0,m,n)
prob.cob = 0
for(i in 1:n){
  err[,i]=(Reps[,i]-theta_vdd)/theta
  min[,i] = Res[[i]][[1]][,3]
  max[,i] = Res[[i]][[1]][,5]
  cd[,i] = Res[[i]][[1]][,'CD']
  l1 = as.numeric(Res[[i]][[1]][,3] < theta_vdd)
  l2 = as.numeric(Res[[i]][[1]][,5] > theta_vdd)
  prob.piv = matrix(round(0.5*(l1+l2),0), ncol=1)
  prob.cob = prob.cob + prob.piv
}
Data = cbind(theta_vdd,
             matrix(apply(Reps, MARGIN=1,mean), ncol=1),
             matrix(apply(min, MARGIN=1,mean), ncol=1),
             matrix(apply(max, MARGIN=1,mean), ncol=1),
             matrix(apply(cd, MARGIN=1,mean), ncol=1),
             matrix(apply(err, MARGIN=1,mean), ncol=1),
             matrix(sqrt(apply(err^2, MARGIN=1,mean)), ncol=1),
             prob.cob/length(Res))
row.names(Data) = c('mu', 'phi', 'sh', 'sa', 'v')
colnames(Data) = c('Theta','mean','2.5%','97.5%','CD','bias', 'reqm', 'prob.cob')
round(Data,4)


