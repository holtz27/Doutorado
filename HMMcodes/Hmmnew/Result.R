remove_nulls = function(lst){
  
  null_count <- sum(sapply(lst, is.null))
  cleaned_lst <- Filter(Negate(is.null), lst)
  cat('Elements null in list: ', null_count)
  
  return(cleaned_lst)
}
################################################################################
mu=0.1
phi=0.98
sigma=0.1
beta=c(0.2,0.07,-0.18)
nu=10
theta = matrix(c(beta,mu,phi,sigma,nu), ncol = 1) 

### Simulations results
mean(npd)
mean(times)
reps_cleaned = remove_nulls(reps)
Reps = matrix(0, 
              nrow=dim(reps_cleaned[[1]]$Results)[1], 
              ncol=length(reps_cleaned)
              )
for(i in 1:length(reps_cleaned)){
  Reps[, i] = reps_cleaned[[i]]$Results[, 'mean']
}

# boxplot
par(mfrow=c(2,4))
boxplot(Reps[1,], main='b0')
boxplot(Reps[2,], main='b1')
boxplot(Reps[3,], main='b2')
boxplot(Reps[4,], main='mu')
boxplot(Reps[5,], main='phi')
boxplot(Reps[6,], main='sigma')
boxplot(Reps[7,], main='nu')
par(mfrow=c(1,1))

indx=numeric(0)
d1=dim(Reps)
indx = which(Reps[7,] >= 100)
if(length(indx)!=0) Reps = Reps[, -indx]
d2=dim(Reps)
(d1[2]-d2[2])/d1[2]
###############################################################################
### Metrics evaluation
rel.err = err = matrix(0, 
                       nrow=dim(reps_cleaned[[1]]$Results)[1], 
                       ncol=dim(Reps)[2]
)
for(i in 1:dim(Reps)[2]){
  err[, i] = Reps[, i] - theta
  rel.err[, i] = err[, i]/theta
}

# Mean, MRB, MRAD, RMSE, Prob. cob
m = apply(Reps, 1, mean)
mrb = apply(rel.err, 1, mean)
mrad = apply(abs(rel.err), 1, mean)
rel.mse = apply(rel.err^2, 1, mean)
p.cob = rep(0, 7)
Indxs = 1:length(reps_cleaned)
if(length(indx)!=0) Indxs=Indxs[-indx]
for(i in Indxs){
  M = as.numeric(reps_cleaned[[i]]$Results[, c('2.5%')] < theta)
  M = rbind(M,as.numeric(reps_cleaned[[i]]$Results[, c('97.5%')] > theta))
  p.cob = p.cob + round(apply(M, 2, mean), 0)
}
p.cob = p.cob/length(Indxs)

Result = theta
Result = cbind(Result, m)
Result = cbind(Result, mrb)
Result = cbind(Result, mrad)
#Result = cbind(Result, mse)
Result = cbind(Result, rel.mse)
Result = cbind(Result, p.cob)
Result = round(Result, 4) #signif(Result, 2)
colnames(Result) = c('Theta', 'Mean', 'MRB', 'MRAD', 'Rel.MSE', 'Cover prob.')
row.names(Result) = c('b0', 'b1', 'b2', 'mu', 'phi', 'sigma', 'nu')
Result

