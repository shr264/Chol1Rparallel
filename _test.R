library(dplyr)
library(knitr)  
library(parallel)

compare = function(p,n,z) {
  require(Matrix)
  require(MASS)
  require(spam)
  require(fillgraph)
  require(concentrationgraph)
  s = 0.5  # fraction of negative coefficients
  a = 0.3  # minimum magnitude of non-zero coefficients
  b = 0.7  # maximum magnitude of non-zero coefficients
  
  plower = p*(p-1)/2
  
  set.seed(12345) #seed for generating L
  ## diagonals
  D = runif(p,2,5)
  
  ## off-diagonals
  T = diag(p)
  T[upper.tri(T)] = 0
  T[lower.tri(T)] = (ifelse(runif(plower)<s, -1, 1) * 
  		   ifelse(runif(plower)<z,  1, 0) * 
  		   runif(plower, a, b))
  
  L = diag(1.0/sqrt(D)) %*% T   # cholesky factor
  omega = t(L) %*% L            # omega
  sigma = solve(omega)
  sigma[abs(sigma)<1e-10] = 0   # set numerical error to zero
  
  density = sum(abs(omega[lower.tri(omega)])>0)/choose(p,2)
  
  set.seed(23456 + 1) #seed for generating data
  X = mvrnorm(n, mu=rep(0, p), Sigma=sigma) # observations
  
  X = scale(X, center = TRUE, scale = FALSE) # centered obs
  
  algo1time = proc.time()[3]
  algo1 = concentrationgraph(Y = X, G = abs(omega)>0)$Shat
  algo1time = proc.time()[3] - algo1time
  algo1norm = norm(algo1-omega,type="F")/norm(omega,type="F")
  
  ipftime = proc.time()[3]
  ipfshat = ipf(X,abs(omega)>0,10^(-5))$Shat
  ipftime = proc.time()[3] - ipftime
  ipfnorm = norm(ipfshat-omega,type="F")/norm(omega,type="F")
  
  ipfalgotime = proc.time()[3]
  ipfalgoshat = ipf(X,abs(omega)>0,10^(-5), algo1)$Shat
  ipfalgotime = proc.time()[3] - ipfalgotime
  ipfalgonorm = norm(ipfalgoshat-omega,type="F")/norm(omega,type="F")

return(list(density = density, algo1time = algo1time,algo1norm = algo1norm, ipftime = ipftime,ipfnorm = ipfnorm, ipfalgotime = ipfalgotime, ipfalgonorm = ipfalgonorm, n = n, p = p))
}

p = 1500
cat('p =', p, '...')
#OM
#no_cores <- detectCores() - 1
slurm_cores <- strtoi(Sys.getenv('SLURM_CPUS_PER_TASK'))
print(paste("Number of SLURM cores: ", slurm_cores, " ", typeof(slurm_cores)))
no_cores <- slurm_cores - 1
#OM
cl <- makeCluster(no_cores)
z = c(0.003, 0.001, 0.0005, 0.0001, 0.00001)
n = c(4*p, 2*p, p, floor(0.75*p))
zvec = rep(z,4)
nvec = rep(n,each=5)
clusterExport(cl, varlist = c("p","compare"))
data = clusterMap(cl,function(n,z) compare(p,n,z), nvec, zvec)
stopCluster(cl)
datamat = matrix(unlist(data),nrow = 20, byrow = TRUE)
colnames(datamat) = c("density", "algo1time" ,"algo1norm", "ipftime" ,"ipfnorm", "ipfalgotime", "ipfalgonorm", "n" , "p")
save(datamat,file=paste("data",toString(p),".RData",sep=""))

p = 1000
cat('p =', p, '...')
cl <- makeCluster(no_cores)
z = c(0.006, 0.003, 0.001, 0.0005, 0.0001)
n = c(4*p, 2*p, p, floor(0.75*p))
zvec = rep(z,4)
nvec = rep(n,each=5)
clusterExport(cl, varlist = c("p","compare"))
data = clusterMap(cl,function(n,z) compare(p,n,z), nvec, zvec)
stopCluster(cl)
datamat = matrix(unlist(data),nrow = 20, byrow = TRUE)
colnames(datamat) = c("density", "algo1time" ,"algo1norm", "ipftime" ,"ipfnorm", "ipfalgotime", "ipfalgonorm", "n" , "p")
save(datamat,file=paste("data",toString(p),".RData",sep=""))

p = 500
cat('p =', p, '...')
cl <- makeCluster(no_cores)
z = c(0.01, 0.008, 0.004, 0.001, 0.0005)
n = c(4*p, 2*p, p, floor(0.75*p))
zvec = rep(z,4)
nvec = rep(n,each=5)
clusterExport(cl, varlist = c("p","compare"))
data = clusterMap(cl,function(n,z) compare(p,n,z), nvec, zvec)
stopCluster(cl)
datamat = matrix(unlist(data),nrow = 20, byrow = TRUE)
colnames(datamat) = c("density", "algo1time" ,"algo1norm", "ipftime" ,"ipfnorm", "ipfalgotime", "ipfalgonorm", "n" , "p")
save(datamat,file=paste("data",toString(p),".RData",sep=""))

p = 200
cat('p =', p, '...')
cl <- makeCluster(no_cores)
z = c(0.04, 0.02, 0.01, 0.005, 0.001)
n = c(4*p, 2*p, p, floor(0.75*p))
zvec = rep(z,4)
nvec = rep(n,each=5)
clusterExport(cl, varlist = c("p","compare"))
data = clusterMap(cl,function(n,z) compare(p,n,z), nvec, zvec)
stopCluster(cl)
datamat = matrix(unlist(data),nrow = 20, byrow = TRUE)
colnames(datamat) = c("density", "algo1time" ,"algo1norm", "ipftime" ,"ipfnorm", "ipfalgotime", "ipfalgonorm", "n" , "p")
save(datamat,file=paste("data",toString(p),".RData",sep=""))

p = 100
cat('p =', p, '...')
cl <- makeCluster(no_cores)
z = c(0.07, 0.05, 0.025, 0.01, 0.001)
n = c(4*p, 2*p, p, floor(0.75*p))
zvec = rep(z,4)
nvec = rep(n,each=5)
clusterExport(cl, varlist = c("p","compare"))
data = clusterMap(cl,function(n,z) compare(p,n,z), nvec, zvec)
stopCluster(cl)
datamat = matrix(unlist(data),nrow = 20, byrow = TRUE)
colnames(datamat) = c("density", "algo1time" ,"algo1norm", "ipftime" ,"ipfnorm", "ipfalgotime", "ipfalgonorm", "n" , "p")
save(datamat,file=paste("data",toString(p),".RData",sep=""))

p = 50
cat('p =', p, '...')
cl <- makeCluster(no_cores)
z = c(0.1, 0.05, 0.025, 0.01, 0.001)
n = c(4*p, 2*p, p, floor(0.75*p))
zvec = rep(z,4)
nvec = rep(n,each=5)
clusterExport(cl, varlist = c("p","compare"))
data = clusterMap(cl,function(n,z) compare(p,n,z), nvec, zvec)
stopCluster(cl)
datamat = matrix(unlist(data),nrow = 20, byrow = TRUE)
colnames(datamat) = c("density", "algo1time" ,"algo1norm", "ipftime" ,"ipfnorm", "ipfalgotime", "ipfalgonorm", "n" , "p")
save(datamat,file=paste("data",toString(p),".RData",sep=""))

p = 20
cl <- makeCluster(no_cores)
z = c(0.1, 0.05, 0.025, 0.01, 0.0001)
n = c(4*p, 2*p, p, floor(0.75*p))
zvec = rep(z,4)
nvec = rep(n,each=5)
clusterExport(cl, varlist = c("p","compare"))
data = clusterMap(cl,function(n,z) compare(p,n,z), nvec, zvec)
stopCluster(cl)
datamat = matrix(unlist(data),nrow = 20, byrow = TRUE)
colnames(datamat) = c("density", "algo1time" ,"algo1norm", "ipftime" ,"ipfnorm", "ipfalgotime", "ipfalgonorm", "n" , "p")
save(datamat,file=paste("data",toString(p),".RData",sep=""))
