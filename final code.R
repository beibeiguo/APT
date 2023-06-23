
rm(list=ls()); 

library(mvtnorm)
library(twang)
library(CBPS)
library(dfcrm)
library(rjags)
library(R2jags)
library(UBCRM)
library(BOIN)
mu <- c(0,0,0,0,0,0)
Sigma <- matrix(0,6,6)
diag(Sigma) <- .25
for (i in 1:6)
  for (j in 1:6) {
    if (i!=j)
      Sigma[i,j] <- .025
  }



outcome <- function(n.pat,Z,S,beta00,beta01,beta1,gamma2,gamma3) {
  x <- rmvnorm(n.pat,mean=mu[1:5],sigma=Sigma[1:5,1:5])
  
  beta0 <- ifelse(Z==0,beta00,beta01)
  logit <- beta0 + x%*%beta1
  prob <- 1/(1+exp(-logit))
  y <- rbinom(n.pat,1,prob)
  list(y=y,x=x)
}




get.tau.h <- function(data.comb,K,ps.model,n.control) {
  data.trt <- subset(data.comb,(Z==1 & U==1))
  x.trt <- data.trt[,2:6]
  y.trt <- data.trt$y
  mean.y.trt <- mean(y.trt)
  n.trt <- nrow(data.trt)
  

  data.concurrent <- subset(data.comb,(Z==0 & U==1))
  x.concurrent <- data.concurrent[,2:6]
  y.concurrent <- data.concurrent$y
  mean.y.concurrent <- mean(y.concurrent)
  n.cc <- length(y.concurrent)
  

  data.nonconcurrent <- subset(data.comb,(Z==0 & U==0))
  S.nonconcurrent.vec <- data.nonconcurrent$S
  y.nonconcurrent.vec <- data.nonconcurrent$y
  x.nonconcurrent.vec <- data.nonconcurrent[,2:6]
  
  S.nonconcurrent <- vector("list",K)
  y.nonconcurrent <- vector("list",K)
  x.nonconcurrent <- vector("list",K)
  mean.y.nonconcurrent <- rep(0,K)

  
  for (k in 1:K) {
    nonconcurrent.ind.k <- which(S.nonconcurrent.vec>=(k-1)*2/K & S.nonconcurrent.vec<2*k/K)
    S.nonconcurrent[[k]] <- S.nonconcurrent.vec[nonconcurrent.ind.k]
    y.nonconcurrent[[k]] <- y.nonconcurrent.vec[nonconcurrent.ind.k]
    x.nonconcurrent[[k]] <- x.nonconcurrent.vec[nonconcurrent.ind.k,]
    mean.y.nonconcurrent[k] <- mean(y.nonconcurrent[[k]])
  }
  
  
  ESS <- rep(0,K)
  
  data.trt <- cbind(y.trt,x.trt,1,1)
  colnames(data.trt) <- c("y","x1","x2","x3","x4","x5","Z","U")
  
  data.concurrent <- cbind(y.concurrent,x.concurrent,0,1)
  colnames(data.concurrent) <- c("y","x1","x2","x3","x4","x5","Z","U")
  
  y.hat <- rep(0,K)
  n.Qk <- rep(0,K)
  y.nonconcurrent.comb <- matrix(0,n.control,K)
  
  for (k in 1:K) {
    data.Qk <- cbind(y.nonconcurrent[[k]],x.nonconcurrent[[k]],0,0)
    colnames(data.Qk) <- c("y","x1","x2","x3","x4","x5","Z","U")
    n.Qk[k] <- nrow(data.Qk)
    
    data.comb.k <- rbind(data.Qk,data.trt,data.concurrent)
    colnames(data.comb.k) <- c("y","x1","x2","x3","x4","x5","Z","U")
    data.comb.k <- data.frame(data.comb.k)
    n.U0 <- nrow(data.Qk)
    
    if (ps.model==1) {
      
      s.x <- predict(glm(U~(x1+x2+x3+x4+x5)^2,data=data.comb.k,family="binomial"),type="response")
      w0.x <- s.x[1:n.U0]/(1-s.x[1:n.U0])
    }
    if (ps.model==2) {
      ps_GBM <- ps(U~x1+x2+x3+x4+x5,data=data.comb.k,estimand="ATT",n.trees=1000,shirinkage=.01,stop.method="es.mean")
      w0.x <- get.weights(ps_GBM,estimand="ATT",stop.method="es.mean")[1:n.U0]
    }
    if (ps.model==3) {
      ps_CBPS <- CBPS(U~x1+x2+x3+x4+x5,data=data.comb.k,ATT=1)
      w0.x <- ps_CBPS$weights[1:n.U0]
    }
    
    ESS[k] <- (sum(w0.x))^2/sum(w0.x^2)
    y.hat[k] <- sum(y.nonconcurrent[[k]]*w0.x)/sum(w0.x)
    
    weighted.y <- ESS[k]*y.nonconcurrent[[k]]*w0.x/sum(w0.x)
    y.nonconcurrent.comb[,k] <- c(weighted.y,rep(0,n.control-n.Qk[k]))
  }
  mean.y.control <- (n.cc*mean.y.concurrent+sum(ESS*y.hat))/(n.cc+sum(ESS))
  tau.hat.c <- mean.y.trt-mean.y.control
  

  
  list(tau.hat.c=tau.hat.c,ESS=ESS,theta0=mean.y.control,y.trt=y.trt,y.concurrent=y.concurrent,y.nonconcurrent=y.nonconcurrent.comb,n.trt=n.trt,n.cc=n.cc,n.Qk=n.Qk)
  
}

model.biomarker <- function() {
  
  for (i in 1:n.trt) {
    y.trt[i] ~ dbin(p.trt,1)
  }
  logit(p.trt) <- beta[K+1]+alpha
  
  for (i in 1:n.cc) {
    y.concurrent[i] ~ dbin(p.cc,1)
  } 
  logit(p.cc) <- beta[K+1]
  
  for (k in 1:K) {
    for (i in 1:n.Qk[k]) {
      
      logit(p.nc[i,k]) <- beta[k]
      prob[i,k] <- (p.nc[i,k]^y.nonconcurrent[i,k])*((1-p.nc[i,k])^(1-y.nonconcurrent[i,k]))
      ones[i,k] ~ dbern(prob[i,k]/C)
    }
  }
  
   
  alpha ~ dnorm(0,1/.65^2) 
  beta[K+1] ~ dnorm(-.97,1/.38^2)
  for (k in 1:K) {
    beta[k] ~ dnorm(beta[k+1],1)
  }
  
  
}

get.post <- function(data.comb,K,n.control,ps.model) {
  ones <- matrix(1,n.control,K)
  get.tau <- get.tau.h(data.comb,K,ps.model,n.control)
  y.trt=get.tau$y.trt
  y.concurrent=get.tau$y.concurrent
  y.nonconcurrent=get.tau$y.nonconcurrent
  n.trt=get.tau$n.trt
  n.cc=get.tau$n.cc
  n.Qk=get.tau$n.Qk
  theta0=get.tau$theta0
  ESS=get.tau$ESS
  dat.list <- list(K=K,y.trt=y.trt,y.concurrent=y.concurrent,y.nonconcurrent=y.nonconcurrent,n.trt=n.trt,n.cc=n.cc,n.Qk=n.Qk,C=10000,ones=ones)
  params <- c("alpha","beta")
  mcmcmodel <- jags(data=dat.list,model.file=model.biomarker,DIC=F,progress.bar="none",
                    param=params,
                    n.chains=1,n.iter=1500,n.burnin=300,n.thin=1)
  
  param.mcmc <- as.mcmc(mcmcmodel)[[1]]
  
  alpha.mcmc <- param.mcmc[,1]
  prob.greater <- mean(alpha.mcmc>0)
  
  
  list(prob.greater=prob.greater,ESS=ESS,theta0=theta0)
}




main <- function(ps.model,K,n.control,n.trt,nsim) {
    beta00 <- -1.5
    beta01 <- -1.5
    beta1 <- rep(0.5,5)
    gamma2 <- 0
    gamma3 <- rep(0,5)
  
  
  true.theta0 <- .202
  
  pval.sim <- rep(0,nsim)
  ESS.sim <- rep(0,nsim)
  theta0.sim <- rep(0,nsim)
  
  for (ite in 1:nsim) {
    
    S.trt <- runif(n.trt,2,3)
    outcome.trt <- outcome(n.trt,Z=1,S.trt,beta00,beta01,beta1,gamma2,gamma3)
    y.trt <- outcome.trt$y
    x.trt <- outcome.trt$x[,1:5]
    data.trt <- cbind(y.trt,x.trt,S.trt,1,1)
    
    S.control <- runif(n.control,0,3)
    outcome.control <- outcome(n.control,Z=0,S.control,beta00,beta01,beta1,gamma2,gamma3)
    y.control <- outcome.control$y
    x.control <- outcome.control$x[,1:5]
    
    concurrent.ind <- which(S.control>2)
    y.concurrent <- y.control[concurrent.ind]
    x.concurrent <- x.control[concurrent.ind,]
    S.concurrent <- S.control[concurrent.ind]
    data.concurrent <- cbind(y.concurrent,x.concurrent,S.concurrent,0,1)
    
    nonconcurrent.ind <- which(S.control<=2)
    y.nonconcurrent <- y.control[nonconcurrent.ind]
    x.nonconcurrent <- x.control[nonconcurrent.ind,]
    S.nonconcurrent <- S.control[nonconcurrent.ind]
    data.nonconcurrent <- cbind(y.nonconcurrent,x.nonconcurrent,S.nonconcurrent,0,0)
    
    data.comb <- rbind(data.trt,data.concurrent,data.nonconcurrent)
    colnames(data.comb) <- c("y","x1","x2","x3","x4","x5","S","Z","U")
    data.comb <- data.frame(data.comb)
    
    
    post <- get.post(data.comb,K,n.control,ps.model)
    pval.sim[ite] <- post$prob.greater
    ESS.sim[ite] <- sum(post$ESS)
    theta0.sim[ite] <- post$theta0
    
    
  }
  bias <- mean(theta0.sim)-true.theta0
  
  typeI <- mean(pval.sim>=.95)
  MSE <- var(theta0.sim)+bias^2
  ESS <- mean(ESS.sim)
  list(typeI=typeI,ESS=ESS,bias=bias,MSE=MSE)
}

set.seed(1)
tmp <- main(ps.model=1,K=2,n.control=150,n.trt=100,nsim=2000) 
write.table(tmp$typeI,"typeI.txt",quote=F,row.names=F) 
write.table(tmp$ESS,"ESS.txt",quote=F,row.names=F) 
write.table(tmp$bias,"bias.txt",quote=F,row.names=F) 
write.table(tmp$MSE,"MSE.txt",quote=F,row.names=F) 
