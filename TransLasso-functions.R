
library(glmnet)
agg.fun<- function(B, X.test,y.test, total.step=10, selection=F){
  if(sum(B==0)==ncol(B)*nrow(B)){
    return(rep(0,nrow(B)))
  }
  p<-nrow(B)
  K<-ncol(B)
  colnames(B)<-NULL
  if(selection){#select beta.hat with smallest prediction error
    khat<-which.min(colSums((y.test-X.test%*%B)^2))
    theta.hat<-rep(0, ncol(B))
    theta.hat[khat] <- 1
    beta=B[,khat]
    beta.ew=NULL
  }else{#Q-aggregation
    theta.hat<- exp(-colSums((y.test-X.test%*%B)^2)/2)
    theta.hat=theta.hat/sum(theta.hat)
    theta.old=theta.hat
    beta<-as.numeric(B%*%theta.hat)
    beta.ew<-beta
    # theta.old=theta.hat
    for(ss in 1:total.step){
      theta.hat<- exp(-colSums((y.test-X.test%*%B)^2)/2+colSums((as.vector(X.test%*%beta)-X.test%*%B)^2)/8)
      theta.hat<-theta.hat/sum(theta.hat)
      beta<- as.numeric(B%*%theta.hat*1/4+3/4*beta)
      if(sum(abs(theta.hat-theta.old))<10^(-3)){break}
      theta.old=theta.hat
    }
  }
  list(theta=theta.hat, beta=beta, beta.ew=beta.ew)
}


###oracle Trans-Lasso
las.kA<-function(X, y, A0, n.vec, lam.const=NULL, l1=T){
  p<-ncol(X)
  size.A0<- length(A0)
  if(size.A0 > 0){
    ind.kA<- ind.set(n.vec, c(1, A0+1))
    ind.1<-1:n.vec[1]
    if(l1){
      y.A<-y[ind.kA]
    }else{ #the l0-method
      y.A<- y[ind.1]
      Sig.hat<-t(X)%*%X/nrow(X)
      for(k in 1:size.A0){
        ind.k<- ind.set(n.vec,k+1)
        lam.k <- sqrt(mean(y[ind.1]^2)/n.vec[1]+mean(y[ind.k]^2)/n.vec[k]) * sqrt(2*log(p))
        delta.hat.k<-lassoshooting(XtX=Sig.hat, 
                                   Xty=t(X[ind.k,])%*%y[ind.k]/n.vec[k+1]-t(X[1:n.vec[1],])%*%y[1:n.vec[1]]/n.vec[1],
                                   lambda=lam.k)$coef
        y.A<-c(y.A, y[ind.k]-X[ind.k,]%*%delta.hat.k)
      }
    }
    if(is.null(lam.const)){
      cv.init<-cv.glmnet(X[ind.kA,], y.A, nfolds=8, lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/length(ind.kA)))
      lam.const <- cv.init$lambda.min/sqrt(2*log(p)/length(ind.kA))
    }
    w.kA <- as.numeric(glmnet(X[ind.kA,], y.A, lambda=lam.const*sqrt(2*log(p)/length(ind.kA)))$beta)
    w.kA<-w.kA*(abs(w.kA)>=lam.const*sqrt(2*log(p)/length(ind.kA)))
    # cv.delta<-cv.glmnet(x=X[ind.1,],y=y[ind.1]-X[ind.1,]%*%w.kA, lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/length(ind.1)))
    #delta.kA<-predict(cv.delta, s='lambda.min', type='coefficients')[-1]
    delta.kA <- as.numeric(glmnet(x=X[ind.1,],y=y[ind.1]-X[ind.1,]%*%w.kA, lambda=lam.const*sqrt(2*log(p)/length(ind.1)))$beta)
    delta.kA<-delta.kA*(abs(delta.kA)>=lam.const*sqrt(2*log(p)/length(ind.1)))
    beta.kA <- w.kA + delta.kA
    lam.const=NA
  }else{
    cv.init<-cv.glmnet(X[1:n.vec[1],], y[1:n.vec[1]], nfolds=8, lambda=seq(1,0.1,length.out=20)*sqrt(2*log(p)/n.vec[1]))
    lam.const<-cv.init$lambda.min/sqrt(2*log(p)/n.vec[1])
    beta.kA <- predict(cv.init, s='lambda.min', type='coefficients')[-1]
    w.kA<-NA
  }
  list(beta.kA=as.numeric(beta.kA),w.kA=w.kA, lam.const=lam.const)
  
}

#Trans Lasso method
Trans.lasso <- function(X, y, n.vec, I.til, l1=T){
  M= length(n.vec)-1
  #step 1
  X0.til<-X[I.til,] #used for aggregation
  y0.til<-y[I.til]
  X<- X[-I.til,]
  y<-y[-I.til]
  #step 2
  Rhat <- rep(0, M+1)
  p<- ncol(X)
  n.vec[1]<- n.vec[1]-length(I.til)
  ind.1<-ind.set(n.vec,1)
  for(k in 2: (M+1)){
    ind.k<-ind.set(n.vec,k)
    Xty.k <- t(X[ind.k,])%*%y[ind.k]/n.vec[k] - t(X[ind.1,])%*%y[ind.1]/ n.vec[1]
    margin.T<-sort(abs(Xty.k),decreasing=T)[1:round(n.vec[1]/3)]
    Rhat[k] <-  sum(margin.T^2)
  }
  Tset<- list()
  k0=0
  kk.list<-unique(rank(Rhat[-1]))
  #cat(rank(Rhat[-1]),'\n')
  for(kk in 1:length(kk.list)){#use Rhat as the selection rule
    Tset[[k0+kk]]<- which(rank(Rhat[-1]) <= kk.list[kk])
  }
  k0=length(Tset)
  Tset<- unique(Tset)
  #cat(length(Tset),'\n')
  
  beta.T<-list()
  init.re<-las.kA(X=X, y=y, A0=NULL, n.vec=n.vec, l1=l1)
  beta.T[[1]] <- init.re$beta.kA
  beta.pool.T<-beta.T ##another method for comparison
  for(kk in 1:length(Tset)){#use pi.hat as selection rule
    T.k <- Tset[[kk]]
    re.k<- las.kA(X=X, y=y, A0=T.k, n.vec=n.vec, l1=l1, lam.const=init.re$lam.const)
    beta.T[[kk+1]] <-re.k$beta.kA
    beta.pool.T[[kk+1]]<-re.k$w.kA
  }
  beta.T<-beta.T[!duplicated((beta.T))]
  beta.T<- as.matrix(as.data.frame(beta.T))
  agg.re1 <- agg.fun(B=beta.T, X.test=X0.til, y.test=y0.til)
  beta.pool.T<-beta.pool.T[!duplicated((beta.pool.T))]
  beta.pool.T<- as.matrix(as.data.frame(beta.pool.T))
  agg.re2<-agg.fun(B=beta.pool.T, X.test=X0.til, y.test=y0.til)
  
  return(list(beta.hat=agg.re1$beta, theta.hat=agg.re1$theta, rank.pi=rank(Rhat[-1]),
              beta.pool=agg.re2$beta, theta.pool=agg.re2$theta))
}


#A method for comparison: Trans-Lasso(l1). It has the same pipeline of the Trans-Lasso 
###but with sparsity index R_k=\|w^{(k)}-\beta\|_1 and a naive aggregation (empirical risk minimization)
Trans.lasso.sp <- function(X, y, n.vec, I.til, l1=T){
  M= length(n.vec)-1
  #step 1
  X0.til<-X[I.til,] #used for aggregation
  y0.til<-y[I.til]
  X<- X[-I.til,]
  y<-y[-I.til]
  #step 2
  Rhat <- rep(0, M+1)
  p<- ncol(X)
  n.vec[1]<- n.vec[1]-length(I.til)
  ind.1<-ind.set(n.vec,1)
  init.re<-las.kA(X=X, y=y, A0=NULL, n.vec=n.vec, l1=l1)
  for(k in 2: (M+1)){
    ind.k <- ind.set(n.vec,k)
    w.init.k<- as.numeric(glmnet(X[ind.k,], y[ind.k], lambda=init.re$lam.const*sqrt(2*log(p)/length(ind.k)))$beta)
    Rhat[k] <-  sum(abs(w.init.k-init.re$beta.kA)) ##\|w^{(k)}-\beta\|_1
  }
  Tset<- list()
  k0=0
  kk.list<-unique(rank(Rhat[-1]))
  #cat(rank(Rhat[-1]),'\n')
  for(kk in 1:length(kk.list)){#use pi.hat as selection rule
    Tset[[k0+kk]]<- which(rank(Rhat[-1]) <= kk.list[kk])
  }
  k0=length(Tset)
  Tset<- unique(Tset)
  # cat(length(Tset),'\n')
  
  beta.T<-list()
  beta.T[[1]] <- init.re$beta.kA
  for(kk in 1:length(Tset)){#use pi.hat as selection rule
    T.k <- Tset[[kk]]
    beta.T[[kk+1]] <-las.kA(X=X, y=y, A0=T.k, lam.const=init.re$lam.const,n.vec=n.vec, l1=l1)$beta.kA
  }
  beta.T<-beta.T[!duplicated((beta.T))]
  beta.T<- as.matrix(as.data.frame(beta.T))
  agg.re <- agg.fun(B=beta.T, X.test=X0.til, y.test=y0.til, selection =T)
  return(list(beta.sp=agg.re$beta, theta.sp=agg.re$theta, rank.pi=rank(Rhat[-1])))
}

#computing the MSE
mse.fun<- function(beta,est, X.test=NULL){
  pred.err<-NA
  est.err<- sum((beta-est)^2)
  
  if(!is.null(X.test)){
    pred.err<-  mean((X.test%*%(beta-est))^2)
  }
  return(list(est.err=est.err, pred.err= pred.err))
}

ind.set<- function(n.vec, k.vec){
  ind.re <- NULL
  for(k in k.vec){
    if(k==1){
      ind.re<-c(ind.re,1: n.vec[1])
    }else{
      ind.re<- c(ind.re, (sum(n.vec[1:(k-1)])+1): sum(n.vec[1:k]))
    }
  }
  ind.re
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

