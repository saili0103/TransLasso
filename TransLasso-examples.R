###An example
##coefficient generating functions used in the simulation
Coef.gen<- function(s, h,q=30, size.A0, M, sig.beta,sig.delta1, sig.delta2, p, exact=T){
  beta0 <- c(rep(sig.beta, s), rep(0, p - s))
  W <- matrix(0, nrow = p, ncol = M)# ten prior estimates
  samp0<- sample(1:p, h, replace=F)
  for(k in 1:M){
    if(k <= size.A0){
      W[,k]<-beta0   ###h^*=0
      if(exact){
        W[samp0,k] <-W[samp0,k] + rep(-sig.delta1,h)
      }else{
        W[1:(p/2),k] <-W[1:(p/2),k] + rnorm(p/2, 0, h/p*2)
      }
      
    }else{
      W[,k]<-0
      samp1 <- sample(1:p, q, replace = F)
      if(exact){
        W[samp1,k] <-W[samp1,k] + rep(-sig.delta2,q)
      }else{
        W[1:(p/2),k] <-W[1:(p/2),k] + rnorm(p/2, 0, sig.delta2)
      }
    }
  }
  return(list(W=W, beta0=beta0))
}


source("~/TransLasso-functions.R")
set.seed(123)
p = 500
s = 16
M = 20
sig.beta = 0.3

sig.z <- 1
n0 <- 150
M = 20
n.vec <- c(n0, rep(100, M))
Sig.X <- diag(1, p)
Niter = 200
l1=T

size.A0 = 12 
h=6
A0 = 1:size.A0
beta0<- 
coef.all <-Coef.gen( s, h = h, q = 2*s, size.A0 = size.A0,  M = M,   sig.beta = sig.beta,
             sig.delta1 = sig.beta, sig.delta2 = sig.beta, p = p, exact=F)
B <- cbind(coef.all$beta0, coef.all$W)
beta0 <- coef.all$beta0

###generate the data
X <- NULL
y <- NULL
for (k in 1:(M + 1)) {
  X <- rbind(X, rmvnorm(n.vec[k], rep(0, p), Sig.X))
  ind.k <- ind.set(n.vec, k)
  y <- c(y, X[ind.k, ] %*% B[, k] + rnorm (n.vec[k], 0, 1))
}
###compute init beta ####
mse.vec<-rep(NA,6)
beta.init <-
  as.numeric(glmnet(X[1:n.vec[1], ], y[1:n.vec[1]], lambda = sqrt(2 * log(p) / n.vec[1]))$beta)
mse.vec[1] = mse.fun(as.numeric(beta.init), beta0)$est.err

######Oracle Trans-Lasso#######
if (size.A0 == 0) {
  beta.kA <- beta.init
} else{
  beta.kA <- las.kA(X, y, A0 = 1:size.A0, n.vec = n.vec, l1=l1)$beta.kA
}
mse.vec[2] = mse.fun(as.numeric(beta.kA), beta0)$est.err
###########Trans-Lasso#############
prop.re1 <- Trans.lasso(X, y, n.vec, I.til = 1:50, l1 = l1)
prop.re2 <- Trans.lasso(X, y, n.vec, I.til = 101:n.vec[1], l1=l1)
if(size.A0 > 0 & size.A0< M){ #Rank.re characterizes the performance of the sparsity index Rk
  Rank.re<- (sum(prop.re1$rank.pi[1:size.A0]<=size.A0) +
                     sum(prop.re2$rank.pi[1:size.A0]<=size.A0))/2/size.A0
}else{ Rank.re <- 1 }
beta.prop <- (prop.re1$beta.hat + prop.re2$beta.hat) / 2
mse.vec[3] = mse.fun(beta.prop, beta0)$est.err

######A method for comparison: it has the same pipeline of the Trans-Lasso 
###but with sparsity index R_k=\|w^{(k)}-\beta\|_1 and a naive aggregation (empirical risk minimization)
prop.sp.re1 <- Trans.lasso.sp(X, y, n.vec, I.til = 1:50, l1 = l1)
prop.sp.re2 <- Trans.lasso.sp(X, y, n.vec, I.til = 101:n.vec[1], l1=l1)
if(size.A0 > 0 & size.A0< M){
  Rank.re.sp <- (sum(prop.sp.re1$rank.pi[1:size.A0]<=size.A0) +
                        sum(prop.sp.re2$rank.pi[1:size.A0]<=size.A0))/2/size.A0
}else{ Rank.re.sp <-1 }
beta.sp <- (prop.sp.re1$beta.sp + prop.sp.re2$beta.sp) / 2
mse.vec[4] = mse.fun(beta.sp, beta0)$est.err    

######another method for comparison: it is the same as Trans-Lasso except 
##that the bias correction step (step 2 of Oracle Trans-Lasso) is omitted
beta.pool<-(prop.re1$beta.pool+prop.re2$beta.pool)/2
mse.vec[5] = mse.fun(beta.pool, beta0)$est.err  
#####naive translasso: simply assumes A0=1:K
beta.all <- las.kA(X, y, A0 = 1:M, n.vec = n.vec, l1=l1)$beta.kA#naive transLasso
mse.vec[6] = mse.fun(as.numeric(beta.all), beta0)$est.err

cat(mse.vec, '\n')
cat(Rank.re, Rank.re.sp, '\n')



