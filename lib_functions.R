library(pre)
library(sandwich)
library(lmtest)
library(MASS) # mvrnorm() function

### Generate the bootstrapped null distribution before running the replications ###

# parameter setting may vary across scenarios
expit <- function(x){
  exp(x)/(1+exp(x))
}

getnullstats <- function(x1.mean, x1.sd, x2.mean, x2.sd, covar, sigsq=NULL, beta0, beta1, beta2, beta3, n, nsamp, f1, family){
  
  # This function generates a vector of nsamp null statistics
  
  # x1.mean, x1.sd, x2.mean, x2.sd: mean and standard deviation of x1 and x2, respectively
  # covar: covariance between x1 and x2
  # sigsq: sigma square, NULL if the outcome is binary
  # n: sample size
  # nsamp: number of bootstrapped samples
  # f1: functional form of x1, can be "linear", "quadratic", "log", or "exp"
  # family: "gaussian" for continuous outcome, "binomial" for binomial outcome
  
  Sigma <- matrix(c(x1.sd^2, covar, covar, x2.sd^2), ncol=2)
  
  x <- mvrnorm(n=n, mu=c(x1.mean, x2.mean), Sigma=Sigma)
  
  if(f1=="log"){
    x1 <- exp(x[,1])
  } else {
    x1 <- x[,1]
  }
  
  x2 <- x[,2]
  
  if(f1=="linear"){ # begin: if 1
    fx1=x1
  } else if(f1=="quadratic"){
    fx1=x1+2*(x1^2)
  } else if(f1=="log"){
    fx1=log(x1)
  } else if(f1=="exp"){
    fx1=exp(x1)
  } # end: if 1
  
  # Generate response variable
  if(family=="gaussian"){ # begin: if 2
    epsilon <- rnorm(n, 0, sqrt(sigsq))
    y <- beta0 + beta1*fx1 + beta2*x2 + beta3*x1*x2 + epsilon
  } else if(family=="binomial"){
    p <- expit(beta0 + beta1*fx1 + beta2*x2 + beta3*x1*x2)
    y <- rbinom(n, 1, p)
  } else {stop("family must be 'gaussian' or 'binomial'!")} # end: if 2
  
  dat <- data.frame(y=y, x1=x1, x2=x2) # data frame to be passed to the function pre()
  if(family=="binomial") {dat$y <- as.factor(dat$y)} # convert the response variable to a factor if binary
  
  # Fit prediction rule ensemble
  rulefit <- do.call("pre", args=list(formula=y~x1+x2, data=dat, tree.unbiased=F, maxdepth=maxdepth_sampler(),
                                      family=family, sampfrac=min(1, (11*sqrt(n)+1)/n), nfolds=round(min(20, max(0, 5200/n-2)))))
  
  
  nullstats <- NULL
  nsamp.tmp <- 10
  
  for(kk in 1:(nsamp/nsamp.tmp)){
    # Compute null-interaction models
    nullmods <- bsnullinteract(rulefit, nsamp=nsamp.tmp, penalty.par.val="lambda.min")
    # Compute null-interaction test statistics
    interact.fit <- interact(rulefit, varnames="x1", nullmods=nullmods, penalty.par.val="lambda.min", plot=F)
    nullstats.tmp <- interact.fit$nullH2$x1
    nullstats <- c(nullstats, nullstats.tmp)
    
  }
  
  return(nullstats) # return a vector of nsamp null statistics
  
}


logit <- function(x) log(x/(1-x))

stdize <- function(x){ # function for standardization
  mu <- mean(x)
  sigma <- sd(x)
  x <- (x-me)/sigma
  return(x)
}

paste0 <- function(..., collapse = NULL) {
  paste(..., sep = "", collapse = collapse)
}

simdata <- function(x1.mean, x1.sd, x2.mean, x2.sd, sigsq=NULL, covar=0, n=n, beta, f1, family="gaussian"){
  # function for simulating the data (x1, x2, and y)
  # f1 is the function form of x1 in the outcome generating model
  # beta is a vector of (beta0, beta1, beta2, beta3)
  # family: "gaussian" for continuous outcome and "binomial" for binary outcome; default is "gaussian"
  
  Sigma <- matrix(c(x1.sd^2, covar, covar, x2.sd^2), ncol=2)
  
  x <- mvrnorm(n=n, mu=c(x1.mean, x2.mean), Sigma=Sigma)
  
  if(f1=="log"){ # begin: if 1
    x1 <- exp(x[,1])
  } else {
    x1 <- x[,1]
  } # end: if 1
  x2 <- x[,2]
  
  ### Function form of x1 ###
  if(f1=="linear"){ # begin: if 2
    fx1 <- x1
  } else if(f1=="quadratic"){
    fx1 <- x1+2*(x1^2)
  } else if(f1=="log"){
    fx1 <- log(x1)
  } else if(f1=="exp"){
    fx1 <- exp(x1)
  } # end: if 2
  
  beta0 <- beta[1]; beta1 <- beta[2]; beta2 <- beta[3]; beta3 <- beta[4]
  
  ### Generate response variable ###
  if(family=="gaussian"){ # begin: if 2
    epsilon <- rnorm(n, 0, sqrt(sigsq))
    y <- beta0 + beta1*fx1 + beta2*x2 + beta3*x1*x2 + epsilon
  } else if(family=="binomial"){
    p <- expit(beta0 + beta1*fx1 + beta2*x2 + beta3*x1*x2)
    y <- rbinom(n, 1, p)
  } else {stop("family must be 'gaussian' or 'binomial'!")} # end: if 2
  
  return(list(x1=x1, x2=x2, fx1=fx1, y=y))
  
}

Wald.test <- function(x1, x2, y, family="gaussian"){
  # family: default is "gaussian"
  
  n <- length(y)
  if(family=="gaussian"){
    mod <- glm(y ~ x1 + x2 + I(x1*x2))
    res <- matrix(residuals(mod), ncol=1) # a vector of estimated residuals (y-y_hat)
    ressq <- matrix(res^2, ncol=1) # a vector of squared residuals (y-y_hat)^2
    RSS <- sum(ressq) # residual sum of squares
    ressq.diag <- diag(as.numeric(ressq), nrow=length(res))
  } else if (family=="binomial"){
    mod <- do.call("glm", args=list(formula=y~x1+x2+I(x1*x2), family=family))
    res <- matrix(residuals(mod, type="response"), ncol=1) # a vector of estimated residuals (y-y_hat)
    ressq <- matrix(res^2, ncol=1) # a vector of squared residuals (y-y_hat)^2
    ressq.diag <- diag(as.numeric(ressq), nrow=length(res))
  }
  
  p <- length(coef(mod)) # number of paramters to be estimated
  x <- matrix(cbind(rep(1,n),x1,x2,x1*x2), ncol=4, byrow=F) # design matrix
  B <- t(x) %*% ressq.diag %*% x
  
  if(family=="gaussian"){
    A <- t(x) %*% x
    var.mb <- (RSS/(n-p))*solve(A) # variance-covariance matrix of betas
    var.san <- (n/(n-p))*(solve(A) %*% B %*% solve(A))
  } else if(family=="binomial"){
    mu.hat <- fitted(mod)
    W <- diag(mu.hat*(1-mu.hat), ncol=length(mu.hat))
    A <- t(x) %*% W %*% x
    var.mb <- solve(A)
    #var.san <- (n/(n-p))*(solve(A) %*% B %*% solve(A))
    var.san <- ((n-p)/n)*(solve(A) %*% B %*% solve(A))
  }
  
  beta.hat <- coef(mod)
  
  Wald.mb <- beta.hat[4]^2/var.mb[4,4]
  Wald.san <- beta.hat[4]^2/var.san[4,4]
  pval.mb <- 1-pchisq(Wald.mb, df=1)
  pval.san <- 1-pchisq(Wald.san, df=1)
  pval.mb.pack <- coef(summary(mod))[4,4] # p-value computed by built-in functions for packages
  pval.san.pack <- coeftest(mod, df=n-p, vcov.=sandwich)['I(x1 * x2)',4]
  return(list(var.sandwich=var.san, var.model=var.mb,
              pval.model=unname(pval.mb), pval.sandwich=unname(pval.san),
              pval.model.pack=unname(pval.mb.pack), pval.sandwich.pack=unname(pval.san.pack)))
  
}

score.test <- function(x1,x2,y,family){
  n <- length(y)
  mod <- do.call(glm, args=list(formula=y ~ x1 + x2 + I(x1*x2), family=family))
  mod.null <- do.call(glm, args=list(formula=y ~ x1 + x2, family=family))
  
  p <- length(coef(mod.null))
  x1x2 <- matrix(x1*x2)
  x <- matrix(cbind(1,x1,x2,x1x2), ncol=4)
  x0 <- matrix(cbind(1,x1,x2), ncol=3)
  
  if(family=="gaussian"){
    
    res <- matrix(residuals(mod.null), ncol=1) # a vector of estimated residuals (y-y_hat)
    ressq <- matrix(res^2, ncol=1) # a vector of squared residuals (y-y_hat)^2
    RSS <- sum(ressq) # residual sum of squares
    W <- diag(as.numeric(ressq), nrow=length(res)) # diagonal matrix of squared residuals
    
    A <- matrix(c(-t(x1x2) %*% x0 %*% solve(t(x0)%*%x0), 1), nrow=1)
    B <- t(x) %*% W %*% x
    
    var.mb <- (RSS/(n-p))*(A %*% t(x) %*% x %*% t(A))
    var.san <- n*(A %*% B %*% t(A))/(n-p)
    
  }
  
  if(family=="binomial"){
    res <- residuals(mod.null, type="response")
    mu.hat <- fitted(mod.null)
    W <- diag(mu.hat*(1-mu.hat), ncol=length(mu.hat))
    U <- diag(res^2, nrow=length(res))
    B <- matrix(c(-t(x1x2) %*% W %*% x0 %*% solve(t(x0) %*% W %*% x0), 1), nrow=1)
    C <- t(x) %*% U %*% x
    var.san <- n*(B %*% C %*% t(B))/(n-p)
    var.mb <- B %*% t(x) %*% W %*% x %*% t(B)
    
  }
  
  
  S <- sum(x1x2*res)
  score.mb <- as.numeric(S^2/var.mb)
  score.san <- as.numeric(S^2/var.san)
  pval.mb <- 1-pchisq(score.mb, df=1)
  pval.san <- 1-pchisq(score.san, df=1)
  pval.mb.pack <- anova(mod,mod.null,test='Rao')[2,6]
  
  return(list(pval.model=unname(pval.mb), pval.sandwich=unname(pval.san),
              pval.model.pack=unname(pval.mb.pack)))
  
}

score.linear = function(x1,x2,y){
  n = length(y)
  mod = glm(y~x1 + x2 + I(x1*x2))
  mod.null = glm(y~x1 + x2)
  x0 = matrix(cbind(rep(1,n),x1,x2),ncol=3)
  x3 = matrix(x1*x2,ncol=1)
  b0_fit = matrix(coef(mod.null),ncol=1)
  res0 = y - matrix(x0%*%b0_fit,ncol=1)
  Vi = sum(res0^2)/(n-3)
  U3 = sum(x3*res0) / Vi
  t1 = matrix(0,nrow=1,ncol=3);t2 = matrix(0,nrow=3,ncol=3)
  D_mb = matrix(0,ncol=4,nrow=4);D_san = matrix(0,ncol=4,nrow=4)
  for(i in 1:n){
    x0i = matrix(x0[i,],nrow=1)
    x3i = x3[i]
    xi = matrix(c(x0i,x3i),nrow=1)
    t1 = t1 + x3i * x0i / Vi
    t2 = t2 + t(x0i) %*% x0i / Vi
    D_mb = D_mb + t(xi) %*% xi / Vi
    D_san = D_san + t(xi) %*% xi * res0[i]^2 / Vi^2
  }
  D_mb = D_mb
  A = matrix(c(-t1 %*% solve(t2),1),nrow=1)
  vs_mb = A %*% D_mb %*% t(A)
  Score_mb = U3^2/vs_mb
  p.score.mb = 1 - pchisq(Score_mb,df = 1) # model-based score test
  #p.score.pack.mb = anova(mod,mod.null,test='Rao')[2,'Pr(>Chi)'] # model-based score test by package
  vs_san = A %*% D_san %*% t(A)
  Score_san = U3^2/vs_san
  p.score.san = 1 - pchisq(Score_san,df = 1) # sandwich score test
  return(list(pval.model = p.score.mb, pval.sandwich = p.score.san))
}

rulefit.test <- function(x1, x2, y, nullstats){
  n <- length(y)
  dat <- data.frame(y=y, x1=x1, x2=x2) # data frame to be passed to the function pre()
  if(length(table(y))==2){family="binomial"} else {family="gaussian"}
  if(family=="binomial") {dat$y <- factor(dat$y)}
  # Fit prediction rule ensemble
  rulefit <- do.call("pre", args=list(formula=y~x1+x2, data=dat, tree.unbiased=F, maxdepth=maxdepth_sampler(),
                                      family=family, sampfrac=min(1, (11*sqrt(n)+1)/n), nfolds=round(min(20, max(0, 5200/n-2)))))
  # observed test statistics
  obsstat <- interact(rulefit, varnames="x1", penalty.par.val="lambda.min", plot=F)
  pval <- mean(obsstat<=nullstats)
  return(list(teststat=obsstat, pval=pval))
}

