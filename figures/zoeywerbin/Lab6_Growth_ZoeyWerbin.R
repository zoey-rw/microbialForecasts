## ------------------------------------------------------------------------
load("data/Lab7.RData")
library(mvtnorm)
library(coda)


## ------------------------------------------------------------------------
theta=.4
sg=20
beta=c(2,70)

fun <- function(L) (beta[1] + beta[2]*(L/(L+theta)))
ci_hi <- function(L) (beta[1] + beta[2]*(L/(L+theta)) + 1.96*sqrt(sg))
ci_low <- function(L) (beta[1] + beta[2]*(L/(L+theta)) - 1.96*sqrt(sg))
plot(grow ~ L)
curve(fun, add = T)
curve(ci_hi, col = "green", add = T)
curve(ci_low,  col = "green", add = T)


## ------------------------------------------------------------------------
b0 <- as.vector(c(0,0))
vinvert <- solve(diag(1000,2))
s1 <- 0.1
s2 <- 0.1


## ------------------------------------------------------------------------
a1 = 1.91
a2 = 10.17


## ------------------------------------------------------------------------
Vb <- vinvert %*% b0


## ------------------------------------------------------------------------
##storage for MCMC
ngibbs <- 10    			## number of updates
bgibbs <- matrix(0.0,nrow=ngibbs,ncol=2) 	## storage for beta
sgibbs <- rep(sg,ngibbs)			## storage for sigma2
tgibbs <- rep(theta,ngibbs)   ## storage for theta


## ------------------------------------------------------------------------
sinv = 1/sg
n <- length(L)
z <- L/(L+theta)
X <- cbind(rep(1,n),z)


## ------------------------------------------------------------------------
## jump
dtnorm <- function(x,mu,sd){
  y = dnorm(x,mu,sd,log=TRUE)-log(pnorm(1,mu,sd)-pnorm(0,mu,sd))
  y[x<0 | x > 1] = -Inf
  return(y)
}
xseq = seq(-0.5,1,length=100)
plot(xseq,exp(dtnorm(xseq,0.25,0.3)),type='l')
lines(xseq,dnorm(xseq,0.25,0.3),col=2)


## ------------------------------------------------------------------------
rtnorm <- function(n,mu,sd){
  x <- rnorm(n,mu,sd)
  sel <- which(x < 0 | x > 1)
  while(length(sel)> 0){
    x[sel] <- rnorm(length(sel),mu,sd)
    sel <- which(x < 0 | x > 1)
  }
  return(x)
}


## ------------------------------------------------------------------------
JumpSD <- 0.1


## ------------------------------------------------------------------------
## sample regression parameters
  bigV    <- solve(sinv*crossprod(X) + vinvert)
  littlev <- sinv*crossprod(X,grow) + Vb
  b <- t(rmvnorm(1,bigV %*% littlev,bigV))


## ------------------------------------------------------------------------
  ## sample variance
  u1 <- s1 + n/2
  u2 <- s2 + 0.5*crossprod(grow-X%*%b)
  sinv <- rgamma(1,u1,u2)
  sg <- 1/sinv


## ------------------------------------------------------------------------
 ##theta
  tnew <- rtnorm(1,theta,JumpSD)  		##propose new theta
  znew <- L/(L+tnew)					## calculate new z
  Xnew <- cbind(rep(1,n),znew)				## calculate new X
  anum <- dmvnorm(grow,Xnew%*%b,diag(sg,n),log=TRUE) + 	##likelihood
	        dbeta(tnew,a1,a2,log=TRUE)			##prior
  jnum <- dtnorm(tnew,theta,JumpSD)				##jump
  aden <- dmvnorm(grow,X%*%b,diag(sg,n),log=TRUE) +	##likelihood
		      dbeta(theta,a1,a2,log=TRUE)			##prior
  jden <- dtnorm(theta,tnew,JumpSD)				##jump
  a <- exp((anum-jnum)-(aden-jden))			## acceptance criteria
  if(a > runif(1)){					## accept with probability a
    theta <- tnew						## update theta if step accepted
    X <- Xnew						## update X if step accepted
  }


## ------------------------------------------------------------------------
# I made this into a function for ease of re-running
print("test")
run_mcmc <- function(ngibbs, JumpSD, thin=1) {
  require(coda)
    accept <- 0
    reject <- 0
    
    bgibbs <- matrix(0.0,nrow=ngibbs,ncol=2) 	## storage for beta
    sgibbs <- rep(sg,ngibbs)			## storage for sigma2
    tgibbs <- rep(theta,ngibbs)   ## storage for theta
    
    ## Gibbs loop
    for(g in 1:ngibbs){
    
      ## sample regression parameters
      bigV    <- solve(sinv*crossprod(X) + vinvert)
      littlev <- sinv*crossprod(X,grow) + Vb
      b <- t(rmvnorm(1,bigV %*% littlev,bigV))
    
      ## sample variance - inverse gamma 
      u1 <- s1 + n/2
      u2 <- s2 + 0.5*crossprod(grow-X%*%b)
      sinv <- rgamma(1,u1,u2)
      sg <- 1/sinv
    
      ## Sample theta
       ##theta
      tnew <- rtnorm(1,theta,JumpSD)  		##propose new theta
      znew <- L/(L+tnew)					## calculate new z
      Xnew <- cbind(rep(1,n),znew)				## calculate new X
      anum <- dmvnorm(grow,Xnew%*%b,diag(sg,n),log=TRUE) + 	##likelihood
    	        dbeta(tnew,a1,a2,log=TRUE)			##prior
      jnum <- dtnorm(tnew,theta,JumpSD)				##jump
      aden <- dmvnorm(grow,X%*%b,diag(sg,n),log=TRUE) +	##likelihood
    		      dbeta(theta,a1,a2,log=TRUE)			##prior
      jden <- dtnorm(theta,tnew,JumpSD)				##jump
      a <- exp((anum-jnum)-(aden-jden))			## acceptance criteria
      if(a > runif(1)){					## accept with probability a
        theta <- tnew						## update theta if step accepted
        X <- Xnew						## update X if step accepted
        accept <- accept + 1
      } else {
        reject <- reject + 1
      }
        
      ## storage
      bgibbs[g,] <- b  ## store the current value of beta vector
      sgibbs[g]  <- sg	## store the current value of the variance
      tgibbs[g]  <- theta
    
      if(g %%1000 == 0) print(g)	##counter to show how many steps have been performed
    }
    acceptance_rate <- accept/(accept+reject) 
    print(paste0("acceptance rate: ", acceptance_rate,"%"))
    return(list(as.mcmc(bgibbs, thin = thin), 
                as.mcmc(sgibbs, thin = thin), 
                as.mcmc(tgibbs, thin = thin), 
                acceptance_rate, JumpSD))
}


## ------------------------------------------------------------------------

# Adjusting the Jump SD to get a good mixing rate

# t1 <- run_mcmc(ngibbs = 10000, JumpSD = .005)
# print(t1[[4]]) # 88% acceptance
# 
# 
# t2 <- run_mcmc(ngibbs = 1000, JumpSD = .2)
# print(t2[[4]]) # 10% acceptance
# 
# 
# t3 <- run_mcmc(ngibbs = 10000, JumpSD = .05)
# print(t3[[4]]) # 27.8% acceptance

t4 <- run_mcmc(ngibbs = 10000, JumpSD = .04)
print(t4[[4]]) # 33% acceptance



## ------------------------------------------------------------------------
t4.chain2 <- run_mcmc(ngibbs = 10000, JumpSD = .04)
t4.chain3 <- run_mcmc(ngibbs = 10000, JumpSD = .04)

# Combine each chain into one mcmc.list
beta.mcmc <- as.mcmc.list(list(t4[[1]], t4.chain2[[1]], t4.chain3[[1]]))
sg.mcmc <- as.mcmc.list(list(t4[[2]], t4.chain2[[2]], t4.chain3[[2]]))
theta.mcmc <- as.mcmc.list(list(t4[[3]], t4.chain2[[3]], t4.chain3[[3]]))

# View trace plots
plot(beta.mcmc)
plot(sg.mcmc)
plot(theta.mcmc)


## ------------------------------------------------------------------------
# Check convergence
gelman.plot(beta.mcmc)
gelman.plot(sg.mcmc)
gelman.plot(theta.mcmc)


## ------------------------------------------------------------------------
# Remove burn in
burnin = 1000                             
beta.mcmc.burn <- window(beta.mcmc, start=burnin)
sg.mcmc.burn <- window(sg.mcmc, start=burnin)
theta.mcmc.burn <- window(theta.mcmc, start=burnin)


# Now that we've removed burn-in, we can view autocorrelation plots
acfplot(beta.mcmc.burn, lag.max = 100)
acfplot(sg.mcmc.burn, lag.max = 100)
acfplot(theta.mcmc.burn, lag.max = 100)

effectiveSize(beta.mcmc.burn) # check effective sample sizes, which are only about 300 for our first beta
effectiveSize(sg.mcmc.burn) # check effective sample sizes, which are ~21k
effectiveSize(theta.mcmc.burn) # check effective sample sizes, which are only about 200 for theta


## ------------------------------------------------------------------------
t4.chain1 <- run_mcmc(ngibbs = 100000, JumpSD = .04)
t4.chain2 <- run_mcmc(ngibbs = 100000, JumpSD = .04)
t4.chain3 <- run_mcmc(ngibbs = 100000, JumpSD = .04)

# Combine each chain into one mcmc.list
beta.mcmc <- as.mcmc.list(list(t4.chain1[[1]], t4.chain2[[1]], t4.chain3[[1]]))
sg.mcmc <- as.mcmc.list(list(t4.chain1[[2]], t4.chain2[[2]], t4.chain3[[2]]))
theta.mcmc <- as.mcmc.list(list(t4.chain1[[3]], t4.chain2[[3]], t4.chain3[[3]]))
# Thin and remove burn-in
beta.thin = window(beta.mcmc,thin=70,start=burnin)
sg.thin =  window(sg.mcmc,thin=70,start=burnin)
theta.thin =  window(theta.mcmc,thin=70,start=burnin)

effectiveSize(beta.thin) # check effective sample sizes
effectiveSize(sg.thin) # check effective sample sizes
effectiveSize(theta.thin) # check effective sample sizes


## ------------------------------------------------------------------------
summary(beta.thin)
summary(sg.thin)
summary(theta.thin)


## ------------------------------------------------------------------------
## credible and prediction intervals
xpred <- seq(0,1,length=30)
npred <- length(xpred)
ypred <- matrix(NA,nrow=ngibbs,ncol=npred)
ycred <- matrix(NA,nrow=ngibbs,ncol=npred)

for(g in 1:ngibbs){
  Ey <- bgibbs[g,1] + bgibbs[g,2] * xpred/(xpred + tgibbs[g])
  ycred[g,] <- Ey
  ypred[g,] <- rnorm(npred,Ey,sqrt(sgibbs[g]))
}
ci <- apply(ycred,2,quantile,c(0.025,0.5,0.975))
pi <- apply(ypred,2,quantile,c(0.025,0.975))

plot(L,grow)
lines(xpred,ci[2,],col=3,lwd=2)  ## median model
lines(xpred,ci[1,],col=3,lty=2)	## model CI
lines(xpred,ci[3,],col=3,lty=2)
lines(xpred,pi[1,],col=4,lty=2)	## model PI
lines(xpred,pi[2,],col=4,lty=2)


## ------------------------------------------------------------------------
theta_regression <- "
model{

  b ~ dmnorm(b0,Vb)  	## multivariate Normal prior on vector of regression params
  S ~ dgamma(s1,s2)    ## prior precision
  theta ~ dbeta(1.91, 10.17) ## prior on theta - shape values provided by Mike

  for(i in 1:n){
	  mu[i] <- b[1] + b[2]*(x[i]/(x[i]+theta))   	## process model
	  y[i]  ~ dnorm(mu[i],S)		        ## data model
  }
}
"
data <- list(x = L, y = grow, n = length(grow))
## specify priors
data$b0 <- as.vector(c(0,0))      ## regression b means
data$Vb <- solve(diag(10000,2))   ## regression b precisions
data$s1 <- 0.1                    ## error prior n/2
data$s2 <- 0.1                    ## error prior SS/2
## initial conditions
inits <- list()
for(i in 1:3){
 inits[[i]] <- list(b = rnorm(2,0,5), S = runif(1,1/200,1/20))
}
j.model   <- rjags::jags.model (file = textConnection(theta_regression),
                             data = data,
                             inits = inits,
                             n.chains = 3)


## ------------------------------------------------------------------------
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("b","S","theta"),
                                n.iter = 5000)


## ------------------------------------------------------------------------
plot(jags.out)
summary(jags.out)

