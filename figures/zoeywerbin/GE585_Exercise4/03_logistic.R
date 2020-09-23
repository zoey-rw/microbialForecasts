##' Logistic model
##'
##' @param theta  parameter vector
##' @param x      vector of x values
##' @return vector of model predictions
pred.logistic <- function(theta,x){
  z = exp(theta[3]+theta[4]*x)
  Ey = theta[1]+theta[2]*z/(1+z) 
}

##' Fit logistic model
##' 
##' @param dat  dataframe of day of year (doy), gcc_mean, gcc_std
##' @param par  vector of initial parameter guess
##' @return  output from numerical optimization
fit.logistic <- function(dat,par){
  
  ## define log likelihood
  lnL.logistic <- function(theta,dat){
    -sum(dnorm(dat$gcc_mean,pred.logistic(theta,dat$doy),dat$gcc_std,log=TRUE))
  }
  
  ## fit by numerical optimization
  optim(par,fn = lnL.logistic,dat=dat)
}
