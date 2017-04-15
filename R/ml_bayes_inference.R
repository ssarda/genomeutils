#' @title Gaussian MLE parameters and Model selection 
#' @description This function accepts a vector with all values of a given unknown distribution
#' It estimates the maximum likelihood parameters by fitting a guassian distribution over it.
#' General-purpose optimization to estimate parameters is based on Nelder–Mead algorithm
#' The Akaike Information Criterion is also output to enable downstream model selection
#' @param X A Vector of values.
#' @param start_mu Initial value to begin optimization at for mean of gaussian distribution
#' @param start_sigma Initial value to begin optimization at for standard deviation of gaussian distribution 
#' @return A vector containing the optimized MLE mean and standard deviation as well as the log likelihood and AIC to enable downstream model selection
#' @export
mle_gaussian_dist = function(X,start_mu,start_sigma)
{
  LL.normal <- function(p,X)
  {
  mu = p['mu']
  sigma = p['sigma']
  -sum(dnorm(X, mean = mu, sd = sigma, log = TRUE))
  }
  
  p0.normal = c(mu=start_mu,sigma=start_sigma)

  mle.normal = optim(p0.normal, LL.normal, X=X,hessian=T)
  
  AIC.normal = 2*length(mle.normal$par) + (2*mle.normal$value)
  
  info.fit = cbind(mu=mle.normal$par[1], sigma=mle.normal$par[2], LogLik.normal=-mle.normal$value, AIC.normal=AIC.normal)
  
  return(info.fit)
}




#' @title Weibull MLE parameters and Model selection 
#' @description This function accepts a vector with all values of a given unknown distribution
#' It estimates the maximum likelihood parameters by fitting a weibull distribution over it.
#' General-purpose optimization to estimate parameters is based on Nelder–Mead algorithm
#' The Akaike Information Criterion is also output to enable downstream model selection
#' @param X A Vector of values.
#' @param start_k Initial value to begin optimization at for shape of weibull distribution
#' @param start_lambda Initial value to begin optimization at for sscale of weibull distribution 
#' @return A vector containing the optimized MLE mean and standard deviation as well as the log likelihood and AIC to enable downstream model selection
#' @export
mle_weibull_dist = function(X,start_k,start_lambda)
{
 
  LL.weibull <- function(p,X)
  {
  k=p['k']
  lambda=p['lambda']
  -sum(dweibull(X, shape = k, scale = lambda, log = TRUE))
  }
  
  p0.weibull = c(k=start_k,lambda=start_lambda)

  mle.weibull = optim(p0.weibull, LL.weibull, X=X,hessian=T)
  
  AIC.weibull = 2*length(mle.weibull$par) + (2*mle.weibull$value)
  

  info.fit = cbind(k=mle.weibull$par[1], lambda=mle.weibull$par[2], LogLik.weibull=-mle.weibull$value, AIC.weibull=AIC.weibull)
  
  return(info.fit)
}




#' @title Implementing binomial GLM (logistic regression) from scratch. 
#' @description This function accepts a covariate X and response Y and uses optimization to estimate the maximum likelihood parameters of the binomial GLM fit.
#' General-purpose optimization to estimate parameters is based on Nelder–Mead algorithm
#' The Akaike and Bayesian Information Criterion is also output to enable downstream model selection
#' @param X A Covariate vector of values.
#' @param Y A Response vector of values.
#' @return A vector containing the optimized MLE parameters of the fit as well as the AIC and BIC to enable downstream model selection
#' @export
logistic.model = function(X,Y)
{
  getP.hat <- function(beta0, beta1, X){
  lP <- (beta0 + beta1 * X)
  exp(lP)/(1 + exp(lP))
  }

  logLike.Binomial <- function(coefs, X, Y){
  beta0 <- coefs['beta0']
  beta1 <- coefs['beta1']
  p.hat <- getP.hat(beta0, beta1, X)
  -sum(dbinom(Y, size = 1, prob = p.hat, log = TRUE))
  }

  coef0 <- c(beta0 = 0, beta1 = -1)
  covariate.fit <- optim(coef0, logLike.Binomial, X=X, Y=Y)

  params.covariate = length(covariate.fit$par)
  AIC.covariate = (2*params.covariate) + (2*covariate.fit$value)
  BIC.covariate = log(n)*params.covariate + (2*covariate.fit$value)
  
  out = c(covariate.fit$par[1],covariate.fit$par[2],AIC.covariate,BIC.covariate)
  return(out)
}



#' @title Bayesian estimation for a mixture of betas.
#' @description This function returns the bayesian posterior estimate for a mixture of beta distributions
#' @param x Number of observations per beta distribution
#' @param p Number of betas
#' @param a A vector of shape 1 parameters for all beta distributions
#' @param b A vector of shape 2 parameters for all beta distributions
#' @return The combined posterior estimate of the mixtures
#' @export
mix.beta<-function(x,p,a,b)
{
   dense<-0*x
   for(j in 1:length(p))
   {
      dense<-dense+p[j]*dbeta(x,a[j],b[j])
   }
	dense
}


#' @title Bayesian estimation for beta binomial distribution
#' @description This function returns the bayesian posterior estimate for a beta binomial distribution
#' @param n Total number of binomial draws
#' @param y The total number of successes
#' @param a Shape 1 parameter of beta distribution
#' @param b Shape 1 parameter of beta distribution
#' @return The likelihood and posterior of the beta binmomial distribution
#' @export
beta_binom<-function(n,y,a=1,b=1)
{
    #likelihood: y|theta~binom(n,theta)
    #prior: theta~beta(a,b)
    #posterior: theta|y~beta(a+y,n-y+b)

    theta<-seq(0.001,0.999,0.001)
    prior<-dbeta(theta,a,b)
    if(n>0){likelihood<-dbinom(rep(y,length(theta)),n,theta)}
    if(n>0){posterior<-dbeta(theta,a+y,n-y+b)}

    #standardize!
    prior<-prior/sum(prior)
    if(n>0){likelihood<-likelihood/sum(likelihood)}
    if(n>0){posterior<-posterior/sum(posterior)}
	return(c(likelihood, posterior))
}