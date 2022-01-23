model {
  
  ## Priors
  CL ~ dnorm(15, 0.04)
  beta0 ~ dnorm(0, 0.25)
  beta1 ~ dnorm(0, 0.25)
  beta2 ~ dnorm(0, 0.25)
  beta3 ~ dnorm(0, 0.25)
  beta4 ~ dnorm(0, 0.25)
  beta5 ~ dnorm(0, 0.25)
  betaN ~ dnorm(0, 0.25)

  ## Likelihood
  for(i in 1:nsite) {
    SR[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + step(N[i]-CL) * betaN * (N[i]-CL)
  }
  
  for(i in 1:40) {
    SRpre[i] <- exp(beta0 + betaN * step((i-1) - CL) * ((i-1) - CL))
  }
  
}