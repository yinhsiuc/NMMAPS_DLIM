BCDLM_rjags <- function(Y, W1_up, W2_up, WI_up, W1_p, W2_p, WI_p, Z, S1inv, S2inv, SIinv, n.sample = 1000, writemodel = 1, n.burnin = 20000, n.thin = 1, bugsfiles = c("Modeldp.BUGS", "BUGS_data.txt", "BUGS_inits.txt"), writedataintits = TRUE, idx = idx, inits = inits) {
  if(writemodel) {
    cat("
        model
        {
    	tau1 ~ dgamma(0.001, 0.001)
    	tau2 ~ dgamma(0.001, 0.001) 
    	tauI ~ dgamma(0.001, 0.001) 
    	alpha ~ dmnorm(alphamu, 0.0001*alphaSigma)
        beta10 ~ dmnorm(beta10mu, 0.0001*beta10Sigma)
        beta20 ~ dmnorm(beta20mu, 0.0001*beta20Sigma)
        betaI0 ~ dmnorm(betaI0mu, 0.0001*betaI0Sigma)

        beta11 ~ dmnorm(beta11mu, tau1*S1inv)
        beta21 ~ dmnorm(beta21mu, tau2*S2inv)
        betaI1 ~ dmnorm(betaI1mu, tauI*SIinv)

        for(iter1 in 1:n){
		muY[iter1] <- exp(Z[iter1, ]%*%alpha + W1_up[iter1, ]%*%beta10 + W1_p[iter1, ]%*%beta11 + W2_up[iter1, ]%*%beta20 + W2_p[iter1, ]%*%beta21 + WI_up[iter1, ]%*%betaI0 + WI_p[iter1, ]%*%betaI1)
        }
        
        for(iter1 in 1:n){
          Y[iter1] ~ dpois(muY[iter1])
        }
        
        }
        ", file=bugsfiles[1]);
  }
  n <- length(Y)
  p <- dim(Z)[2]
  alphamu <- rep(0, p)
  alphaSigma <- diag(p)
  up1 <- dim(W1_up)[2]
  up2 <- dim(W2_up)[2]
  upI <- dim(WI_up)[2]
  p1 <- dim(W1_p)[2]
  p2 <- dim(W2_p)[2]
  pI <- dim(WI_p)[2]
  beta10mu <- rep(0, up1)
  beta20mu <- rep(0, up2)
  betaI0mu <- rep(0, upI)
  beta11mu <- rep(0, p1)
  beta21mu <- rep(0, p2)
  betaI1mu <- rep(0, pI)
  beta10Sigma <- diag(up1)
  beta20Sigma <- diag(up2)
  betaI0Sigma <- diag(upI)
  pl <- pairlist(n = n, Y = Y, W1_up = W1_up, W2_up = W2_up, WI_up = WI_up, W1_p = W1_p, W2_p = W2_p, WI_p = WI_p, Z = Z, S1inv = S1inv, S2inv = S2inv, SIinv = SIinv, alphamu = alphamu, alphaSigma = alphaSigma, beta10mu = beta10mu, beta20mu = beta20mu, betaI0mu = betaI0mu, beta10Sigma = beta10Sigma, beta20Sigma = beta20Sigma, betaI0Sigma = betaI0Sigma, beta11mu = beta11mu, beta21mu = beta21mu, betaI1mu = betaI1mu) 
  jags <- jags.model(bugsfiles[1], data = pl, n.chains = 2, inits = inits)
  para <- c("beta10", "beta11", "beta20", "beta21", "betaI0", "betaI1")
  update(jags, n.burnin)
  PS <- (coda.samples(jags, para, n.iter = n.sample*n.thin, thin = n.thin)[[1]])
  return(PS) 
  }
