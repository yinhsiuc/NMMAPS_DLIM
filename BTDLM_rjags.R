BTDLM <- function(Y, X1, X2, XI, Z, C1, C2, corinv_array, n.sample = 1000, writemodel = 1, n.burnin = 20000, n.thin = 1, bugsfiles = c("Modeldp.BUGS", "BUGS_data.txt", "BUGS_inits.txt"), writedatainits=TRUE, idx, inits = inits) {
  if(writemodel) {
    cat("
        model
        {
	alpha ~ dmnorm(alphamu, 0.0001*alphaSigma)
	theta1 ~ dmnorm(mutheta1, 0.0001*I1)
	theta2 ~ dmnorm(mutheta2, 0.0001*I2) 
	beta1 <- C1%*%theta1
	beta2 <- C2%*%theta2
	psiind ~ dcat(ppsi[])
	eta ~ dmnorm(mueta, 0.0001*corinv_array[, , psiind])
	for(d1 in 1:nlag1) {
		for(d2 in 1:nlag2) {
			betaI[nlag2*(d1 - 1) + d2] <- beta1[d1]*beta2[d2]*eta[nlag2*(d1 - 1) + d2]
		}
	}
        for(iter1 in 1:n){
		muY[iter1] <- exp(Z[iter1, ]%*%alpha + X1[iter1, ]%*%beta1 + X2[iter1, ]%*%beta2 + XI[iter1, ]%*%betaI)
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
  nlag1 <- dim(X1)[2]
  nlag2 <- dim(X2)[2] 
  dtheta1 <- dim(C1)[2]
  dtheta2 <- dim(C2)[2]
  mutheta1 <- rep(0, dtheta1)
  mutheta2 <- rep(0, dtheta2)
  penaltyscale <- 1
  I1 <- diag(dtheta1) 
  I2 <- diag(dtheta2)
  npsi <- dim(corinv_array)[3]
  ppsi <- rep(1 / npsi, npsi)
  mueta <- rep(0, nlag1*nlag2)
  pl <- pairlist(n = n, Y = Y, X1 = X1, X2 = X2, XI = XI, Z = Z, ppsi = ppsi, corinv_array = corinv_array, nlag1 = nlag1, nlag2 = nlag2, mutheta1 = mutheta1, mutheta2 = mutheta2, I1 = I1, I2 = I2, C1 = C1, C2 = C2, mueta = mueta, alphamu = alphamu, alphaSigma = alphaSigma)
  jags <- jags.model(bugsfiles[1], data = pl, n.chains = 2, inits = inits)
  para <- c("beta1", "beta2", "betaI")
  update(jags, n.burnin)
  PS <- (coda.samples(jags, para, n.iter = n.sample*n.thin, thin = n.thin)[[1]])
  return(PS) 
  }
