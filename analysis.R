rm(list = ls())
load("./NMMAPS_2D_Chi.RData")
source("function.R")
source("BTDLM_rjags.R")
source("BCDLM_rjags.R")
require(rjags)

## UDLIM
mod_UDLM <- glm(Y ~ -1 + Z + X1 + X2 + XI, family = poisson(link = "log"), control = glm.control(epsilon = 1e-10, maxit = 1000),
               na.action = na.omit)
coefmat[, 1] <- coef(mod_UDLM)[-(1:ncov)]
vcovarr[, , 1] <- vcov(mod_UDLM)[-(1:ncov), -(1:ncov)]

## CDLIM
mod_CDLM <- glm(Y ~ -1 + Z + W1 + W2 + WI, family = poisson(link = "log"), control = glm.control(epsilon = 1e-10, maxit = 1000),
               na.action = na.omit)
coefmat[, 2] <- C%*%coef(mod_CDLM)[-(1:ncov)]
vcovarr[, , 2] <- C%*%vcov(mod_CDLM)[-(1:ncov), -(1:ncov)]%*%t(C)

## TDLIM 
Zoffset <- Z%*%coef(mod_CDLM)[1:ncov]
theta1init <- coef(mod_CDLM)[(ncov + 1):(ncov + 1 + df1)]
theta2init <- coef(mod_CDLM)[(ncov + 2 + df1):(ncov + 2 + df1 + df2)]
mod_TDLM <- Tukey(Y, Zoffset, W1, W2, WI, theta1init, theta2init)
coefmat[, 3] <- C%*%mod_TDLM

coefTDLMboot <- matrix(NA, nboot, (df1 + 1) + (df2 + 1) + (df1 + 1)*(df2 + 1))
for(i in 1:nboot) {
	bootsample <- sample(1:length(Y), replace = TRUE)
	Yboot <- Y[bootsample]
	Zboot <- Z[bootsample, ]
	W1boot <- W1[bootsample, ]
	W2boot <- W2[bootsample, ]
	WIboot <- WI[bootsample, ]
	CDLMboot <- glm(Yboot ~ -1 + Zboot + W1boot + W2boot + WIboot, family = poisson(link = "log"), 
		control = glm.control(epsilon = 1e-10, maxit = 1000), na.action = na.omit)
	Zoffsetboot <- Zboot%*%coef(CDLMboot)[1:ncov]
	theta1initboot <- coef(CDLMboot)[(ncov + 1):(ncov + 1 + df1)]
	theta2initboot <- coef(CDLMboot)[(ncov + 2 + df1):(ncov + 2 + df1 + df2)]
	TDLMboot <- Tukey(Yboot, Zoffsetboot, W1boot, W2boot, WIboot, theta1initboot, theta2initboot)
	coefTDLMboot[i, ] <- TDLMboot
}

# BTDLIM
bugsfiles <- paste("method", methodindex, c("Modeldp.BUGS", "BUGS_data.txt", "BUGS_inits.txt"), sep = "")
theta1Tukey <- mod_TDLM[1:(df1 + 1)]
theta2Tukey <- mod_TDLM[(df1 + 2):(df1 + df2 + 2)]
etainitBTDLM <- rep((mod_TDLM[-(1:(df1 + df2 + 2))] / (theta1Tukey %x% theta2Tukey))[1], nlag1*nlag2)
mod_BTDLM <- BTDLM(Y = Y, X1 = X1, X2 = X2, XI = XI, Z = Z, C1 = C1, C2 = C2, corinv_array = corinv_array, n.thin = 10, n.burnin = 1000, idx = idx, inits = list(theta1 = theta1Tukey, theta2 = theta2Tukey, eta = etainitBTDLM, alpha = alphainit)) 
coefmat[, 4] <- apply(mod_BDLM[, (ncov + 1):(ncov + nlag1 + nlag2 + nlag1*nlag2)], 2, mean)

# BCDLIM
bugsfiles <- paste("method", methodindex, c("Modeldp.BUGS", "BUGS_data.txt", "BUGS_inits.txt"), sep = "")
alphainit <- coef(mod_UDLM)[1:ncov]
beta10init <- coef(mod_UDLM)[(ncov + 1):(ncov + df1 + 1)]
beta11init <- coef(mod_UDLM)[(ncov + df1 + 2):(ncov + nlag1)]
beta20init<- coef(mod_UDLM)[(ncov + nlag1 + 1):(ncov + nlag1 + df2 + 1)]
beta21init<- coef(mod_UDLM)[(ncov + nlag1 + df2 + 2):(ncov + nlag1 + nlag2)]
betaI0init<- coef(mod_UDLM)[(ncov + nlag1 + nlag2 + 1):(ncov + nlag1 + nlag2 + (df1 + 1)*(df2 + 1))]
betaI1init<- coef(mod_UDLM)[-(1:(ncov + nlag1 + nlag2 + (df1 + 1)*(df2 + 1)))]

mod_BCDLM <- BCDLM_rjags(Y = Y, W1_up = W1_up, W2_up = W2_up, WI_up = WI_up, W1_p = W1_p, W2_p = W2_p, WI_p = WI_p, Z = Z, S1inv = S1inv, S2inv = S2inv, SIinv = SIinv, n.thin = 10, n.burnin = 1000, idx = idx, inits = list(alpha = alphainit, beta10 = beta10init, beta11 = beta11init, beta20 = beta20init, beta21 = beta21init, betaI0 = betaI0init, betaI1 = betaI1init), bugsfiles = bugsfiles)
coefmat[, 5] <- apply(mod_BCDLM[, (ncov + 1):(ncov + nlag1 + nlag2 + nlag1*nlag2)], 2, mean)
