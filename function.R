Tukey <- function(Y, Zoffset, W1, W2, WI, theta1, theta2) {
	W1offset <- W1%*%theta1
	W2offset <- W2%*%theta2
	logliknew <- -1e10
	delta <- 1
	while(delta > 1e-5) {

		WIstar <- WI%*%(theta1 %x% theta2)
		mod_eta <- glm(Y ~ -1 + offset(Zoffset) + offset(W1offset) + offset(W2offset) + WIstar, family = poisson(link = "log"), 
			control = glm.control(epsilon = 1e-10, maxit = 1000), na.action = na.omit)
		eta <- coef(mod_eta)
		W2star <- W2 + eta*t(apply(WI, 1, function(x) apply(matrix(x*rep(theta1, each = df2 + 1), df2 + 1, df1 + 1), 1, sum)))
		mod_theta2 <- glm(Y ~ -1 + offset(Zoffset) + offset(W1offset) + W2star, family = poisson(link = "log"), 
			control = glm.control(epsilon = 1e-10, maxit = 1000), na.action = na.omit)
		theta2 <- coef(mod_theta2)
		W2offset <- W2%*%theta2

		W1star <- W1 + eta*t(apply(WI, 1, function(x) apply(matrix(x*rep(theta2, df1 + 1), df2 + 1, df1 + 1), 2, sum)))
		mod_theta1 <- glm(Y ~ -1 + offset(Zoffset) + offset(W2offset) + W1star, family = poisson(link = "log"), 
			control = glm.control(epsilon = 1e-10, maxit = 1000), na.action = na.omit)
		theta1 <- coef(mod_theta1)
		W1offset <- W1%*%theta1

		loglikold <- logliknew
		logliknew <- as.numeric(logLik(mod_theta1))
		delta <- logliknew - loglikold
	}
	return(c(theta1, theta2, eta*(theta1 %x% theta2)))
}
ConditionCor <- function(MarginalCor, ind) {
	R11 <- MarginalCor[which(ind != 0), which(ind != 0)]
	R12 <- MarginalCor[which(ind != 0), which(ind == 0)]
	R22 <- MarginalCor[which(ind == 0), which(ind == 0)]
	return(R11 - R12%*%solve(R22)%*%t(R12))
}
