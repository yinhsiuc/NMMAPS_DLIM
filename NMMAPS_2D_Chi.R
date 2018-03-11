rm(list = ls())
require(dlnm)
require(Matrix)
require(splines)
require(fields)
require(magic)

## read-in data
load("./NMMAP-CHI.RData")
source("function.R")

## data cleaning
nmonth <- 5114
time <- rep(seq(-nmonth/2 + 0.5, length = nmonth), 3)
addlag_pm <- addlag_o3 <- data.frame(matrix(NA, nmonth, 7))
colnames(addlag_pm) <- paste("l", 8:14, "pm10tmean", sep = "")
colnames(addlag_o3) <- paste("l", 8:14, "o3tmean", sep = "")
addlag_pm[8:nmonth, 1:7] <- chi.dat[1:(nmonth-7), seq(208, 244, 6)]
addlag_o3[8:nmonth, 1:7] <- chi.dat[1:(nmonth-7), seq(213, 249, 6)]
addlags_pm <- rbind(addlag_pm, addlag_pm, addlag_pm)
addlags_o3 <- rbind(addlag_o3, addlag_o3, addlag_o3)
age2ind <- 1*(chi.dat$agecat == 2)
age3ind <- 1*(chi.dat$agecat == 3)
dat <- cbind(time, chi.dat[, c(3:4, 8, 13, 17, 198:199, 32, seq(208, 244, 6))], addlags_pm, chi.dat[, c(82, seq(213, 249, 6))], addlags_o3, age2ind, age3ind)
dat$dow <- as.factor(dat$dow)
dat$agecat <- as.factor(dat$agecat)
dat1 <- dat[complete.cases(dat), ]
n <- dim(dat1)[1]
T <- dat1$time - min(dat$time)

## construct design matrices
df1 <- df2 <- 4
L1 <- L2 <- 14
nlag1 <- L1 + 1
nlag2 <- L2 + 1
C1 <- onebasis(0:L1, fun = "bs", degree = df1, intercept = TRUE)
C2 <- onebasis(0:L2, fun = "bs", degree = df2, intercept = TRUE)
CI <- C1 %x% C2
class(C1) <- "matrix"
class(C2) <- "matrix"
C <- as.matrix(bdiag(C1, C2, CI))

scalefac <- 1000
X1 <- as.matrix(dat1[, 9:(9+L1)]) / scalefac
minX1 <- min(X1)
X1 <- X1 - min(X1)
X2 <- as.matrix(dat1[, 24:(24+L2)]) / scalefac
minX2 <- min(X2)
X2 <- X2 - min(X2)
XI <- t(apply(cbind(X1, X2), 1, function(x) x[1:(1+L1)] %x% x[-(1:(1 + L1))]))
colnames(XI) <- paste("int", 1:dim(XI)[2], sep = "")

W1 <- X1%*%C1
W2 <- X2%*%C2
WI <- XI%*%CI

colnames(W1) <- paste("b1", 1:(df1 + 1), sep = "")
colnames(W2) <- paste("b2", 1:(df2 + 1), sep = "")
colnames(WI) <- paste("bi", 1:((df1 + 1)*(df2 + 1)), sep = "")

Zformu <- formula("death ~ dow + agecat + ns(tmpd, 6) + ns(dptp, df = 3) + ns(rmtmpd, df = 6) + ns(rmdptp, df = 3) + ns(time, df = 98) + I(ns(time, df = 14)*age2ind) + I(ns(time, df = 14)*age3ind)")
Z <- model.matrix(Zformu, dat1)
Y <- dat1$death

## summary
ncov <- dim(Z)[2]

# TDLM
nboot <- 500

## for Gaussian process
M1 <- 20
M2 <- 20
a1 <- a2 <- aI <- 2
b1 <- b2 <- bI <- 1
nu1 <- nu2 <- nuI <- 3
psi1 <- psi2 <- psiI <- 0.6

grid_dim1 <- grid_dim2 <- 8
dim_all <- (M1 + 1) + (M2 + 1) + grid_dim1 + grid_dim2 + (M1 + 1)*(M2 + 1) + grid_dim1*grid_dim2
an1 <- a1 + grid_dim1 / 2
an2 <- a2 + grid_dim2 / 2
anI <- aI + grid_dim1*grid_dim2 / 2

# GPindbeta1 <- 2:(grid_dim1 + 1)
# GPindbeta2 <- (grid_dim1 + 2):(grid_dim1 + grid_dim2 + 1)
# GPindbetaI <- (grid_dim1 + grid_dim2 + 2):(grid_dim1 + grid_dim2 + grid_dim1*grid_dim2 + 1)

l2 <- cbind(rep(0:M1, each = M2 + 1), rep(0:M2, times = M1 + 1))
l2_inner <- apply(l2, 1, function(x) 1*((x[1] <= L1) && (x[2] <= L2)))

gridwidth1 <- L1 / (grid_dim1 + 1)
gridwidth2 <- L2 / (grid_dim2 + 1)
grid1 <- seq(gridwidth1 / 2, L1 - gridwidth1 / 2,, grid_dim1)
grid2 <- seq(gridwidth2 / 2, L2 - gridwidth2 / 2,, grid_dim2)
gridI <- cbind(rep(grid1, each = grid_dim2), rep(grid2, times = grid_dim1))

pos1 <- c(0:M1, grid1)
Dis1 <- abs(outer(pos1, pos1, "-"))
pos2 <- c(0:M2, grid2)
Dis2 <- abs(outer(pos2, pos2, "-"))
posI <- rbind(l2, gridI)
nposI <- dim(posI)[1]
DisI <- outer(1:nposI, 1:nposI, function(x, y) sqrt((posI[x, 1] - posI[y, 1])^2 + (posI[x, 2] - posI[y, 2])^2))

matern_cor1 <- Matern(Dis1, smoothness = nu1, range = psi1, phi = 1)
matern_cor2 <- Matern(Dis2, smoothness = nu2, range = psi2, phi = 1)
matern_corI <- Matern(DisI, smoothness = nuI, range = psiI, phi = 1)
ind1 <- c(rep(1, (L1 + 1)), rep(0, (M1 - L1)), rep(2, grid_dim1)) 
ind2 <- c(rep(1, (L2 + 1)), rep(0, (M2 - L2)), rep(2, grid_dim2)) 
indI <- c(l2_inner, rep(2, grid_dim1*grid_dim2))

cor_con1 <- ConditionCor(matern_cor1, ind1)
cor_con2 <- ConditionCor(matern_cor2, ind2)
cor_conI <- ConditionCor(matern_corI, indI)

cor_con1star <- cor_con1[-(1:nlag1), -(1:nlag1)]
cor_con2star <- cor_con2[-(1:nlag2), -(1:nlag2)]
cor_conIstar <- cor_conI[-(1:(nlag1*nlag2)), -(1:(nlag1*nlag2))]
cor_mar1star <- cor_con1[1:nlag1, -(1:nlag1)]
cor_mar2star <- cor_con2[1:nlag2, -(1:nlag2)]
cor_marIstar <- cor_conI[1:(nlag1*nlag2), -(1:(nlag1*nlag2))]

invcor_con1star <- solve(cor_con1star)
invcor_con2star <- solve(cor_con2star)
invcor_conIstar <- solve(cor_conIstar)

B1 <- cor_mar1star%*%invcor_con1star
B2 <- cor_mar2star%*%invcor_con2star
BI <- cor_marIstar%*%invcor_conIstar

B <- adiag(B1, B2, BI)
# nGP <- ncol(B)

# BTDLM
rho_grid <- exp(seq(-0.002, -0.01, -0.002)) 
corinv_array <- array(NA, dim = c(nlag1*nlag2, nlag1*nlag2, length(rho_grid)))
coordinate1 <- c(rep(0:L1, each = nlag2))
coordinate2 <- c(rep(0:L2, times = nlag1)) 
Dist <- outer(1:(length(coordinate1)), 1:(length(coordinate1)), function(x, y) sqrt((coordinate1[x] - coordinate1[y])^2 + (coordinate2[x] - coordinate2[y])^2))
for(i in seq_along(rho_grid)) {
	cormat <- rho_grid[i]^Dist
	corinv_array[, , i] <- solve(cormat)  
}
# BCDLM
### newly added
nknots1 <- nlag1 - df1
knotslen1 <- L1  / nknots1
C1full <- onebasis(0:L1, "bs", degree = df1, knots = seq(knotslen1, knotslen1*(nknots1 - 1), knotslen1), intercept = TRUE)
nknots2 <- nlag2 - df2
knotslen2 <- L2  / nknots2
C2full <- onebasis(0:L2, "bs", degree = df2, knots = seq(knotslen2, knotslen2*(nknots2 - 1), knotslen2), intercept = TRUE)
class(C1full) <- "matrix"
class(C2full) <- "matrix"
CIfull <- C1full%x%C2full
ProjC1 <- diag(nlag1) - C1%*%solve(crossprod(C1))%*%t(C1)
C1full <- cbind(C1, svd(ProjC1%*%C1full)$u[, 1:(nlag1 - df1 - 1)])
ProjC2 <- diag(nlag2) - C2%*%solve(crossprod(C2))%*%t(C2)
C2full <- cbind(C2, svd(ProjC2%*%C2full)$u[, 1:(nlag2 - df2 - 1)])
ProjCI <- diag(nlag1*nlag2) - CI%*%solve(crossprod(CI))%*%t(CI)
CIfull <- cbind(CI, svd(ProjCI%*%CIfull)$u[, 1:(nlag1*nlag2 - (df1 + 1)*(df2 + 1))])
Cfull <- adiag(C1full, C2full, CIfull) 

W1full <- X1%*%C1full
W2full <- X2%*%C2full
WIfull <- XI%*%CIfull
W1_up <- W1full[, 1:(df1 + 1)]
W1_p <- W1full[, -(1:(df1 + 1))]
W2_up <- W2full[, 1:(df2 + 1)]
W2_p <- W2full[, -(1:(df2 + 1))]
WI_up <- WIfull[, 1:((df1 + 1)*(df2 + 1))]
WI_p <- WIfull[, -(1:((df1 + 1)*(df2 + 1)))]

SZ <- diag(ncov)
S1 <- diag(rep(0:1, times = c(df1 + 1, nlag1 - df1 - 1)))
S2 <- diag(rep(0:1, times = c(df2 + 1, nlag2 - df2 - 1)))
SI <- diag(rep(0:1, times = c((df1 + 1)*(df2 + 1), nlag1*nlag2 - (df1 + 1)*(df2 + 1))))

S1inv <- diag(nlag1 - df1 - 1)
S2inv <- diag(nlag2 - df2 - 1)
SIinv <- diag(nlag1*nlag2 - (df1 + 1)*(df2 + 1))

save.image(file = "./NMMAPS_2D_Chi.RData")
