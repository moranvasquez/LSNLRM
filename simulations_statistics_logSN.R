#==================================================================
# DESCRIPTION: This file reproduces the tables in the manuscript 
# related to the simulations of the LSNLRM. It is required to load 
# the objects saved in the workspace 'simulations_logSN.RData'
#==================================================================
load("simulations_logSN.RData")
#=====================
# Required packages
#=====================
require(sn)

# 3%, 5%, 25%, 50%, 75%, 95% and 97% quantiles of the 
# standard skew-normal distribution (with true alpha)
q03_stdsn <- qsn(0.03, xi=0, omega=1, alpha=true_alpha_logsn, solver="RFB")
q05_stdsn <- qsn(0.05, xi=0, omega=1, alpha=true_alpha_logsn, solver="RFB")
q25_stdsn <- qsn(0.25, xi=0, omega=1, alpha=true_alpha_logsn, solver="RFB")
q50_stdsn <- qsn(0.50, xi=0, omega=1, alpha=true_alpha_logsn, solver="RFB")  
q75_stdsn <- qsn(0.75, xi=0, omega=1, alpha=true_alpha_logsn, solver="RFB")
q95_stdsn <- qsn(0.95, xi=0, omega=1, alpha=true_alpha_logsn, solver="RFB")
q97_stdsn <- qsn(0.97, xi=0, omega=1, alpha=true_alpha_logsn, solver="RFB")

# Female and male children data (real and simulated) 
# n = 50, 100, 500, 1000
X_50_logsn <- as.data.frame(X_50_logsn)  # n=50
X_50_f_logsn <- X_50_logsn[X_50_logsn$X1_50_logsn == 0, ]  # female
X_50_m_logsn <- X_50_logsn[X_50_logsn$X1_50_logsn == 1, ]  # male
fix.X_50_f_logsn <- apply(X_50_f_logsn, 2, mean)
fix.X_50_m_logsn <- apply(X_50_m_logsn, 2, mean)
X_100_logsn <- as.data.frame(X_100_logsn)  # n=100
X_100_f_logsn <- X_100_logsn[X_100_logsn$X1_100_logsn == 0, ]  # female
X_100_m_logsn <- X_100_logsn[X_100_logsn$X1_100_logsn == 1, ]  # male
fix.X_100_f_logsn <- apply(X_100_f_logsn, 2, mean)
fix.X_100_m_logsn <- apply(X_100_m_logsn, 2, mean)
X_500_logsn <- as.data.frame(X_500_logsn)  # n=500
X_500_f_logsn <- X_500_logsn[X_500_logsn$X1_500_logsn == 0, ]  # female
X_500_m_logsn <- X_500_logsn[X_500_logsn$X1_500_logsn == 1, ]  # male
fix.X_500_f_logsn <- apply(X_500_f_logsn, 2, mean)
fix.X_500_m_logsn <- apply(X_500_m_logsn, 2, mean)
X_1000_logsn <- as.data.frame(X_1000_logsn)  # n=1000
X_1000_f_logsn <- X_1000_logsn[X_1000_logsn$X1_1000_logsn == 0, ]  # female
X_1000_m_logsn <- X_1000_logsn[X_1000_logsn$X1_1000_logsn == 1, ]  # male
fix.X_1000_f_logsn <- apply(X_1000_f_logsn, 2, mean)
fix.X_1000_m_logsn <- apply(X_1000_m_logsn, 2, mean)

# True quartiles (by gender)
true_betas_logsn <- c(true_beta1_logsn,true_beta2_logsn,true_beta3_logsn)
true_q03_f_logsn <- exp(sum(true_betas_logsn*fix.X_f) + true_omega_logsn*q03_stdsn)
true_q05_f_logsn <- exp(sum(true_betas_logsn*fix.X_f) + true_omega_logsn*q05_stdsn)
true_q25_f_logsn <- exp(sum(true_betas_logsn*fix.X_f) + true_omega_logsn*q25_stdsn)
true_q50_f_logsn <- exp(sum(true_betas_logsn*fix.X_f) + true_omega_logsn*q50_stdsn)
true_q75_f_logsn <- exp(sum(true_betas_logsn*fix.X_f) + true_omega_logsn*q75_stdsn)
true_q95_f_logsn <- exp(sum(true_betas_logsn*fix.X_f) + true_omega_logsn*q95_stdsn)
true_q97_f_logsn <- exp(sum(true_betas_logsn*fix.X_f) + true_omega_logsn*q97_stdsn)
true_q03_m_logsn <- exp(sum(true_betas_logsn*fix.X_m) + true_omega_logsn*q03_stdsn)
true_q05_m_logsn <- exp(sum(true_betas_logsn*fix.X_m) + true_omega_logsn*q05_stdsn)
true_q25_m_logsn <- exp(sum(true_betas_logsn*fix.X_m) + true_omega_logsn*q25_stdsn)
true_q50_m_logsn <- exp(sum(true_betas_logsn*fix.X_m) + true_omega_logsn*q50_stdsn)
true_q75_m_logsn <- exp(sum(true_betas_logsn*fix.X_m) + true_omega_logsn*q75_stdsn)
true_q95_m_logsn <- exp(sum(true_betas_logsn*fix.X_m) + true_omega_logsn*q95_stdsn)
true_q97_m_logsn <- exp(sum(true_betas_logsn*fix.X_m) + true_omega_logsn*q97_stdsn)

#=======================
# Sample size n = 50
#=======================

# Estimates of the regression coefficients using the 
# classical QR
beta_qr_q03_Y_50 <- lapply(simulation_logsn_50, "[[", 7)
beta_qr_q05_Y_50 <- lapply(simulation_logsn_50, "[[", 8)
beta_qr_q25_Y_50 <- lapply(simulation_logsn_50, "[[", 9)
beta_qr_q50_Y_50 <- lapply(simulation_logsn_50, "[[", 10)
beta_qr_q75_Y_50 <- lapply(simulation_logsn_50, "[[", 11)
beta_qr_q95_Y_50 <- lapply(simulation_logsn_50, "[[", 12)
beta_qr_q97_Y_50 <- lapply(simulation_logsn_50, "[[", 13)

# Estimates of the 3%, 5%, 25%, 50%, 75%, 95% and 97% 
# quantiles (by gender) using linear quantile regression
estim_q03_qr_50_f_Y <- c()
estim_q05_qr_50_f_Y <- c()
estim_q25_qr_50_f_Y <- c()
estim_q50_qr_50_f_Y <- c()
estim_q75_qr_50_f_Y <- c()
estim_q95_qr_50_f_Y <- c()
estim_q97_qr_50_f_Y <- c()
estim_q03_qr_50_m_Y <- c()
estim_q05_qr_50_m_Y <- c()
estim_q25_qr_50_m_Y <- c()
estim_q50_qr_50_m_Y <- c()
estim_q75_qr_50_m_Y <- c()
estim_q95_qr_50_m_Y <- c()
estim_q97_qr_50_m_Y <- c()
for (i in 1:N) {
  estim_q03_qr_50_f_Y[i] <- sum(beta_qr_q03_Y_50[[i]] * fix.X_50_f_logsn)
  estim_q05_qr_50_f_Y[i] <- sum(beta_qr_q05_Y_50[[i]] * fix.X_50_f_logsn)
  estim_q25_qr_50_f_Y[i] <- sum(beta_qr_q25_Y_50[[i]] * fix.X_50_f_logsn)
  estim_q50_qr_50_f_Y[i] <- sum(beta_qr_q50_Y_50[[i]] * fix.X_50_f_logsn)
  estim_q75_qr_50_f_Y[i] <- sum(beta_qr_q75_Y_50[[i]] * fix.X_50_f_logsn)
  estim_q95_qr_50_f_Y[i] <- sum(beta_qr_q95_Y_50[[i]] * fix.X_50_f_logsn)
  estim_q97_qr_50_f_Y[i] <- sum(beta_qr_q97_Y_50[[i]] * fix.X_50_f_logsn)
  estim_q03_qr_50_m_Y[i] <- sum(beta_qr_q03_Y_50[[i]] * fix.X_50_m_logsn)
  estim_q05_qr_50_m_Y[i] <- sum(beta_qr_q05_Y_50[[i]] * fix.X_50_m_logsn)
  estim_q25_qr_50_m_Y[i] <- sum(beta_qr_q25_Y_50[[i]] * fix.X_50_m_logsn)
  estim_q50_qr_50_m_Y[i] <- sum(beta_qr_q50_Y_50[[i]] * fix.X_50_m_logsn)
  estim_q75_qr_50_m_Y[i] <- sum(beta_qr_q75_Y_50[[i]] * fix.X_50_m_logsn)
  estim_q95_qr_50_m_Y[i] <- sum(beta_qr_q95_Y_50[[i]] * fix.X_50_m_logsn)
  estim_q97_qr_50_m_Y[i] <- sum(beta_qr_q97_Y_50[[i]] * fix.X_50_m_logsn)
}

# Estimates of the quantiles (by gender) by using the 
# log-skew-normal linear regression model
estim_betas_logsn_50 <- cbind(estim_beta1_logsn_50,estim_beta2_logsn_50,estim_beta3_logsn_50)
estim_q03_logsn_50_f <- c()
estim_q05_logsn_50_f <- c()
estim_q25_logsn_50_f <- c()
estim_q50_logsn_50_f <- c()
estim_q75_logsn_50_f <- c()
estim_q95_logsn_50_f <- c()
estim_q97_logsn_50_f <- c()
estim_q03_logsn_50_m <- c()
estim_q05_logsn_50_m <- c()
estim_q25_logsn_50_m <- c()
estim_q50_logsn_50_m <- c()
estim_q75_logsn_50_m <- c()
estim_q95_logsn_50_m <- c()
estim_q97_logsn_50_m <- c()
for (i in 1:N) {
  estim_q03_logsn_50_f[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_f_logsn) + estim_omega_logsn_50[i]*estim_q03_stdsn_50[i])
  estim_q05_logsn_50_f[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_f_logsn) + estim_omega_logsn_50[i]*estim_q05_stdsn_50[i])
  estim_q25_logsn_50_f[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_f_logsn) + estim_omega_logsn_50[i]*estim_q25_stdsn_50[i])
  estim_q50_logsn_50_f[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_f_logsn) + estim_omega_logsn_50[i]*estim_q50_stdsn_50[i])
  estim_q75_logsn_50_f[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_f_logsn) + estim_omega_logsn_50[i]*estim_q75_stdsn_50[i])
  estim_q95_logsn_50_f[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_f_logsn) + estim_omega_logsn_50[i]*estim_q95_stdsn_50[i])
  estim_q97_logsn_50_f[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_f_logsn) + estim_omega_logsn_50[i]*estim_q97_stdsn_50[i])
  estim_q03_logsn_50_m[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_m_logsn) + estim_omega_logsn_50[i]*estim_q03_stdsn_50[i])
  estim_q05_logsn_50_m[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_m_logsn) + estim_omega_logsn_50[i]*estim_q05_stdsn_50[i])
  estim_q25_logsn_50_m[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_m_logsn) + estim_omega_logsn_50[i]*estim_q25_stdsn_50[i])
  estim_q50_logsn_50_m[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_m_logsn) + estim_omega_logsn_50[i]*estim_q50_stdsn_50[i])
  estim_q75_logsn_50_m[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_m_logsn) + estim_omega_logsn_50[i]*estim_q75_stdsn_50[i])
  estim_q95_logsn_50_m[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_m_logsn) + estim_omega_logsn_50[i]*estim_q95_stdsn_50[i])
  estim_q97_logsn_50_m[i] <- exp(sum(estim_betas_logsn_50[i,]*fix.X_50_m_logsn) + estim_omega_logsn_50[i]*estim_q97_stdsn_50[i])
}

# Medians
median_beta1_logsn_50 <- median(estim_beta1_logsn_50)
median_beta2_logsn_50 <- median(estim_beta2_logsn_50)
median_beta3_logsn_50 <- median(estim_beta3_logsn_50)
median_omega_logsn_50 <- median(estim_omega_logsn_50)
median_alpha_logsn_50 <- median(estim_alpha_logsn_50)

median_q03_logsn_50_f_Y <- median(estim_q03_logsn_50_f)
median_q05_logsn_50_f_Y <- median(estim_q05_logsn_50_f)
median_q25_logsn_50_f_Y <- median(estim_q25_logsn_50_f)
median_q50_logsn_50_f_Y <- median(estim_q50_logsn_50_f)  
median_q75_logsn_50_f_Y <- median(estim_q75_logsn_50_f)
median_q95_logsn_50_f_Y <- median(estim_q95_logsn_50_f)
median_q97_logsn_50_f_Y <- median(estim_q97_logsn_50_f)
median_q03_logsn_50_m_Y <- median(estim_q03_logsn_50_m)
median_q05_logsn_50_m_Y <- median(estim_q05_logsn_50_m)
median_q25_logsn_50_m_Y <- median(estim_q25_logsn_50_m)
median_q50_logsn_50_m_Y <- median(estim_q50_logsn_50_m)  
median_q75_logsn_50_m_Y <- median(estim_q75_logsn_50_m)
median_q95_logsn_50_m_Y <- median(estim_q95_logsn_50_m)
median_q97_logsn_50_m_Y <- median(estim_q97_logsn_50_m)

median_q03_qr_50_f_Y <- median(estim_q03_qr_50_f_Y)
median_q05_qr_50_f_Y <- median(estim_q05_qr_50_f_Y)
median_q25_qr_50_f_Y <- median(estim_q25_qr_50_f_Y)
median_q50_qr_50_f_Y <- median(estim_q50_qr_50_f_Y)
median_q75_qr_50_f_Y <- median(estim_q75_qr_50_f_Y)
median_q95_qr_50_f_Y <- median(estim_q95_qr_50_f_Y)
median_q97_qr_50_f_Y <- median(estim_q97_qr_50_f_Y)
median_q03_qr_50_m_Y <- median(estim_q03_qr_50_m_Y)
median_q05_qr_50_m_Y <- median(estim_q05_qr_50_m_Y)
median_q25_qr_50_m_Y <- median(estim_q25_qr_50_m_Y)
median_q50_qr_50_m_Y <- median(estim_q50_qr_50_m_Y)
median_q75_qr_50_m_Y <- median(estim_q75_qr_50_m_Y)
median_q95_qr_50_m_Y <- median(estim_q95_qr_50_m_Y)
median_q97_qr_50_m_Y <- median(estim_q97_qr_50_m_Y)

# Median Absolute Deviations (MAD)
mad_beta1_logsn_50 <- mad(estim_beta1_logsn_50, constant = 1)
mad_beta2_logsn_50 <- mad(estim_beta2_logsn_50, constant = 1)
mad_beta3_logsn_50 <- mad(estim_beta3_logsn_50, constant = 1)
mad_omega_logsn_50 <- mad(estim_omega_logsn_50, constant = 1)
mad_alpha_logsn_50 <- mad(estim_alpha_logsn_50, constant = 1)
mad_q03_logsn_50_f_Y <- mad(estim_q03_logsn_50_f, constant = 1)
mad_q05_logsn_50_f_Y <- mad(estim_q05_logsn_50_f, constant = 1)
mad_q25_logsn_50_f_Y <- mad(estim_q25_logsn_50_f, constant = 1)
mad_q50_logsn_50_f_Y <- mad(estim_q50_logsn_50_f, constant = 1)
mad_q75_logsn_50_f_Y <- mad(estim_q75_logsn_50_f, constant = 1)
mad_q95_logsn_50_f_Y <- mad(estim_q95_logsn_50_f, constant = 1)
mad_q97_logsn_50_f_Y <- mad(estim_q97_logsn_50_f, constant = 1)
mad_q03_logsn_50_m_Y <- mad(estim_q03_logsn_50_m, constant = 1)
mad_q05_logsn_50_m_Y <- mad(estim_q05_logsn_50_m, constant = 1)
mad_q25_logsn_50_m_Y <- mad(estim_q25_logsn_50_m, constant = 1)
mad_q50_logsn_50_m_Y <- mad(estim_q50_logsn_50_m, constant = 1)
mad_q75_logsn_50_m_Y <- mad(estim_q75_logsn_50_m, constant = 1)
mad_q95_logsn_50_m_Y <- mad(estim_q95_logsn_50_m, constant = 1)
mad_q97_logsn_50_m_Y <- mad(estim_q97_logsn_50_m, constant = 1)
mad_q03_qr_50_f_Y <- mad(estim_q03_qr_50_f_Y, constant = 1)
mad_q05_qr_50_f_Y <- mad(estim_q05_qr_50_f_Y, constant = 1)
mad_q25_qr_50_f_Y <- mad(estim_q25_qr_50_f_Y, constant = 1)
mad_q50_qr_50_f_Y <- mad(estim_q50_qr_50_f_Y, constant = 1)
mad_q75_qr_50_f_Y <- mad(estim_q75_qr_50_f_Y, constant = 1)
mad_q95_qr_50_f_Y <- mad(estim_q95_qr_50_f_Y, constant = 1)
mad_q97_qr_50_f_Y <- mad(estim_q97_qr_50_f_Y, constant = 1)
mad_q03_qr_50_m_Y <- mad(estim_q03_qr_50_m_Y, constant = 1)
mad_q05_qr_50_m_Y <- mad(estim_q05_qr_50_m_Y, constant = 1)
mad_q25_qr_50_m_Y <- mad(estim_q25_qr_50_m_Y, constant = 1)
mad_q50_qr_50_m_Y <- mad(estim_q50_qr_50_m_Y, constant = 1)
mad_q75_qr_50_m_Y <- mad(estim_q75_qr_50_m_Y, constant = 1)
mad_q95_qr_50_m_Y <- mad(estim_q95_qr_50_m_Y, constant = 1)
mad_q97_qr_50_m_Y <- mad(estim_q97_qr_50_m_Y, constant = 1)

#=======================
# Sample size n = 100
#=======================

# Estimates of the regression coefficients using the 
# classical QR
beta_qr_q03_Y_100 <- lapply(simulation_logsn_100, "[[", 7)
beta_qr_q05_Y_100 <- lapply(simulation_logsn_100, "[[", 8)
beta_qr_q25_Y_100 <- lapply(simulation_logsn_100, "[[", 9)
beta_qr_q50_Y_100 <- lapply(simulation_logsn_100, "[[", 10)
beta_qr_q75_Y_100 <- lapply(simulation_logsn_100, "[[", 11)
beta_qr_q95_Y_100 <- lapply(simulation_logsn_100, "[[", 12)
beta_qr_q97_Y_100 <- lapply(simulation_logsn_100, "[[", 13)

# Estimates of the 3%, 5%, 25%, 50%, 75%, 95% and 97% 
# quantiles (by gender) using linear quantile regression
estim_q03_qr_100_f_Y <- c()
estim_q05_qr_100_f_Y <- c()
estim_q25_qr_100_f_Y <- c()
estim_q50_qr_100_f_Y <- c()
estim_q75_qr_100_f_Y <- c()
estim_q95_qr_100_f_Y <- c()
estim_q97_qr_100_f_Y <- c()
estim_q03_qr_100_m_Y <- c()
estim_q05_qr_100_m_Y <- c()
estim_q25_qr_100_m_Y <- c()
estim_q50_qr_100_m_Y <- c()
estim_q75_qr_100_m_Y <- c()
estim_q95_qr_100_m_Y <- c()
estim_q97_qr_100_m_Y <- c()
for (i in 1:N) {
  estim_q03_qr_100_f_Y[i] <- sum(beta_qr_q03_Y_100[[i]] * fix.X_100_f_logsn)
  estim_q05_qr_100_f_Y[i] <- sum(beta_qr_q05_Y_100[[i]] * fix.X_100_f_logsn)
  estim_q25_qr_100_f_Y[i] <- sum(beta_qr_q25_Y_100[[i]] * fix.X_100_f_logsn)
  estim_q50_qr_100_f_Y[i] <- sum(beta_qr_q50_Y_100[[i]] * fix.X_100_f_logsn)
  estim_q75_qr_100_f_Y[i] <- sum(beta_qr_q75_Y_100[[i]] * fix.X_100_f_logsn)
  estim_q95_qr_100_f_Y[i] <- sum(beta_qr_q95_Y_100[[i]] * fix.X_100_f_logsn)
  estim_q97_qr_100_f_Y[i] <- sum(beta_qr_q97_Y_100[[i]] * fix.X_100_f_logsn)
  estim_q03_qr_100_m_Y[i] <- sum(beta_qr_q03_Y_100[[i]] * fix.X_100_m_logsn)
  estim_q05_qr_100_m_Y[i] <- sum(beta_qr_q05_Y_100[[i]] * fix.X_100_m_logsn)
  estim_q25_qr_100_m_Y[i] <- sum(beta_qr_q25_Y_100[[i]] * fix.X_100_m_logsn)
  estim_q50_qr_100_m_Y[i] <- sum(beta_qr_q50_Y_100[[i]] * fix.X_100_m_logsn)
  estim_q75_qr_100_m_Y[i] <- sum(beta_qr_q75_Y_100[[i]] * fix.X_100_m_logsn)
  estim_q95_qr_100_m_Y[i] <- sum(beta_qr_q95_Y_100[[i]] * fix.X_100_m_logsn)
  estim_q97_qr_100_m_Y[i] <- sum(beta_qr_q97_Y_100[[i]] * fix.X_100_m_logsn)
}

# Estimates of the quantiles (by gender) by using the 
# log-skew-normal linear regression model
estim_betas_logsn_100 <- cbind(estim_beta1_logsn_100,estim_beta2_logsn_100,estim_beta3_logsn_100)
estim_q03_logsn_100_f <- c()
estim_q05_logsn_100_f <- c()
estim_q25_logsn_100_f <- c()
estim_q50_logsn_100_f <- c()
estim_q75_logsn_100_f <- c()
estim_q95_logsn_100_f <- c()
estim_q97_logsn_100_f <- c()
estim_q03_logsn_100_m <- c()
estim_q05_logsn_100_m <- c()
estim_q25_logsn_100_m <- c()
estim_q50_logsn_100_m <- c()
estim_q75_logsn_100_m <- c()
estim_q95_logsn_100_m <- c()
estim_q97_logsn_100_m <- c()
for (i in 1:N) {
  estim_q03_logsn_100_f[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_f_logsn) + estim_omega_logsn_100[i]*estim_q03_stdsn_100[i])
  estim_q05_logsn_100_f[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_f_logsn) + estim_omega_logsn_100[i]*estim_q05_stdsn_100[i])
  estim_q25_logsn_100_f[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_f_logsn) + estim_omega_logsn_100[i]*estim_q25_stdsn_100[i])
  estim_q50_logsn_100_f[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_f_logsn) + estim_omega_logsn_100[i]*estim_q50_stdsn_100[i])
  estim_q75_logsn_100_f[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_f_logsn) + estim_omega_logsn_100[i]*estim_q75_stdsn_100[i])
  estim_q95_logsn_100_f[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_f_logsn) + estim_omega_logsn_100[i]*estim_q95_stdsn_100[i])
  estim_q97_logsn_100_f[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_f_logsn) + estim_omega_logsn_100[i]*estim_q97_stdsn_100[i])
  estim_q03_logsn_100_m[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_m_logsn) + estim_omega_logsn_100[i]*estim_q03_stdsn_100[i])
  estim_q05_logsn_100_m[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_m_logsn) + estim_omega_logsn_100[i]*estim_q05_stdsn_100[i])
  estim_q25_logsn_100_m[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_m_logsn) + estim_omega_logsn_100[i]*estim_q25_stdsn_100[i])
  estim_q50_logsn_100_m[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_m_logsn) + estim_omega_logsn_100[i]*estim_q50_stdsn_100[i])
  estim_q75_logsn_100_m[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_m_logsn) + estim_omega_logsn_100[i]*estim_q75_stdsn_100[i])
  estim_q95_logsn_100_m[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_m_logsn) + estim_omega_logsn_100[i]*estim_q95_stdsn_100[i])
  estim_q97_logsn_100_m[i] <- exp(sum(estim_betas_logsn_100[i,]*fix.X_100_m_logsn) + estim_omega_logsn_100[i]*estim_q97_stdsn_100[i])
}

# Medians
median_beta1_logsn_100 <- median(estim_beta1_logsn_100)
median_beta2_logsn_100 <- median(estim_beta2_logsn_100)
median_beta3_logsn_100 <- median(estim_beta3_logsn_100)
median_omega_logsn_100 <- median(estim_omega_logsn_100)
median_alpha_logsn_100 <- median(estim_alpha_logsn_100)

median_q03_logsn_100_f_Y <- median(estim_q03_logsn_100_f)
median_q05_logsn_100_f_Y <- median(estim_q05_logsn_100_f)
median_q25_logsn_100_f_Y <- median(estim_q25_logsn_100_f)
median_q50_logsn_100_f_Y <- median(estim_q50_logsn_100_f)  
median_q75_logsn_100_f_Y <- median(estim_q75_logsn_100_f)
median_q95_logsn_100_f_Y <- median(estim_q95_logsn_100_f)
median_q97_logsn_100_f_Y <- median(estim_q97_logsn_100_f)
median_q03_logsn_100_m_Y <- median(estim_q03_logsn_100_m)
median_q05_logsn_100_m_Y <- median(estim_q05_logsn_100_m)
median_q25_logsn_100_m_Y <- median(estim_q25_logsn_100_m)
median_q50_logsn_100_m_Y <- median(estim_q50_logsn_100_m)  
median_q75_logsn_100_m_Y <- median(estim_q75_logsn_100_m)
median_q95_logsn_100_m_Y <- median(estim_q95_logsn_100_m)
median_q97_logsn_100_m_Y <- median(estim_q97_logsn_100_m)

median_q03_qr_100_f_Y <- median(estim_q03_qr_100_f_Y)
median_q05_qr_100_f_Y <- median(estim_q05_qr_100_f_Y)
median_q25_qr_100_f_Y <- median(estim_q25_qr_100_f_Y)
median_q50_qr_100_f_Y <- median(estim_q50_qr_100_f_Y)
median_q75_qr_100_f_Y <- median(estim_q75_qr_100_f_Y)
median_q95_qr_100_f_Y <- median(estim_q95_qr_100_f_Y)
median_q97_qr_100_f_Y <- median(estim_q97_qr_100_f_Y)
median_q03_qr_100_m_Y <- median(estim_q03_qr_100_m_Y)
median_q05_qr_100_m_Y <- median(estim_q05_qr_100_m_Y)
median_q25_qr_100_m_Y <- median(estim_q25_qr_100_m_Y)
median_q50_qr_100_m_Y <- median(estim_q50_qr_100_m_Y)
median_q75_qr_100_m_Y <- median(estim_q75_qr_100_m_Y)
median_q95_qr_100_m_Y <- median(estim_q95_qr_100_m_Y)
median_q97_qr_100_m_Y <- median(estim_q97_qr_100_m_Y)

# Median Absolute Deviations (MAD)
mad_beta1_logsn_100 <- mad(estim_beta1_logsn_100, constant = 1)
mad_beta2_logsn_100 <- mad(estim_beta2_logsn_100, constant = 1)
mad_beta3_logsn_100 <- mad(estim_beta3_logsn_100, constant = 1)
mad_omega_logsn_100 <- mad(estim_omega_logsn_100, constant = 1)
mad_alpha_logsn_100 <- mad(estim_alpha_logsn_100, constant = 1)
mad_q03_logsn_100_f_Y <- mad(estim_q03_logsn_100_f, constant = 1)
mad_q05_logsn_100_f_Y <- mad(estim_q05_logsn_100_f, constant = 1)
mad_q25_logsn_100_f_Y <- mad(estim_q25_logsn_100_f, constant = 1)
mad_q50_logsn_100_f_Y <- mad(estim_q50_logsn_100_f, constant = 1)
mad_q75_logsn_100_f_Y <- mad(estim_q75_logsn_100_f, constant = 1)
mad_q95_logsn_100_f_Y <- mad(estim_q95_logsn_100_f, constant = 1)
mad_q97_logsn_100_f_Y <- mad(estim_q97_logsn_100_f, constant = 1)
mad_q03_logsn_100_m_Y <- mad(estim_q03_logsn_100_m, constant = 1)
mad_q05_logsn_100_m_Y <- mad(estim_q05_logsn_100_m, constant = 1)
mad_q25_logsn_100_m_Y <- mad(estim_q25_logsn_100_m, constant = 1)
mad_q50_logsn_100_m_Y <- mad(estim_q50_logsn_100_m, constant = 1)
mad_q75_logsn_100_m_Y <- mad(estim_q75_logsn_100_m, constant = 1)
mad_q95_logsn_100_m_Y <- mad(estim_q95_logsn_100_m, constant = 1)
mad_q97_logsn_100_m_Y <- mad(estim_q97_logsn_100_m, constant = 1)
mad_q03_qr_100_f_Y <- mad(estim_q03_qr_100_f_Y, constant = 1)
mad_q05_qr_100_f_Y <- mad(estim_q05_qr_100_f_Y, constant = 1)
mad_q25_qr_100_f_Y <- mad(estim_q25_qr_100_f_Y, constant = 1)
mad_q50_qr_100_f_Y <- mad(estim_q50_qr_100_f_Y, constant = 1)
mad_q75_qr_100_f_Y <- mad(estim_q75_qr_100_f_Y, constant = 1)
mad_q95_qr_100_f_Y <- mad(estim_q95_qr_100_f_Y, constant = 1)
mad_q97_qr_100_f_Y <- mad(estim_q97_qr_100_f_Y, constant = 1)
mad_q03_qr_100_m_Y <- mad(estim_q03_qr_100_m_Y, constant = 1)
mad_q05_qr_100_m_Y <- mad(estim_q05_qr_100_m_Y, constant = 1)
mad_q25_qr_100_m_Y <- mad(estim_q25_qr_100_m_Y, constant = 1)
mad_q50_qr_100_m_Y <- mad(estim_q50_qr_100_m_Y, constant = 1)
mad_q75_qr_100_m_Y <- mad(estim_q75_qr_100_m_Y, constant = 1)
mad_q95_qr_100_m_Y <- mad(estim_q95_qr_100_m_Y, constant = 1)
mad_q97_qr_100_m_Y <- mad(estim_q97_qr_100_m_Y, constant = 1)

#=======================
# Sample size n = 500
#=======================

# Estimates of the regression coefficients using the 
# classical QR
beta_qr_q03_Y_500 <- lapply(simulation_logsn_500, "[[", 7)
beta_qr_q05_Y_500 <- lapply(simulation_logsn_500, "[[", 8)
beta_qr_q25_Y_500 <- lapply(simulation_logsn_500, "[[", 9)
beta_qr_q50_Y_500 <- lapply(simulation_logsn_500, "[[", 10)
beta_qr_q75_Y_500 <- lapply(simulation_logsn_500, "[[", 11)
beta_qr_q95_Y_500 <- lapply(simulation_logsn_500, "[[", 12)
beta_qr_q97_Y_500 <- lapply(simulation_logsn_500, "[[", 13)

# Estimates of the 3%, 5%, 25%, 50%, 75%, 95% and 97% 
# quantiles (by gender) using linear quantile regression
estim_q03_qr_500_f_Y <- c()
estim_q05_qr_500_f_Y <- c()
estim_q25_qr_500_f_Y <- c()
estim_q50_qr_500_f_Y <- c()
estim_q75_qr_500_f_Y <- c()
estim_q95_qr_500_f_Y <- c()
estim_q97_qr_500_f_Y <- c()
estim_q03_qr_500_m_Y <- c()
estim_q05_qr_500_m_Y <- c()
estim_q25_qr_500_m_Y <- c()
estim_q50_qr_500_m_Y <- c()
estim_q75_qr_500_m_Y <- c()
estim_q95_qr_500_m_Y <- c()
estim_q97_qr_500_m_Y <- c()
for (i in 1:N) {
  estim_q03_qr_500_f_Y[i] <- sum(beta_qr_q03_Y_500[[i]] * fix.X_500_f_logsn)
  estim_q05_qr_500_f_Y[i] <- sum(beta_qr_q05_Y_500[[i]] * fix.X_500_f_logsn)
  estim_q25_qr_500_f_Y[i] <- sum(beta_qr_q25_Y_500[[i]] * fix.X_500_f_logsn)
  estim_q50_qr_500_f_Y[i] <- sum(beta_qr_q50_Y_500[[i]] * fix.X_500_f_logsn)
  estim_q75_qr_500_f_Y[i] <- sum(beta_qr_q75_Y_500[[i]] * fix.X_500_f_logsn)
  estim_q95_qr_500_f_Y[i] <- sum(beta_qr_q95_Y_500[[i]] * fix.X_500_f_logsn)
  estim_q97_qr_500_f_Y[i] <- sum(beta_qr_q97_Y_500[[i]] * fix.X_500_f_logsn)
  estim_q03_qr_500_m_Y[i] <- sum(beta_qr_q03_Y_500[[i]] * fix.X_500_m_logsn)
  estim_q05_qr_500_m_Y[i] <- sum(beta_qr_q05_Y_500[[i]] * fix.X_500_m_logsn)
  estim_q25_qr_500_m_Y[i] <- sum(beta_qr_q25_Y_500[[i]] * fix.X_500_m_logsn)
  estim_q50_qr_500_m_Y[i] <- sum(beta_qr_q50_Y_500[[i]] * fix.X_500_m_logsn)
  estim_q75_qr_500_m_Y[i] <- sum(beta_qr_q75_Y_500[[i]] * fix.X_500_m_logsn)
  estim_q95_qr_500_m_Y[i] <- sum(beta_qr_q95_Y_500[[i]] * fix.X_500_m_logsn)
  estim_q97_qr_500_m_Y[i] <- sum(beta_qr_q97_Y_500[[i]] * fix.X_500_m_logsn)
}

# Estimates of the quantiles (by gender) by using the 
# log-skew-normal linear regression model
estim_betas_logsn_500 <- cbind(estim_beta1_logsn_500,estim_beta2_logsn_500,estim_beta3_logsn_500)
estim_q03_logsn_500_f <- c()
estim_q05_logsn_500_f <- c()
estim_q25_logsn_500_f <- c()
estim_q50_logsn_500_f <- c()
estim_q75_logsn_500_f <- c()
estim_q95_logsn_500_f <- c()
estim_q97_logsn_500_f <- c()
estim_q03_logsn_500_m <- c()
estim_q05_logsn_500_m <- c()
estim_q25_logsn_500_m <- c()
estim_q50_logsn_500_m <- c()
estim_q75_logsn_500_m <- c()
estim_q95_logsn_500_m <- c()
estim_q97_logsn_500_m <- c()
for (i in 1:N) {
  estim_q03_logsn_500_f[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_f_logsn) + estim_omega_logsn_500[i]*estim_q03_stdsn_500[i])
  estim_q05_logsn_500_f[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_f_logsn) + estim_omega_logsn_500[i]*estim_q05_stdsn_500[i])
  estim_q25_logsn_500_f[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_f_logsn) + estim_omega_logsn_500[i]*estim_q25_stdsn_500[i])
  estim_q50_logsn_500_f[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_f_logsn) + estim_omega_logsn_500[i]*estim_q50_stdsn_500[i])
  estim_q75_logsn_500_f[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_f_logsn) + estim_omega_logsn_500[i]*estim_q75_stdsn_500[i])
  estim_q95_logsn_500_f[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_f_logsn) + estim_omega_logsn_500[i]*estim_q95_stdsn_500[i])
  estim_q97_logsn_500_f[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_f_logsn) + estim_omega_logsn_500[i]*estim_q97_stdsn_500[i])
  estim_q03_logsn_500_m[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_m_logsn) + estim_omega_logsn_500[i]*estim_q03_stdsn_500[i])
  estim_q05_logsn_500_m[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_m_logsn) + estim_omega_logsn_500[i]*estim_q05_stdsn_500[i])
  estim_q25_logsn_500_m[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_m_logsn) + estim_omega_logsn_500[i]*estim_q25_stdsn_500[i])
  estim_q50_logsn_500_m[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_m_logsn) + estim_omega_logsn_500[i]*estim_q50_stdsn_500[i])
  estim_q75_logsn_500_m[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_m_logsn) + estim_omega_logsn_500[i]*estim_q75_stdsn_500[i])
  estim_q95_logsn_500_m[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_m_logsn) + estim_omega_logsn_500[i]*estim_q95_stdsn_500[i])
  estim_q97_logsn_500_m[i] <- exp(sum(estim_betas_logsn_500[i,]*fix.X_500_m_logsn) + estim_omega_logsn_500[i]*estim_q97_stdsn_500[i])
}

# Medians
median_beta1_logsn_500 <- median(estim_beta1_logsn_500)
median_beta2_logsn_500 <- median(estim_beta2_logsn_500)
median_beta3_logsn_500 <- median(estim_beta3_logsn_500)
median_omega_logsn_500 <- median(estim_omega_logsn_500)
median_alpha_logsn_500 <- median(estim_alpha_logsn_500)

median_q03_logsn_500_f_Y <- median(estim_q03_logsn_500_f)
median_q05_logsn_500_f_Y <- median(estim_q05_logsn_500_f)
median_q25_logsn_500_f_Y <- median(estim_q25_logsn_500_f)
median_q50_logsn_500_f_Y <- median(estim_q50_logsn_500_f)  
median_q75_logsn_500_f_Y <- median(estim_q75_logsn_500_f)
median_q95_logsn_500_f_Y <- median(estim_q95_logsn_500_f)
median_q97_logsn_500_f_Y <- median(estim_q97_logsn_500_f)
median_q03_logsn_500_m_Y <- median(estim_q03_logsn_500_m)
median_q05_logsn_500_m_Y <- median(estim_q05_logsn_500_m)
median_q25_logsn_500_m_Y <- median(estim_q25_logsn_500_m)
median_q50_logsn_500_m_Y <- median(estim_q50_logsn_500_m)  
median_q75_logsn_500_m_Y <- median(estim_q75_logsn_500_m)
median_q95_logsn_500_m_Y <- median(estim_q95_logsn_500_m)
median_q97_logsn_500_m_Y <- median(estim_q97_logsn_500_m)

median_q03_qr_500_f_Y <- median(estim_q03_qr_500_f_Y)
median_q05_qr_500_f_Y <- median(estim_q05_qr_500_f_Y)
median_q25_qr_500_f_Y <- median(estim_q25_qr_500_f_Y)
median_q50_qr_500_f_Y <- median(estim_q50_qr_500_f_Y)
median_q75_qr_500_f_Y <- median(estim_q75_qr_500_f_Y)
median_q95_qr_500_f_Y <- median(estim_q95_qr_500_f_Y)
median_q97_qr_500_f_Y <- median(estim_q97_qr_500_f_Y)
median_q03_qr_500_m_Y <- median(estim_q03_qr_500_m_Y)
median_q05_qr_500_m_Y <- median(estim_q05_qr_500_m_Y)
median_q25_qr_500_m_Y <- median(estim_q25_qr_500_m_Y)
median_q50_qr_500_m_Y <- median(estim_q50_qr_500_m_Y)
median_q75_qr_500_m_Y <- median(estim_q75_qr_500_m_Y)
median_q95_qr_500_m_Y <- median(estim_q95_qr_500_m_Y)
median_q97_qr_500_m_Y <- median(estim_q97_qr_500_m_Y)

# Median Absolute Deviations (MAD)
mad_beta1_logsn_500 <- mad(estim_beta1_logsn_500, constant = 1)
mad_beta2_logsn_500 <- mad(estim_beta2_logsn_500, constant = 1)
mad_beta3_logsn_500 <- mad(estim_beta3_logsn_500, constant = 1)
mad_omega_logsn_500 <- mad(estim_omega_logsn_500, constant = 1)
mad_alpha_logsn_500 <- mad(estim_alpha_logsn_500, constant = 1)
mad_q03_logsn_500_f_Y <- mad(estim_q03_logsn_500_f, constant = 1)
mad_q05_logsn_500_f_Y <- mad(estim_q05_logsn_500_f, constant = 1)
mad_q25_logsn_500_f_Y <- mad(estim_q25_logsn_500_f, constant = 1)
mad_q50_logsn_500_f_Y <- mad(estim_q50_logsn_500_f, constant = 1)
mad_q75_logsn_500_f_Y <- mad(estim_q75_logsn_500_f, constant = 1)
mad_q95_logsn_500_f_Y <- mad(estim_q95_logsn_500_f, constant = 1)
mad_q97_logsn_500_f_Y <- mad(estim_q97_logsn_500_f, constant = 1)
mad_q03_logsn_500_m_Y <- mad(estim_q03_logsn_500_m, constant = 1)
mad_q05_logsn_500_m_Y <- mad(estim_q05_logsn_500_m, constant = 1)
mad_q25_logsn_500_m_Y <- mad(estim_q25_logsn_500_m, constant = 1)
mad_q50_logsn_500_m_Y <- mad(estim_q50_logsn_500_m, constant = 1)
mad_q75_logsn_500_m_Y <- mad(estim_q75_logsn_500_m, constant = 1)
mad_q95_logsn_500_m_Y <- mad(estim_q95_logsn_500_m, constant = 1)
mad_q97_logsn_500_m_Y <- mad(estim_q97_logsn_500_m, constant = 1)
mad_q03_qr_500_f_Y <- mad(estim_q03_qr_500_f_Y, constant = 1)
mad_q05_qr_500_f_Y <- mad(estim_q05_qr_500_f_Y, constant = 1)
mad_q25_qr_500_f_Y <- mad(estim_q25_qr_500_f_Y, constant = 1)
mad_q50_qr_500_f_Y <- mad(estim_q50_qr_500_f_Y, constant = 1)
mad_q75_qr_500_f_Y <- mad(estim_q75_qr_500_f_Y, constant = 1)
mad_q95_qr_500_f_Y <- mad(estim_q95_qr_500_f_Y, constant = 1)
mad_q97_qr_500_f_Y <- mad(estim_q97_qr_500_f_Y, constant = 1)
mad_q03_qr_500_m_Y <- mad(estim_q03_qr_500_m_Y, constant = 1)
mad_q05_qr_500_m_Y <- mad(estim_q05_qr_500_m_Y, constant = 1)
mad_q25_qr_500_m_Y <- mad(estim_q25_qr_500_m_Y, constant = 1)
mad_q50_qr_500_m_Y <- mad(estim_q50_qr_500_m_Y, constant = 1)
mad_q75_qr_500_m_Y <- mad(estim_q75_qr_500_m_Y, constant = 1)
mad_q95_qr_500_m_Y <- mad(estim_q95_qr_500_m_Y, constant = 1)
mad_q97_qr_500_m_Y <- mad(estim_q97_qr_500_m_Y, constant = 1)

#=======================
# Sample size n = 1000
#=======================

# Estimates of the regression coefficients using the 
# classical QR
beta_qr_q03_Y_1000 <- lapply(simulation_logsn_1000, "[[", 7)
beta_qr_q05_Y_1000 <- lapply(simulation_logsn_1000, "[[", 8)
beta_qr_q25_Y_1000 <- lapply(simulation_logsn_1000, "[[", 9)
beta_qr_q50_Y_1000 <- lapply(simulation_logsn_1000, "[[", 10)
beta_qr_q75_Y_1000 <- lapply(simulation_logsn_1000, "[[", 11)
beta_qr_q95_Y_1000 <- lapply(simulation_logsn_1000, "[[", 12)
beta_qr_q97_Y_1000 <- lapply(simulation_logsn_1000, "[[", 13)

# Estimates of the 3%, 5%, 25%, 50%, 75%, 95% and 97% 
# quantiles (by gender) using linear quantile regression
estim_q03_qr_1000_f_Y <- c()
estim_q05_qr_1000_f_Y <- c()
estim_q25_qr_1000_f_Y <- c()
estim_q50_qr_1000_f_Y <- c()
estim_q75_qr_1000_f_Y <- c()
estim_q95_qr_1000_f_Y <- c()
estim_q97_qr_1000_f_Y <- c()
estim_q03_qr_1000_m_Y <- c()
estim_q05_qr_1000_m_Y <- c()
estim_q25_qr_1000_m_Y <- c()
estim_q50_qr_1000_m_Y <- c()
estim_q75_qr_1000_m_Y <- c()
estim_q95_qr_1000_m_Y <- c()
estim_q97_qr_1000_m_Y <- c()
for (i in 1:N) {
  estim_q03_qr_1000_f_Y[i] <- sum(beta_qr_q03_Y_1000[[i]] * fix.X_1000_f_logsn)
  estim_q05_qr_1000_f_Y[i] <- sum(beta_qr_q05_Y_1000[[i]] * fix.X_1000_f_logsn)
  estim_q25_qr_1000_f_Y[i] <- sum(beta_qr_q25_Y_1000[[i]] * fix.X_1000_f_logsn)
  estim_q50_qr_1000_f_Y[i] <- sum(beta_qr_q50_Y_1000[[i]] * fix.X_1000_f_logsn)
  estim_q75_qr_1000_f_Y[i] <- sum(beta_qr_q75_Y_1000[[i]] * fix.X_1000_f_logsn)
  estim_q95_qr_1000_f_Y[i] <- sum(beta_qr_q95_Y_1000[[i]] * fix.X_1000_f_logsn)
  estim_q97_qr_1000_f_Y[i] <- sum(beta_qr_q97_Y_1000[[i]] * fix.X_1000_f_logsn)
  estim_q03_qr_1000_m_Y[i] <- sum(beta_qr_q03_Y_1000[[i]] * fix.X_1000_m_logsn)
  estim_q05_qr_1000_m_Y[i] <- sum(beta_qr_q05_Y_1000[[i]] * fix.X_1000_m_logsn)
  estim_q25_qr_1000_m_Y[i] <- sum(beta_qr_q25_Y_1000[[i]] * fix.X_1000_m_logsn)
  estim_q50_qr_1000_m_Y[i] <- sum(beta_qr_q50_Y_1000[[i]] * fix.X_1000_m_logsn)
  estim_q75_qr_1000_m_Y[i] <- sum(beta_qr_q75_Y_1000[[i]] * fix.X_1000_m_logsn)
  estim_q95_qr_1000_m_Y[i] <- sum(beta_qr_q95_Y_1000[[i]] * fix.X_1000_m_logsn)
  estim_q97_qr_1000_m_Y[i] <- sum(beta_qr_q97_Y_1000[[i]] * fix.X_1000_m_logsn)
}

# Estimates of the quantiles (by gender) by using the 
# log-skew-normal linear regression model
estim_betas_logsn_1000 <- cbind(estim_beta1_logsn_1000,estim_beta2_logsn_1000,estim_beta3_logsn_1000)
estim_q03_logsn_1000_f <- c()
estim_q05_logsn_1000_f <- c()
estim_q25_logsn_1000_f <- c()
estim_q50_logsn_1000_f <- c()
estim_q75_logsn_1000_f <- c()
estim_q95_logsn_1000_f <- c()
estim_q97_logsn_1000_f <- c()
estim_q03_logsn_1000_m <- c()
estim_q05_logsn_1000_m <- c()
estim_q25_logsn_1000_m <- c()
estim_q50_logsn_1000_m <- c()
estim_q75_logsn_1000_m <- c()
estim_q95_logsn_1000_m <- c()
estim_q97_logsn_1000_m <- c()
for (i in 1:N) {
  estim_q03_logsn_1000_f[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_f_logsn) + estim_omega_logsn_1000[i]*estim_q03_stdsn_1000[i])
  estim_q05_logsn_1000_f[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_f_logsn) + estim_omega_logsn_1000[i]*estim_q05_stdsn_1000[i])
  estim_q25_logsn_1000_f[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_f_logsn) + estim_omega_logsn_1000[i]*estim_q25_stdsn_1000[i])
  estim_q50_logsn_1000_f[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_f_logsn) + estim_omega_logsn_1000[i]*estim_q50_stdsn_1000[i])
  estim_q75_logsn_1000_f[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_f_logsn) + estim_omega_logsn_1000[i]*estim_q75_stdsn_1000[i])
  estim_q95_logsn_1000_f[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_f_logsn) + estim_omega_logsn_1000[i]*estim_q95_stdsn_1000[i])
  estim_q97_logsn_1000_f[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_f_logsn) + estim_omega_logsn_1000[i]*estim_q97_stdsn_1000[i])
  estim_q03_logsn_1000_m[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_m_logsn) + estim_omega_logsn_1000[i]*estim_q03_stdsn_1000[i])
  estim_q05_logsn_1000_m[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_m_logsn) + estim_omega_logsn_1000[i]*estim_q05_stdsn_1000[i])
  estim_q25_logsn_1000_m[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_m_logsn) + estim_omega_logsn_1000[i]*estim_q25_stdsn_1000[i])
  estim_q50_logsn_1000_m[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_m_logsn) + estim_omega_logsn_1000[i]*estim_q50_stdsn_1000[i])
  estim_q75_logsn_1000_m[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_m_logsn) + estim_omega_logsn_1000[i]*estim_q75_stdsn_1000[i])
  estim_q95_logsn_1000_m[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_m_logsn) + estim_omega_logsn_1000[i]*estim_q95_stdsn_1000[i])
  estim_q97_logsn_1000_m[i] <- exp(sum(estim_betas_logsn_1000[i,]*fix.X_1000_m_logsn) + estim_omega_logsn_1000[i]*estim_q97_stdsn_1000[i])
}

# Medians
median_beta1_logsn_1000 <- median(estim_beta1_logsn_1000)
median_beta2_logsn_1000 <- median(estim_beta2_logsn_1000)
median_beta3_logsn_1000 <- median(estim_beta3_logsn_1000)
median_omega_logsn_1000 <- median(estim_omega_logsn_1000)
median_alpha_logsn_1000 <- median(estim_alpha_logsn_1000)

median_q03_logsn_1000_f_Y <- median(estim_q03_logsn_1000_f)
median_q05_logsn_1000_f_Y <- median(estim_q05_logsn_1000_f)
median_q25_logsn_1000_f_Y <- median(estim_q25_logsn_1000_f)
median_q50_logsn_1000_f_Y <- median(estim_q50_logsn_1000_f)  
median_q75_logsn_1000_f_Y <- median(estim_q75_logsn_1000_f)
median_q95_logsn_1000_f_Y <- median(estim_q95_logsn_1000_f)
median_q97_logsn_1000_f_Y <- median(estim_q97_logsn_1000_f)
median_q03_logsn_1000_m_Y <- median(estim_q03_logsn_1000_m)
median_q05_logsn_1000_m_Y <- median(estim_q05_logsn_1000_m)
median_q25_logsn_1000_m_Y <- median(estim_q25_logsn_1000_m)
median_q50_logsn_1000_m_Y <- median(estim_q50_logsn_1000_m)  
median_q75_logsn_1000_m_Y <- median(estim_q75_logsn_1000_m)
median_q95_logsn_1000_m_Y <- median(estim_q95_logsn_1000_m)
median_q97_logsn_1000_m_Y <- median(estim_q97_logsn_1000_m)

median_q03_qr_1000_f_Y <- median(estim_q03_qr_1000_f_Y)
median_q05_qr_1000_f_Y <- median(estim_q05_qr_1000_f_Y)
median_q25_qr_1000_f_Y <- median(estim_q25_qr_1000_f_Y)
median_q50_qr_1000_f_Y <- median(estim_q50_qr_1000_f_Y)
median_q75_qr_1000_f_Y <- median(estim_q75_qr_1000_f_Y)
median_q95_qr_1000_f_Y <- median(estim_q95_qr_1000_f_Y)
median_q97_qr_1000_f_Y <- median(estim_q97_qr_1000_f_Y)
median_q03_qr_1000_m_Y <- median(estim_q03_qr_1000_m_Y)
median_q05_qr_1000_m_Y <- median(estim_q05_qr_1000_m_Y)
median_q25_qr_1000_m_Y <- median(estim_q25_qr_1000_m_Y)
median_q50_qr_1000_m_Y <- median(estim_q50_qr_1000_m_Y)
median_q75_qr_1000_m_Y <- median(estim_q75_qr_1000_m_Y)
median_q95_qr_1000_m_Y <- median(estim_q95_qr_1000_m_Y)
median_q97_qr_1000_m_Y <- median(estim_q97_qr_1000_m_Y)

# Median Absolute Deviations (MAD)
mad_beta1_logsn_1000 <- mad(estim_beta1_logsn_1000, constant = 1)
mad_beta2_logsn_1000 <- mad(estim_beta2_logsn_1000, constant = 1)
mad_beta3_logsn_1000 <- mad(estim_beta3_logsn_1000, constant = 1)
mad_omega_logsn_1000 <- mad(estim_omega_logsn_1000, constant = 1)
mad_alpha_logsn_1000 <- mad(estim_alpha_logsn_1000, constant = 1)
mad_q03_logsn_1000_f_Y <- mad(estim_q03_logsn_1000_f, constant = 1)
mad_q05_logsn_1000_f_Y <- mad(estim_q05_logsn_1000_f, constant = 1)
mad_q25_logsn_1000_f_Y <- mad(estim_q25_logsn_1000_f, constant = 1)
mad_q50_logsn_1000_f_Y <- mad(estim_q50_logsn_1000_f, constant = 1)
mad_q75_logsn_1000_f_Y <- mad(estim_q75_logsn_1000_f, constant = 1)
mad_q95_logsn_1000_f_Y <- mad(estim_q95_logsn_1000_f, constant = 1)
mad_q97_logsn_1000_f_Y <- mad(estim_q97_logsn_1000_f, constant = 1)
mad_q03_logsn_1000_m_Y <- mad(estim_q03_logsn_1000_m, constant = 1)
mad_q05_logsn_1000_m_Y <- mad(estim_q05_logsn_1000_m, constant = 1)
mad_q25_logsn_1000_m_Y <- mad(estim_q25_logsn_1000_m, constant = 1)
mad_q50_logsn_1000_m_Y <- mad(estim_q50_logsn_1000_m, constant = 1)
mad_q75_logsn_1000_m_Y <- mad(estim_q75_logsn_1000_m, constant = 1)
mad_q95_logsn_1000_m_Y <- mad(estim_q95_logsn_1000_m, constant = 1)
mad_q97_logsn_1000_m_Y <- mad(estim_q97_logsn_1000_m, constant = 1)
mad_q03_qr_1000_f_Y <- mad(estim_q03_qr_1000_f_Y, constant = 1)
mad_q05_qr_1000_f_Y <- mad(estim_q05_qr_1000_f_Y, constant = 1)
mad_q25_qr_1000_f_Y <- mad(estim_q25_qr_1000_f_Y, constant = 1)
mad_q50_qr_1000_f_Y <- mad(estim_q50_qr_1000_f_Y, constant = 1)
mad_q75_qr_1000_f_Y <- mad(estim_q75_qr_1000_f_Y, constant = 1)
mad_q95_qr_1000_f_Y <- mad(estim_q95_qr_1000_f_Y, constant = 1)
mad_q97_qr_1000_f_Y <- mad(estim_q97_qr_1000_f_Y, constant = 1)
mad_q03_qr_1000_m_Y <- mad(estim_q03_qr_1000_m_Y, constant = 1)
mad_q05_qr_1000_m_Y <- mad(estim_q05_qr_1000_m_Y, constant = 1)
mad_q25_qr_1000_m_Y <- mad(estim_q25_qr_1000_m_Y, constant = 1)
mad_q50_qr_1000_m_Y <- mad(estim_q50_qr_1000_m_Y, constant = 1)
mad_q75_qr_1000_m_Y <- mad(estim_q75_qr_1000_m_Y, constant = 1)
mad_q95_qr_1000_m_Y <- mad(estim_q95_qr_1000_m_Y, constant = 1)
mad_q97_qr_1000_m_Y <- mad(estim_q97_qr_1000_m_Y, constant = 1)

# Tables of simulation studies

# Table 1 (median and MAD of parameters obtained from simulation 
# study of LSNLRM)
options(scipen = 999)
names_param <- c("beta_1", "beta_2", "beta_3", "omega", "alpha")
true_param_logsn <- c(true_beta1_logsn,
                      true_beta2_logsn, 
                      true_beta3_logsn, 
                      true_omega_logsn, 
                      true_alpha_logsn)
median_param_logsn_50 <- c(median_beta1_logsn_50, 
                           median_beta2_logsn_50, 
                           median_beta3_logsn_50, 
                           median_omega_logsn_50, 
                           median_alpha_logsn_50)
mad_param_logsn_50 <- c(mad_beta1_logsn_50, 
                        mad_beta2_logsn_50, 
                        mad_beta3_logsn_50, 
                        mad_omega_logsn_50, 
                        mad_alpha_logsn_50)
table_param_logsn_50 <- round(data.frame(true_param_logsn, 
                                         median_param_logsn_50, 
                                         mad_param_logsn_50, 
                                         row.names = names_param),4)
names(table_param_logsn_50) <- c("True", "Median", "MAD")

options(scipen = 999)
names_param <- c("beta_1", "beta_2", "beta_3", "omega", "alpha")
true_param_logsn <- c(true_beta1_logsn,
                      true_beta2_logsn, 
                      true_beta3_logsn, 
                      true_omega_logsn, 
                      true_alpha_logsn)
median_param_logsn_100 <- c(median_beta1_logsn_100, 
                            median_beta2_logsn_100, 
                            median_beta3_logsn_100, 
                            median_omega_logsn_100, 
                            median_alpha_logsn_100)
mad_param_logsn_100 <- c(mad_beta1_logsn_100, 
                         mad_beta2_logsn_100, 
                         mad_beta3_logsn_100, 
                         mad_omega_logsn_100, 
                         mad_alpha_logsn_100)
table_param_logsn_100 <- round(data.frame(true_param_logsn, 
                                          median_param_logsn_100, 
                                          mad_param_logsn_100, 
                                          row.names = names_param),4)
names(table_param_logsn_100) <- c("True", "Median", "MAD")

options(scipen = 999)
names_param <- c("beta_1", "beta_2", "beta_3", "omega", "alpha")
true_param_logsn <- c(true_beta1_logsn,
                      true_beta2_logsn, 
                      true_beta3_logsn, 
                      true_omega_logsn, 
                      true_alpha_logsn)
median_param_logsn_500 <- c(median_beta1_logsn_500, 
                            median_beta2_logsn_500, 
                            median_beta3_logsn_500, 
                            median_omega_logsn_500, 
                            median_alpha_logsn_500)
mad_param_logsn_500 <- c(mad_beta1_logsn_500, 
                         mad_beta2_logsn_500, 
                         mad_beta3_logsn_500, 
                         mad_omega_logsn_500, 
                         mad_alpha_logsn_500)
table_param_logsn_500 <- round(data.frame(true_param_logsn, 
                                          median_param_logsn_500, 
                                          mad_param_logsn_500, 
                                          row.names = names_param),4)
names(table_param_logsn_500) <- c("True", "Median", "MAD")

options(scipen = 999)
names_param <- c("beta_1", "beta_2", "beta_3", "omega", "alpha")
true_param_logsn <- c(true_beta1_logsn,
                      true_beta2_logsn, 
                      true_beta3_logsn, 
                      true_omega_logsn, 
                      true_alpha_logsn)
median_param_logsn_1000 <- c(median_beta1_logsn_1000, 
                             median_beta2_logsn_1000, 
                             median_beta3_logsn_1000, 
                             median_omega_logsn_1000, 
                             median_alpha_logsn_1000)
mad_param_logsn_1000 <- c(mad_beta1_logsn_1000, 
                          mad_beta2_logsn_1000, 
                          mad_beta3_logsn_1000, 
                          mad_omega_logsn_1000, 
                          mad_alpha_logsn_1000)
table_param_logsn_1000 <- round(data.frame(true_param_logsn, 
                                           median_param_logsn_1000, 
                                           mad_param_logsn_1000, 
                                           row.names = names_param),4)
names(table_param_logsn_1000) <- c("True", "Median", "MAD")

table_sim_study_param_logsn <- list(table_param_logsn_50, 
                                    table_param_logsn_100, 
                                    table_param_logsn_500, 
                                    table_param_logsn_1000)

names(table_sim_study_param_logsn) <- c("n=50", "n=100", "n=500", "n=1000")

# Table 2 (median and MAD of 3%, 5%, 25%, 50%, 75%, 95% and 97% quantiles 
# obtained from simulation study of LSNLRM and classical QR)
options(scipen = 999)
name_quantiles_f <- paste("Y_", c(t(outer(c("0.03", "0.05", "0.25", "0.50", "0.75", "0.95", "0.97"), "f", FUN = "paste", sep = "_"))), sep = "")
name_quantiles_m <- paste("Y_", c(t(outer(c("0.03", "0.05", "0.25", "0.50", "0.75", "0.95", "0.97"), "m", FUN = "paste", sep = "_"))), sep = "")
name_quantiles <- c(name_quantiles_f, name_quantiles_m)
true_quantiles_logsn <- c(true_q03_f_logsn, 
                          true_q05_f_logsn, 
                          true_q25_f_logsn, 
                          true_q50_f_logsn, 
                          true_q75_f_logsn, 
                          true_q95_f_logsn, 
                          true_q97_f_logsn,
                          true_q03_m_logsn, 
                          true_q05_m_logsn, 
                          true_q25_m_logsn, 
                          true_q50_m_logsn, 
                          true_q75_m_logsn, 
                          true_q95_m_logsn, 
                          true_q97_m_logsn)

median_quantiles_logsn_50 <- c(median_q03_logsn_50_f_Y, 
                               median_q05_logsn_50_f_Y, 
                               median_q25_logsn_50_f_Y, 
                               median_q50_logsn_50_f_Y, 
                               median_q75_logsn_50_f_Y, 
                               median_q95_logsn_50_f_Y, 
                               median_q97_logsn_50_f_Y, 
                               median_q03_logsn_50_m_Y, 
                               median_q05_logsn_50_m_Y, 
                               median_q25_logsn_50_m_Y, 
                               median_q50_logsn_50_m_Y, 
                               median_q75_logsn_50_m_Y, 
                               median_q95_logsn_50_m_Y, 
                               median_q97_logsn_50_m_Y)
median_quantiles_qr_50 <- c(median_q03_qr_50_f_Y, 
                            median_q05_qr_50_f_Y, 
                            median_q25_qr_50_f_Y, 
                            median_q50_qr_50_f_Y, 
                            median_q75_qr_50_f_Y, 
                            median_q95_qr_50_f_Y, 
                            median_q97_qr_50_f_Y, 
                            median_q03_qr_50_m_Y, 
                            median_q05_qr_50_m_Y, 
                            median_q25_qr_50_m_Y, 
                            median_q50_qr_50_m_Y, 
                            median_q75_qr_50_m_Y, 
                            median_q95_qr_50_m_Y, 
                            median_q97_qr_50_m_Y)
mad_quantiles_logsn_50 <- c(mad_q03_logsn_50_f_Y, 
                            mad_q05_logsn_50_f_Y, 
                            mad_q25_logsn_50_f_Y, 
                            mad_q50_logsn_50_f_Y, 
                            mad_q75_logsn_50_f_Y, 
                            mad_q95_logsn_50_f_Y, 
                            mad_q97_logsn_50_f_Y, 
                            mad_q03_logsn_50_m_Y, 
                            mad_q05_logsn_50_m_Y, 
                            mad_q25_logsn_50_m_Y, 
                            mad_q50_logsn_50_m_Y, 
                            mad_q75_logsn_50_m_Y, 
                            mad_q95_logsn_50_m_Y, 
                            mad_q97_logsn_50_m_Y)
mad_quantiles_qr_50 <- c(mad_q03_qr_50_f_Y, 
                         mad_q05_qr_50_f_Y, 
                         mad_q25_qr_50_f_Y, 
                         mad_q50_qr_50_f_Y, 
                         mad_q75_qr_50_f_Y, 
                         mad_q95_qr_50_f_Y, 
                         mad_q97_qr_50_f_Y, 
                         mad_q03_qr_50_m_Y, 
                         mad_q05_qr_50_m_Y, 
                         mad_q25_qr_50_m_Y, 
                         mad_q50_qr_50_m_Y, 
                         mad_q75_qr_50_m_Y, 
                         mad_q95_qr_50_m_Y, 
                         mad_q97_qr_50_m_Y)
table_quant_logsn_qr_50 <- data.frame(true_quantiles_logsn, 
                                      median_quantiles_logsn_50, 
                                      mad_quantiles_logsn_50, 
                                      median_quantiles_qr_50, 
                                      mad_quantiles_qr_50,
                                      row.names = name_quantiles)
table_quant_logsn_qr_50 <- signif(round(table_quant_logsn_qr_50, 4), 5)
names(table_quant_logsn_qr_50) <- c("True", "Median_logsn", "MAD_logsn", "Median_qr", "MAD_qr")

options(scipen = 999)
name_quantiles_f <- paste("Y_", c(t(outer(c("0.03", "0.05", "0.25", "0.50", "0.75", "0.95", "0.97"), "f", FUN = "paste", sep = "_"))), sep = "")
name_quantiles_m <- paste("Y_", c(t(outer(c("0.03", "0.05", "0.25", "0.50", "0.75", "0.95", "0.97"), "m", FUN = "paste", sep = "_"))), sep = "")
name_quantiles <- c(name_quantiles_f, name_quantiles_m)
true_quantiles_logsn <- c(true_q03_f_logsn, 
                          true_q05_f_logsn, 
                          true_q25_f_logsn, 
                          true_q50_f_logsn, 
                          true_q75_f_logsn, 
                          true_q95_f_logsn, 
                          true_q97_f_logsn,
                          true_q03_m_logsn, 
                          true_q05_m_logsn, 
                          true_q25_m_logsn, 
                          true_q50_m_logsn, 
                          true_q75_m_logsn, 
                          true_q95_m_logsn, 
                          true_q97_m_logsn)

median_quantiles_logsn_100 <- c(median_q03_logsn_100_f_Y, 
                                median_q05_logsn_100_f_Y, 
                                median_q25_logsn_100_f_Y, 
                                median_q50_logsn_100_f_Y, 
                                median_q75_logsn_100_f_Y, 
                                median_q95_logsn_100_f_Y, 
                                median_q97_logsn_100_f_Y, 
                                median_q03_logsn_100_m_Y, 
                                median_q05_logsn_100_m_Y, 
                                median_q25_logsn_100_m_Y, 
                                median_q50_logsn_100_m_Y, 
                                median_q75_logsn_100_m_Y, 
                                median_q95_logsn_100_m_Y, 
                                median_q97_logsn_100_m_Y)
median_quantiles_qr_100 <- c(median_q03_qr_100_f_Y, 
                             median_q05_qr_100_f_Y, 
                             median_q25_qr_100_f_Y, 
                             median_q50_qr_100_f_Y, 
                             median_q75_qr_100_f_Y, 
                             median_q95_qr_100_f_Y, 
                             median_q97_qr_100_f_Y, 
                             median_q03_qr_100_m_Y, 
                             median_q05_qr_100_m_Y, 
                             median_q25_qr_100_m_Y, 
                             median_q50_qr_100_m_Y, 
                             median_q75_qr_100_m_Y, 
                             median_q95_qr_100_m_Y, 
                             median_q97_qr_100_m_Y)
mad_quantiles_logsn_100 <- c(mad_q03_logsn_100_f_Y, 
                             mad_q05_logsn_100_f_Y, 
                             mad_q25_logsn_100_f_Y, 
                             mad_q50_logsn_100_f_Y, 
                             mad_q75_logsn_100_f_Y, 
                             mad_q95_logsn_100_f_Y, 
                             mad_q97_logsn_100_f_Y, 
                             mad_q03_logsn_100_m_Y, 
                             mad_q05_logsn_100_m_Y, 
                             mad_q25_logsn_100_m_Y, 
                             mad_q50_logsn_100_m_Y, 
                             mad_q75_logsn_100_m_Y, 
                             mad_q95_logsn_100_m_Y, 
                             mad_q97_logsn_100_m_Y)
mad_quantiles_qr_100 <- c(mad_q03_qr_100_f_Y, 
                          mad_q05_qr_100_f_Y, 
                          mad_q25_qr_100_f_Y, 
                          mad_q50_qr_100_f_Y, 
                          mad_q75_qr_100_f_Y, 
                          mad_q95_qr_100_f_Y, 
                          mad_q97_qr_100_f_Y, 
                          mad_q03_qr_100_m_Y, 
                          mad_q05_qr_100_m_Y, 
                          mad_q25_qr_100_m_Y, 
                          mad_q50_qr_100_m_Y, 
                          mad_q75_qr_100_m_Y, 
                          mad_q95_qr_100_m_Y, 
                          mad_q97_qr_100_m_Y)
table_quant_logsn_qr_100 <- data.frame(true_quantiles_logsn, 
                                       median_quantiles_logsn_100, 
                                       mad_quantiles_logsn_100, 
                                       median_quantiles_qr_100, 
                                       mad_quantiles_qr_100,
                                       row.names = name_quantiles)
table_quant_logsn_qr_100 <- signif(round(table_quant_logsn_qr_100, 4), 5)
names(table_quant_logsn_qr_100) <- c("True", "Median_logsn", "MAD_logsn", "Median_qr", "MAD_qr")

options(scipen = 999)
name_quantiles_f <- paste("Y_", c(t(outer(c("0.03", "0.05", "0.25", "0.50", "0.75", "0.95", "0.97"), "f", FUN = "paste", sep = "_"))), sep = "")
name_quantiles_m <- paste("Y_", c(t(outer(c("0.03", "0.05", "0.25", "0.50", "0.75", "0.95", "0.97"), "m", FUN = "paste", sep = "_"))), sep = "")
name_quantiles <- c(name_quantiles_f, name_quantiles_m)
true_quantiles_logsn <- c(true_q03_f_logsn, 
                          true_q05_f_logsn, 
                          true_q25_f_logsn, 
                          true_q50_f_logsn, 
                          true_q75_f_logsn, 
                          true_q95_f_logsn, 
                          true_q97_f_logsn,
                          true_q03_m_logsn, 
                          true_q05_m_logsn, 
                          true_q25_m_logsn, 
                          true_q50_m_logsn, 
                          true_q75_m_logsn, 
                          true_q95_m_logsn, 
                          true_q97_m_logsn)

median_quantiles_logsn_500 <- c(median_q03_logsn_500_f_Y, 
                                median_q05_logsn_500_f_Y, 
                                median_q25_logsn_500_f_Y, 
                                median_q50_logsn_500_f_Y, 
                                median_q75_logsn_500_f_Y, 
                                median_q95_logsn_500_f_Y, 
                                median_q97_logsn_500_f_Y, 
                                median_q03_logsn_500_m_Y, 
                                median_q05_logsn_500_m_Y, 
                                median_q25_logsn_500_m_Y, 
                                median_q50_logsn_500_m_Y, 
                                median_q75_logsn_500_m_Y, 
                                median_q95_logsn_500_m_Y, 
                                median_q97_logsn_500_m_Y)
median_quantiles_qr_500 <- c(median_q03_qr_500_f_Y, 
                             median_q05_qr_500_f_Y, 
                             median_q25_qr_500_f_Y, 
                             median_q50_qr_500_f_Y, 
                             median_q75_qr_500_f_Y, 
                             median_q95_qr_500_f_Y, 
                             median_q97_qr_500_f_Y, 
                             median_q03_qr_500_m_Y, 
                             median_q05_qr_500_m_Y, 
                             median_q25_qr_500_m_Y, 
                             median_q50_qr_500_m_Y, 
                             median_q75_qr_500_m_Y, 
                             median_q95_qr_500_m_Y, 
                             median_q97_qr_500_m_Y)
mad_quantiles_logsn_500 <- c(mad_q03_logsn_500_f_Y, 
                             mad_q05_logsn_500_f_Y, 
                             mad_q25_logsn_500_f_Y, 
                             mad_q50_logsn_500_f_Y, 
                             mad_q75_logsn_500_f_Y, 
                             mad_q95_logsn_500_f_Y, 
                             mad_q97_logsn_500_f_Y, 
                             mad_q03_logsn_500_m_Y, 
                             mad_q05_logsn_500_m_Y, 
                             mad_q25_logsn_500_m_Y, 
                             mad_q50_logsn_500_m_Y, 
                             mad_q75_logsn_500_m_Y, 
                             mad_q95_logsn_500_m_Y, 
                             mad_q97_logsn_500_m_Y)
mad_quantiles_qr_500 <- c(mad_q03_qr_500_f_Y, 
                          mad_q05_qr_500_f_Y, 
                          mad_q25_qr_500_f_Y, 
                          mad_q50_qr_500_f_Y, 
                          mad_q75_qr_500_f_Y, 
                          mad_q95_qr_500_f_Y, 
                          mad_q97_qr_500_f_Y, 
                          mad_q03_qr_500_m_Y, 
                          mad_q05_qr_500_m_Y, 
                          mad_q25_qr_500_m_Y, 
                          mad_q50_qr_500_m_Y, 
                          mad_q75_qr_500_m_Y, 
                          mad_q95_qr_500_m_Y, 
                          mad_q97_qr_500_m_Y)
table_quant_logsn_qr_500 <- data.frame(true_quantiles_logsn, 
                                       median_quantiles_logsn_500, 
                                       mad_quantiles_logsn_500, 
                                       median_quantiles_qr_500, 
                                       mad_quantiles_qr_500,
                                       row.names = name_quantiles)
table_quant_logsn_qr_500 <- signif(round(table_quant_logsn_qr_500, 4), 5)
names(table_quant_logsn_qr_500) <- c("True", "Median_logsn", "MAD_logsn", "Median_qr", "MAD_qr")

options(scipen = 999)
name_quantiles_f <- paste("Y_", c(t(outer(c("0.03", "0.05", "0.25", "0.50", "0.75", "0.95", "0.97"), "f", FUN = "paste", sep = "_"))), sep = "")
name_quantiles_m <- paste("Y_", c(t(outer(c("0.03", "0.05", "0.25", "0.50", "0.75", "0.95", "0.97"), "m", FUN = "paste", sep = "_"))), sep = "")
name_quantiles <- c(name_quantiles_f, name_quantiles_m)
true_quantiles_logsn <- c(true_q03_f_logsn, 
                          true_q05_f_logsn, 
                          true_q25_f_logsn, 
                          true_q50_f_logsn, 
                          true_q75_f_logsn, 
                          true_q95_f_logsn, 
                          true_q97_f_logsn,
                          true_q03_m_logsn, 
                          true_q05_m_logsn, 
                          true_q25_m_logsn, 
                          true_q50_m_logsn, 
                          true_q75_m_logsn, 
                          true_q95_m_logsn, 
                          true_q97_m_logsn)

median_quantiles_logsn_1000 <- c(median_q03_logsn_1000_f_Y, 
                                 median_q05_logsn_1000_f_Y, 
                                 median_q25_logsn_1000_f_Y, 
                                 median_q50_logsn_1000_f_Y, 
                                 median_q75_logsn_1000_f_Y, 
                                 median_q95_logsn_1000_f_Y, 
                                 median_q97_logsn_1000_f_Y, 
                                 median_q03_logsn_1000_m_Y, 
                                 median_q05_logsn_1000_m_Y, 
                                 median_q25_logsn_1000_m_Y, 
                                 median_q50_logsn_1000_m_Y, 
                                 median_q75_logsn_1000_m_Y, 
                                 median_q95_logsn_1000_m_Y, 
                                 median_q97_logsn_1000_m_Y)
median_quantiles_qr_1000 <- c(median_q03_qr_1000_f_Y, 
                              median_q05_qr_1000_f_Y, 
                              median_q25_qr_1000_f_Y, 
                              median_q50_qr_1000_f_Y, 
                              median_q75_qr_1000_f_Y, 
                              median_q95_qr_1000_f_Y, 
                              median_q97_qr_1000_f_Y, 
                              median_q03_qr_1000_m_Y, 
                              median_q05_qr_1000_m_Y, 
                              median_q25_qr_1000_m_Y, 
                              median_q50_qr_1000_m_Y, 
                              median_q75_qr_1000_m_Y, 
                              median_q95_qr_1000_m_Y, 
                              median_q97_qr_1000_m_Y)
mad_quantiles_logsn_1000 <- c(mad_q03_logsn_1000_f_Y, 
                              mad_q05_logsn_1000_f_Y, 
                              mad_q25_logsn_1000_f_Y, 
                              mad_q50_logsn_1000_f_Y, 
                              mad_q75_logsn_1000_f_Y, 
                              mad_q95_logsn_1000_f_Y, 
                              mad_q97_logsn_1000_f_Y, 
                              mad_q03_logsn_1000_m_Y, 
                              mad_q05_logsn_1000_m_Y, 
                              mad_q25_logsn_1000_m_Y, 
                              mad_q50_logsn_1000_m_Y, 
                              mad_q75_logsn_1000_m_Y, 
                              mad_q95_logsn_1000_m_Y, 
                              mad_q97_logsn_1000_m_Y)
mad_quantiles_qr_1000 <- c(mad_q03_qr_1000_f_Y, 
                           mad_q05_qr_1000_f_Y, 
                           mad_q25_qr_1000_f_Y, 
                           mad_q50_qr_1000_f_Y, 
                           mad_q75_qr_1000_f_Y, 
                           mad_q95_qr_1000_f_Y, 
                           mad_q97_qr_1000_f_Y, 
                           mad_q03_qr_1000_m_Y, 
                           mad_q05_qr_1000_m_Y, 
                           mad_q25_qr_1000_m_Y, 
                           mad_q50_qr_1000_m_Y, 
                           mad_q75_qr_1000_m_Y, 
                           mad_q95_qr_1000_m_Y, 
                           mad_q97_qr_1000_m_Y)
table_quant_logsn_qr_1000 <- data.frame(true_quantiles_logsn, 
                                        median_quantiles_logsn_1000, 
                                        mad_quantiles_logsn_1000, 
                                        median_quantiles_qr_1000, 
                                        mad_quantiles_qr_1000,
                                        row.names = name_quantiles)
table_quant_logsn_qr_1000 <- signif(round(table_quant_logsn_qr_1000, 4), 5)
names(table_quant_logsn_qr_1000) <- c("True", "Median_logsn", "MAD_logsn", "Median_qr", "MAD_qr")

