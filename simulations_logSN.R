#====================================================================
# DESCRIPTION: This file conducts simulation studies associated 
# with the LSNLRM. We generate 10000 Monte Carlo replicates of the 
# response for the sample sizes n = 50, 100, 500, 1000.
# This code generates the 'simulations_logSN.RData' file, needed 
# to reproduce Tables 1 and 2 of the paper
#====================================================================

require(sn)
require(tidyverse)
require(quantreg)

#======================
# Required functions
#======================
qsn1 <- Vectorize(qsn, vectorize.args = c("alpha"))

# Loading and pre-procesing the data set for obtaining the real parameters:
data.children <- read_delim("nutricion.csv",delim=";")
colnames(data.children)[colnames(data.children) == "nutricion.peso"] <- "weight"
colnames(data.children)[colnames(data.children) == "nutricion.sexo"] <- "gender"
colnames(data.children)[colnames(data.children) == "nutricion.ano"] <- "year"
colnames(data.children)[colnames(data.children) == "nutricion.edad_dias"] <- "age"
colnames(data.children)[colnames(data.children) == "nutricion.comuna"] <- "commune"
data.children$gender <- recode(data.children$gender, "F" = 0, "M" = 1)
data.children$age <- data.children$age/365
data.children <- data.children %>% filter(year==2018 ,age>=2, age<=5, commune=="Buenos Aires") 
logW <- log(data.children$weight)
n <- nrow(data.children)
X <- data.frame(rep(1, n), data.children$gender, data.children$age)  # real data
names(X) <- c("X0", "X1", "X2")
attach(X)
X_f <- X[X$X1 == 0, ]  # female
X_m <- X[X$X1 == 1, ]  # male
fix.X_f <- apply(X_f, 2, mean)
fix.X_m <- apply(X_m, 2, mean)

# Fit of log-skew-normal linear regression model to children data
fit.logsn.sim <- selm.fit(x=as.matrix(X), y=logW, family = "SN",
                          selm.control=list(method="MLE",
                                            info.type="observed",
                                            opt.method="nlminb",
                                            opt.control=list(abs.tol=1e-10,
                                                             iter.max=2000)))

# True parameters
true_beta1_logsn <- round(as.numeric(fit.logsn.sim$param$dp[1]), 4)
true_beta2_logsn <- round(as.numeric(fit.logsn.sim$param$dp[2]), 4)
true_beta3_logsn <- round(as.numeric(fit.logsn.sim$param$dp[3]), 4)
true_omega_logsn <- round(as.numeric(fit.logsn.sim$param$dp[4]), 4)
true_alpha_logsn <- round(as.numeric(fit.logsn.sim$param$dp[5]), 4)

# Number of Monte Carlo replicates
N <- 10000

#=====================
# Sample size n=50
#=====================

options(warn = -1)
set.seed(1564)
n_50 <- 50
X1_50_logsn <- rbinom(n_50, size = 1, prob = 0.5)
X2_50_logsn <- runif(n_50, min=2, max=5)
X_50_logsn <- cbind(rep(1, n_50), X1_50_logsn, X2_50_logsn)
M <- 1:1e+06
index_50 <- sample(M, N)
alt.index_50 <- setdiff(M, index_50)

estim_logsn_50 <- function(k){
  output <- list(NA)
  ind <- 1
  con <- 0
  while (ind == 1){
    ind1 <- 2
    con1 <- 0
    pos1 <- 0 # location error QR
    pos2 <- 0 # location error logsn-reg
    while (ind1 == 2){
      set.seed(index_50[k])
      logY <- c(NA)
      for (i in 1:n_50) {
        logY[i] <- rsn(1, xi = true_beta1_logsn+true_beta2_logsn*X1_50_logsn[i]+true_beta3_logsn*X2_50_logsn[i], omega = true_omega_logsn, alpha = true_alpha_logsn)
      }
      fit_sim <- try(selm.fit(x=as.matrix(X_50_logsn), y=logY, family = "SN",
                              selm.control=list(method="MLE",
                                                info.type="observed",
                                                opt.method="nlminb",
                                                opt.control=list(abs.tol=1e-10,
                                                                 iter.max=2000))))
      fit_qr_q03_Y <- try(rq(exp(logY) ~ X1_50_logsn + X2_50_logsn, tau = 0.03))
      fit_qr_q05_Y <- try(rq(exp(logY) ~ X1_50_logsn + X2_50_logsn, tau = 0.05))
      fit_qr_q25_Y <- try(rq(exp(logY) ~ X1_50_logsn + X2_50_logsn, tau = 0.25))
      fit_qr_q50_Y <- try(rq(exp(logY) ~ X1_50_logsn + X2_50_logsn, tau = 0.50))
      fit_qr_q75_Y <- try(rq(exp(logY) ~ X1_50_logsn + X2_50_logsn, tau = 0.75))
      fit_qr_q95_Y <- try(rq(exp(logY) ~ X1_50_logsn + X2_50_logsn, tau = 0.95))
      fit_qr_q97_Y <- try(rq(exp(logY) ~ X1_50_logsn + X2_50_logsn, tau = 0.97))
      if (class(fit_sim) == "try-error" | class(fit_qr_q03_Y) == "try-error" | class(fit_qr_q05_Y) == "try-error" | 
          class(fit_qr_q25_Y) == "try-error" | class(fit_qr_q50_Y) == "try-error" | class(fit_qr_q75_Y) == "try-error" | 
          class(fit_qr_q95_Y) == "try-error" | class(fit_qr_q97_Y) == "try-error") {
        pos1 <- pos1 + 1
        set.seed(alt.index_50[pos1+3524])
      }
      ind1 <- ifelse(class(fit_sim) == "try-error", 2, 3)
      ifelse(ind1 == 3, con1 <- con1, con1 <- con1 + 1)
    }
    if (fit_sim$opt.method$convergence != 0) {
      ind <- 1
      pos2 <- pos2 + 1
      set.seed(alt.index_50[pos2+157895])
    } else {
      ind <- 0
    }
    ifelse(ind == 0, con <- con, con <- con + 1)
  }
  output <- list(as.numeric(fit_sim$param$dp[1]),as.numeric(fit_sim$param$dp[2]),as.numeric(fit_sim$param$dp[3]),
                 as.numeric(fit_sim$param$dp[4]),as.numeric(fit_sim$param$dp[5]),fit_sim$param.var$dp,
                 fit_qr_q03_Y$coefficients, fit_qr_q05_Y$coefficients, fit_qr_q25_Y$coefficients, fit_qr_q50_Y$coefficients, fit_qr_q75_Y$coefficients, 
                 fit_qr_q95_Y$coefficients, fit_qr_q97_Y$coefficients, k, con, con1, pos1, pos2)
  names(output) <- c("beta1", "beta2", "beta3", "omega", "alpha", "acov", "beta_q03_Y", "beta_q05_Y", "beta_q25_Y",
                     "beta_q50_Y", "beta_q75_Y", "beta_q95_Y", "beta_q97_Y", "index", "con.div", "error", "location_3524", "location_157895")
  return(output)
}

simulation_logsn_50 <- list()
div_logsn_50 <- 0
error_logsn_50 <- 0
loc_qr_50 <- 0
loc_logsn_50 <- 0

for (i in 1:N){
  simulation_logsn_50[[i]] <- estim_logsn_50(i)
  div_logsn_50 <- div_logsn_50 + simulation_logsn_50[[i]]$con.div
  error_logsn_50 <- error_logsn_50 + simulation_logsn_50[[i]]$error
  loc_qr_50 <- loc_qr_50 + simulation_logsn_50[[i]]$location_3524
  loc_logsn_50 <- loc_logsn_50 + simulation_logsn_50[[i]]$location_157895
  print(c(n_50, i, div_logsn_50, error_logsn_50, loc_qr_50, loc_logsn_50))
}

estim_beta1_logsn_50 <- unlist(lapply(simulation_logsn_50, "[[", 1))
estim_beta2_logsn_50 <- unlist(lapply(simulation_logsn_50, "[[", 2))
estim_beta3_logsn_50 <- unlist(lapply(simulation_logsn_50, "[[", 3))
estim_omega_logsn_50 <- unlist(lapply(simulation_logsn_50, "[[", 4))
estim_alpha_logsn_50 <- unlist(lapply(simulation_logsn_50, "[[", 5))

# 3%, 5%, 25%, 50%, 75%, 95% and 97% cuantiles of the univariate 
# standard skew-normal distribution
estim_q03_stdsn_50 <- qsn1(0.03, xi=0, omega=1, alpha=estim_alpha_logsn_50[1:N], solver="RFB")
estim_q05_stdsn_50 <- qsn1(0.05, xi=0, omega=1, alpha=estim_alpha_logsn_50[1:N], solver="RFB")
estim_q25_stdsn_50 <- qsn1(0.25, xi=0, omega=1, alpha=estim_alpha_logsn_50[1:N], solver="RFB")
estim_q50_stdsn_50 <- qsn1(0.50, xi=0, omega=1, alpha=estim_alpha_logsn_50[1:N], solver="RFB")
estim_q75_stdsn_50 <- qsn1(0.75, xi=0, omega=1, alpha=estim_alpha_logsn_50[1:N], solver="RFB")
estim_q95_stdsn_50 <- qsn1(0.95, xi=0, omega=1, alpha=estim_alpha_logsn_50[1:N], solver="RFB")
estim_q97_stdsn_50 <- qsn1(0.97, xi=0, omega=1, alpha=estim_alpha_logsn_50[1:N], solver="RFB")

#=====================
# Sample size n=100
#=====================

options(warn = -1)
set.seed(0231)
n_100 <- 100
X1_100_logsn <- rbinom(n_100, size = 1, prob = 0.5)
X2_100_logsn <- runif(n_100, min=2, max=5)
X_100_logsn <- cbind(rep(1, n_100), X1_100_logsn, X2_100_logsn)
M <- 1:1e+06
index_100 <- sample(M, N)
alt.index_100 <- setdiff(M, index_100)

estim_logsn_100 <- function(k){
  output <- list(NA)
  ind <- 1
  con <- 0
  while (ind == 1){
    ind1 <- 2
    con1 <- 0
    pos1 <- 0 # location error QR
    pos2 <- 0 # location error logsn-reg
    while (ind1 == 2){
      set.seed(index_100[k])
      logY <- c(NA)
      for (i in 1:n_100) {
        logY[i] <- rsn(1, xi = true_beta1_logsn+true_beta2_logsn*X1_100_logsn[i]+true_beta3_logsn*X2_100_logsn[i], omega = true_omega_logsn, alpha = true_alpha_logsn)
      }
      fit_sim <- try(selm.fit(x=as.matrix(X_100_logsn), y=logY, family = "SN",
                              selm.control=list(method="MLE",
                                                info.type="observed",
                                                opt.method="nlminb",
                                                opt.control=list(abs.tol=1e-10,
                                                                 iter.max=2000))))
      fit_qr_q03_Y <- try(rq(exp(logY) ~ X1_100_logsn + X2_100_logsn, tau = 0.03))
      fit_qr_q05_Y <- try(rq(exp(logY) ~ X1_100_logsn + X2_100_logsn, tau = 0.05))
      fit_qr_q25_Y <- try(rq(exp(logY) ~ X1_100_logsn + X2_100_logsn, tau = 0.25))
      fit_qr_q50_Y <- try(rq(exp(logY) ~ X1_100_logsn + X2_100_logsn, tau = 0.50))
      fit_qr_q75_Y <- try(rq(exp(logY) ~ X1_100_logsn + X2_100_logsn, tau = 0.75))
      fit_qr_q95_Y <- try(rq(exp(logY) ~ X1_100_logsn + X2_100_logsn, tau = 0.95))
      fit_qr_q97_Y <- try(rq(exp(logY) ~ X1_100_logsn + X2_100_logsn, tau = 0.97))
      if (class(fit_sim) == "try-error" | class(fit_qr_q03_Y) == "try-error" | class(fit_qr_q05_Y) == "try-error" | 
          class(fit_qr_q25_Y) == "try-error" | class(fit_qr_q50_Y) == "try-error" | class(fit_qr_q75_Y) == "try-error" | 
          class(fit_qr_q95_Y) == "try-error" | class(fit_qr_q97_Y) == "try-error") {
        pos1 <- pos1 + 1
        set.seed(alt.index_100[pos1+3524])
              }
      ind1 <- ifelse(class(fit_sim) == "try-error", 2, 3)
      ifelse(ind1 == 3, con1 <- con1, con1 <- con1 + 1)
    }
    if (fit_sim$opt.method$convergence != 0) {
      ind <- 1
      pos2 <- pos2 + 1
      set.seed(alt.index_100[pos2+157895])
  } else {
      ind <- 0
    }
    ifelse(ind == 0, con <- con, con <- con + 1)
  }
  output <- list(as.numeric(fit_sim$param$dp[1]),as.numeric(fit_sim$param$dp[2]),as.numeric(fit_sim$param$dp[3]),
                 as.numeric(fit_sim$param$dp[4]),as.numeric(fit_sim$param$dp[5]),fit_sim$param.var$dp,
                 fit_qr_q03_Y$coefficients, fit_qr_q05_Y$coefficients, fit_qr_q25_Y$coefficients, fit_qr_q50_Y$coefficients, fit_qr_q75_Y$coefficients, 
                 fit_qr_q95_Y$coefficients, fit_qr_q97_Y$coefficients, k, con, con1, pos1, pos2)
  names(output) <- c("beta1", "beta2", "beta3", "omega", "alpha", "acov", "beta_q03_Y", "beta_q05_Y", "beta_q25_Y",
                     "beta_q50_Y", "beta_q75_Y", "beta_q95_Y", "beta_q97_Y", "index", "con.div", "error", "location_3524", "location_157895")
  return(output)
}

simulation_logsn_100 <- list()
div_logsn_100 <- 0
error_logsn_100 <- 0
loc_qr_100 <- 0
loc_logsn_100 <- 0

for (i in 1:N){
  simulation_logsn_100[[i]] <- estim_logsn_100(i)
  div_logsn_100 <- div_logsn_100 + simulation_logsn_100[[i]]$con.div
  error_logsn_100 <- error_logsn_100 + simulation_logsn_100[[i]]$error
  loc_qr_100 <- loc_qr_100 + simulation_logsn_100[[i]]$location_3524
  loc_logsn_100 <- loc_logsn_100 + simulation_logsn_100[[i]]$location_157895
  print(c(n_100, i, div_logsn_100, error_logsn_100, loc_qr_100, loc_logsn_100))
}

estim_beta1_logsn_100 <- unlist(lapply(simulation_logsn_100, "[[", 1))
estim_beta2_logsn_100 <- unlist(lapply(simulation_logsn_100, "[[", 2))
estim_beta3_logsn_100 <- unlist(lapply(simulation_logsn_100, "[[", 3))
estim_omega_logsn_100 <- unlist(lapply(simulation_logsn_100, "[[", 4))
estim_alpha_logsn_100 <- unlist(lapply(simulation_logsn_100, "[[", 5))

# 3%, 5%, 25%, 50%, 75%, 95% and 97% cuantiles of the univariate 
# standard skew-normal distribution
estim_q03_stdsn_100 <- qsn1(0.03, xi=0, omega=1, alpha=estim_alpha_logsn_100[1:N], solver="RFB")
estim_q05_stdsn_100 <- qsn1(0.05, xi=0, omega=1, alpha=estim_alpha_logsn_100[1:N], solver="RFB")
estim_q25_stdsn_100 <- qsn1(0.25, xi=0, omega=1, alpha=estim_alpha_logsn_100[1:N], solver="RFB")
estim_q50_stdsn_100 <- qsn1(0.50, xi=0, omega=1, alpha=estim_alpha_logsn_100[1:N], solver="RFB")
estim_q75_stdsn_100 <- qsn1(0.75, xi=0, omega=1, alpha=estim_alpha_logsn_100[1:N], solver="RFB")
estim_q95_stdsn_100 <- qsn1(0.95, xi=0, omega=1, alpha=estim_alpha_logsn_100[1:N], solver="RFB")
estim_q97_stdsn_100 <- qsn1(0.97, xi=0, omega=1, alpha=estim_alpha_logsn_100[1:N], solver="RFB")

#=====================
# Sample size n=500
#=====================

options(warn = -1)
set.seed(7789)
n_500 <- 500
X1_500_logsn <- rbinom(n_500, size = 1, prob = 0.5)
X2_500_logsn <- runif(n_500, min=2, max=5)
X_500_logsn <- cbind(rep(1, n_500), X1_500_logsn, X2_500_logsn)
M <- 1:1e+06
index_500 <- sample(M, N)
alt.index_500 <- setdiff(M, index_500)

estim_logsn_500 <- function(k){
  output <- list(NA)
  ind <- 1
  con <- 0
  while (ind == 1){
    ind1 <- 2
    con1 <- 0
    pos1 <- 0 # location error QR
    pos2 <- 0 # location error logsn-reg
    while (ind1 == 2){
      set.seed(index_500[k])
      logY <- c(NA)
      for (i in 1:n_500) {
        logY[i] <- rsn(1, xi = true_beta1_logsn+true_beta2_logsn*X1_500_logsn[i]+true_beta3_logsn*X2_500_logsn[i], omega = true_omega_logsn, alpha = true_alpha_logsn)
      }
      fit_sim <- try(selm.fit(x=as.matrix(X_500_logsn), y=logY, family = "SN",
                              selm.control=list(method="MLE",
                                                info.type="observed",
                                                opt.method="nlminb",
                                                opt.control=list(abs.tol=1e-10,
                                                                 iter.max=2000))))
      fit_qr_q03_Y <- try(rq(exp(logY) ~ X1_500_logsn + X2_500_logsn, tau = 0.03))
      fit_qr_q05_Y <- try(rq(exp(logY) ~ X1_500_logsn + X2_500_logsn, tau = 0.05))
      fit_qr_q25_Y <- try(rq(exp(logY) ~ X1_500_logsn + X2_500_logsn, tau = 0.25))
      fit_qr_q50_Y <- try(rq(exp(logY) ~ X1_500_logsn + X2_500_logsn, tau = 0.50))
      fit_qr_q75_Y <- try(rq(exp(logY) ~ X1_500_logsn + X2_500_logsn, tau = 0.75))
      fit_qr_q95_Y <- try(rq(exp(logY) ~ X1_500_logsn + X2_500_logsn, tau = 0.95))
      fit_qr_q97_Y <- try(rq(exp(logY) ~ X1_500_logsn + X2_500_logsn, tau = 0.97))
      if (class(fit_sim) == "try-error" | class(fit_qr_q03_Y) == "try-error" | class(fit_qr_q05_Y) == "try-error" | 
          class(fit_qr_q25_Y) == "try-error" | class(fit_qr_q50_Y) == "try-error" | class(fit_qr_q75_Y) == "try-error" | 
          class(fit_qr_q95_Y) == "try-error" | class(fit_qr_q97_Y) == "try-error") {
        pos1 <- pos1 + 1
        set.seed(alt.index_500[pos1+3524])
      }
      ind1 <- ifelse(class(fit_sim) == "try-error", 2, 3)
      ifelse(ind1 == 3, con1 <- con1, con1 <- con1 + 1)
    }
    if (fit_sim$opt.method$convergence != 0) {
      ind <- 1
      pos2 <- pos2 + 1
      set.seed(alt.index_500[pos2+157895])
    } else {
      ind <- 0
    }
    ifelse(ind == 0, con <- con, con <- con + 1)
  }
  output <- list(as.numeric(fit_sim$param$dp[1]),as.numeric(fit_sim$param$dp[2]),as.numeric(fit_sim$param$dp[3]),
                 as.numeric(fit_sim$param$dp[4]),as.numeric(fit_sim$param$dp[5]),fit_sim$param.var$dp,
                 fit_qr_q03_Y$coefficients, fit_qr_q05_Y$coefficients, fit_qr_q25_Y$coefficients, fit_qr_q50_Y$coefficients, fit_qr_q75_Y$coefficients, 
                 fit_qr_q95_Y$coefficients, fit_qr_q97_Y$coefficients, k, con, con1, pos1, pos2)
  names(output) <- c("beta1", "beta2", "beta3", "omega", "alpha", "acov", "beta_q03_Y", "beta_q05_Y", "beta_q25_Y",
                     "beta_q50_Y", "beta_q75_Y", "beta_q95_Y", "beta_q97_Y", "index", "con.div", "error", "location_3524", "location_157895")
  return(output)
}

simulation_logsn_500 <- list()
div_logsn_500 <- 0
error_logsn_500 <- 0
loc_qr_500 <- 0
loc_logsn_500 <- 0

for (i in 1:N){
  simulation_logsn_500[[i]] <- estim_logsn_500(i)
  div_logsn_500 <- div_logsn_500 + simulation_logsn_500[[i]]$con.div
  error_logsn_500 <- error_logsn_500 + simulation_logsn_500[[i]]$error
  loc_qr_500 <- loc_qr_500 + simulation_logsn_500[[i]]$location_3524
  loc_logsn_500 <- loc_logsn_500 + simulation_logsn_500[[i]]$location_157895
  print(c(n_500, i, div_logsn_500, error_logsn_500, loc_qr_500, loc_logsn_500))
}

estim_beta1_logsn_500 <- unlist(lapply(simulation_logsn_500, "[[", 1))
estim_beta2_logsn_500 <- unlist(lapply(simulation_logsn_500, "[[", 2))
estim_beta3_logsn_500 <- unlist(lapply(simulation_logsn_500, "[[", 3))
estim_omega_logsn_500 <- unlist(lapply(simulation_logsn_500, "[[", 4))
estim_alpha_logsn_500 <- unlist(lapply(simulation_logsn_500, "[[", 5))

# 3%, 5%, 25%, 50%, 75%, 95% and 97% cuantiles of the univariate 
# standard skew-normal distribution
estim_q03_stdsn_500 <- qsn1(0.03, xi=0, omega=1, alpha=estim_alpha_logsn_500[1:N], solver="RFB")
estim_q05_stdsn_500 <- qsn1(0.05, xi=0, omega=1, alpha=estim_alpha_logsn_500[1:N], solver="RFB")
estim_q25_stdsn_500 <- qsn1(0.25, xi=0, omega=1, alpha=estim_alpha_logsn_500[1:N], solver="RFB")
estim_q50_stdsn_500 <- qsn1(0.50, xi=0, omega=1, alpha=estim_alpha_logsn_500[1:N], solver="RFB")
estim_q75_stdsn_500 <- qsn1(0.75, xi=0, omega=1, alpha=estim_alpha_logsn_500[1:N], solver="RFB")
estim_q95_stdsn_500 <- qsn1(0.95, xi=0, omega=1, alpha=estim_alpha_logsn_500[1:N], solver="RFB")
estim_q97_stdsn_500 <- qsn1(0.97, xi=0, omega=1, alpha=estim_alpha_logsn_500[1:N], solver="RFB")

#=====================
# Sample size n=1000
#=====================

options(warn = -1)
set.seed(8426)
n_1000 <- 1000
X1_1000_logsn <- rbinom(n_1000, size = 1, prob = 0.5)
X2_1000_logsn <- runif(n_1000, min=2, max=5)
X_1000_logsn <- cbind(rep(1, n_1000), X1_1000_logsn, X2_1000_logsn)
M <- 1:1e+06
index_1000 <- sample(M, N)
alt.index_1000 <- setdiff(M, index_1000)

estim_logsn_1000 <- function(k){
  output <- list(NA)
  ind <- 1
  con <- 0
  while (ind == 1){
    ind1 <- 2
    con1 <- 0
    pos1 <- 0 # location error QR
    pos2 <- 0 # location error logsn-reg
    while (ind1 == 2){
      set.seed(index_1000[k])
      logY <- c(NA)
      for (i in 1:n_1000) {
        logY[i] <- rsn(1, xi = true_beta1_logsn+true_beta2_logsn*X1_1000_logsn[i]+true_beta3_logsn*X2_1000_logsn[i], omega = true_omega_logsn, alpha = true_alpha_logsn)
      }
      fit_sim <- try(selm.fit(x=as.matrix(X_1000_logsn), y=logY, family = "SN",
                              selm.control=list(method="MLE",
                                                info.type="observed",
                                                opt.method="nlminb",
                                                opt.control=list(abs.tol=1e-10,
                                                                 iter.max=2000))))
      fit_qr_q03_Y <- try(rq(exp(logY) ~ X1_1000_logsn + X2_1000_logsn, tau = 0.03))
      fit_qr_q05_Y <- try(rq(exp(logY) ~ X1_1000_logsn + X2_1000_logsn, tau = 0.05))
      fit_qr_q25_Y <- try(rq(exp(logY) ~ X1_1000_logsn + X2_1000_logsn, tau = 0.25))
      fit_qr_q50_Y <- try(rq(exp(logY) ~ X1_1000_logsn + X2_1000_logsn, tau = 0.50))
      fit_qr_q75_Y <- try(rq(exp(logY) ~ X1_1000_logsn + X2_1000_logsn, tau = 0.75))
      fit_qr_q95_Y <- try(rq(exp(logY) ~ X1_1000_logsn + X2_1000_logsn, tau = 0.95))
      fit_qr_q97_Y <- try(rq(exp(logY) ~ X1_1000_logsn + X2_1000_logsn, tau = 0.97))
      if (class(fit_sim) == "try-error" | class(fit_qr_q03_Y) == "try-error" | class(fit_qr_q05_Y) == "try-error" | 
          class(fit_qr_q25_Y) == "try-error" | class(fit_qr_q50_Y) == "try-error" | class(fit_qr_q75_Y) == "try-error" | 
          class(fit_qr_q95_Y) == "try-error" | class(fit_qr_q97_Y) == "try-error") {
        pos1 <- pos1 + 1
        set.seed(alt.index_1000[pos1+3524])
      }
      ind1 <- ifelse(class(fit_sim) == "try-error", 2, 3)
      ifelse(ind1 == 3, con1 <- con1, con1 <- con1 + 1)
    }
    if (fit_sim$opt.method$convergence != 0) {
      ind <- 1
      pos2 <- pos2 + 1
      set.seed(alt.index_1000[pos2+157895])
    } else {
      ind <- 0
    }
    ifelse(ind == 0, con <- con, con <- con + 1)
  }
  output <- list(as.numeric(fit_sim$param$dp[1]),as.numeric(fit_sim$param$dp[2]),as.numeric(fit_sim$param$dp[3]),
                 as.numeric(fit_sim$param$dp[4]),as.numeric(fit_sim$param$dp[5]),fit_sim$param.var$dp,
                 fit_qr_q03_Y$coefficients, fit_qr_q05_Y$coefficients, fit_qr_q25_Y$coefficients, fit_qr_q50_Y$coefficients, fit_qr_q75_Y$coefficients, 
                 fit_qr_q95_Y$coefficients, fit_qr_q97_Y$coefficients, k, con, con1, pos1, pos2)
  names(output) <- c("beta1", "beta2", "beta3", "omega", "alpha", "acov", "beta_q03_Y", "beta_q05_Y", "beta_q25_Y",
                     "beta_q50_Y", "beta_q75_Y", "beta_q95_Y", "beta_q97_Y", "index", "con.div", "error", "location_3524", "location_157895")
  return(output)
}

simulation_logsn_1000 <- list()
div_logsn_1000 <- 0
error_logsn_1000 <- 0
loc_qr_1000 <- 0
loc_logsn_1000 <- 0

for (i in 1:N){
  simulation_logsn_1000[[i]] <- estim_logsn_1000(i)
  div_logsn_1000 <- div_logsn_1000 + simulation_logsn_1000[[i]]$con.div
  error_logsn_1000 <- error_logsn_1000 + simulation_logsn_1000[[i]]$error
  loc_qr_1000 <- loc_qr_1000 + simulation_logsn_1000[[i]]$location_3524
  loc_logsn_1000 <- loc_logsn_1000 + simulation_logsn_1000[[i]]$location_157895
  print(c(n_1000, i, div_logsn_1000, error_logsn_1000, loc_qr_1000, loc_logsn_1000))
}

estim_beta1_logsn_1000 <- unlist(lapply(simulation_logsn_1000, "[[", 1))
estim_beta2_logsn_1000 <- unlist(lapply(simulation_logsn_1000, "[[", 2))
estim_beta3_logsn_1000 <- unlist(lapply(simulation_logsn_1000, "[[", 3))
estim_omega_logsn_1000 <- unlist(lapply(simulation_logsn_1000, "[[", 4))
estim_alpha_logsn_1000 <- unlist(lapply(simulation_logsn_1000, "[[", 5))

# 3%, 5%, 25%, 50%, 75%, 95% and 97% cuantiles of the univariate 
# standard skew-normal distribution
estim_q03_stdsn_1000 <- qsn1(0.03, xi=0, omega=1, alpha=estim_alpha_logsn_1000[1:N], solver="RFB")
estim_q05_stdsn_1000 <- qsn1(0.05, xi=0, omega=1, alpha=estim_alpha_logsn_1000[1:N], solver="RFB")
estim_q25_stdsn_1000 <- qsn1(0.25, xi=0, omega=1, alpha=estim_alpha_logsn_1000[1:N], solver="RFB")
estim_q50_stdsn_1000 <- qsn1(0.50, xi=0, omega=1, alpha=estim_alpha_logsn_1000[1:N], solver="RFB")
estim_q75_stdsn_1000 <- qsn1(0.75, xi=0, omega=1, alpha=estim_alpha_logsn_1000[1:N], solver="RFB")
estim_q95_stdsn_1000 <- qsn1(0.95, xi=0, omega=1, alpha=estim_alpha_logsn_1000[1:N], solver="RFB")
estim_q97_stdsn_1000 <- qsn1(0.97, xi=0, omega=1, alpha=estim_alpha_logsn_1000[1:N], solver="RFB")

save.image("simulations_logSN.RData")
