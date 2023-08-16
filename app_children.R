#============================================================
# DESCRIPTION: This file fits the LSNLRM to children data
#============================================================

#=====================
# Required packages
#=====================
require(tidyverse)
require(sn)

#======================================
# Reading and preparing the data set
#======================================
data.children <- read_delim("nutricion.csv",delim=";")
colnames(data.children)[colnames(data.children) == "nutricion.peso"] <- "weight"
colnames(data.children)[colnames(data.children) == "nutricion.sexo"] <- "gender"
colnames(data.children)[colnames(data.children) == "nutricion.ano"] <- "year"
colnames(data.children)[colnames(data.children) == "nutricion.edad_dias"] <- "age"
colnames(data.children)[colnames(data.children) == "nutricion.comuna"] <- "commune"
data.children$gender <- recode(data.children$gender, "F" = 0, "M" = 1)
data.children$age <- data.children$age/365
data.children$age_cut <- as.factor(cut(data.children$age,
                             include.lowest = T,
                             breaks = seq(2,5,0.5),
                             labels=c("2 - 2.5","2.5 - 3","3 - 3.5",
                                      "3.5 - 4","4 - 4.5","4.5 - 5")))

data.children <- data.children %>% filter(year==2018 ,age>=2, age<=5, commune=="Buenos Aires") 

logW <- log(data.children$weight)
n <- nrow(data.children)
data.female <- data.children[data.children$gender == 0, ]
data.male <- data.children[data.children$gender == 1, ]
weight.female <- data.female$weight
weight.male <- data.male$weight
age.female <- data.female$age
age.male <- data.male$age
X <- data.frame(rep(1, n), data.children$gender, data.children$age)

#========================
# Descriptive analysis 
#========================
data.boxplot <- data.frame(data.children$weight,data.children$gender,data.children$age_cut)
str(data.boxplot)
names(data.boxplot) <- c("weight", "gender", "age_cut")

#pdf("boxplot_logsn.pdf", height = 8.2, width = 16.4)
par(mfrow = c(1, 1))
cols <- c("white", 'gray')
boxplot(weight~gender+age_cut,data=data.boxplot, at = c(1:2, 4:5, 7:8, 10:11, 13:14, 16:17),
        col = cols, xlab = "", ylab = "", xaxt = "n", yaxt="n")
legend("topleft", fill = cols, legend = c("Female","Male"), 
       horiz = F, bty="n",cex=2)
at_1.box <- levels(data.children$age_cut)
at_2.box <- seq(from = 9, to = 31, by = 11)
axis(side = 1, at = c(1.5,4.5,7.5,10.5,13.5,16.5),
     labels=at_1.box, cex.axis = 2.2, line = 0.1, tick = F)
axis(side = 2, at = at_2.box, cex.axis = 2.2, line = -0.8, tick = F)
title(xlab = "Age", line = 3, cex.lab = 2.5)
title(ylab = "Weight", line = 2.4, cex.lab = 2.5)
#dev.off()
#===============================================
# fit of the log-normal linear regression model 
#===============================================
fit.logsn <- selm.fit(x=as.matrix(X), y=logW, family = "SN",
                      selm.control=list(method="MLE",
                                        info.type="observed",
                                        opt.method="nlminb",
                                        opt.control=list(abs.tol=1e-10,
                                                         iter.max=2000)))

fit.logn <- selm.fit(x=as.matrix(X), y=logW, 
                     family = "SN", 
                     fixed.param = list(alpha=0),
                      selm.control=list(method="MLE",
                                        info.type="observed",
                                        opt.method="nlminb",
                                        opt.control=list(abs.tol=1e-10,
                                                         iter.max=2000)))

#==============================================
# Computation of the AIC and BIC for each model
#==============================================
aic.logsn <- 2*as.numeric(fit.logsn$size[3])-2*fit.logsn$logL
aic.logn <- 2*as.numeric(fit.logn$size[3])-2*fit.logn$logL

bic.logsn <- as.numeric(fit.logsn$size[3])*log(n)-2*fit.logsn$logL
bic.logn  <- as.numeric(fit.logn$size[3])*log(n)-2*fit.logn$logL

#=========================================================
# Quantile-quantile plots of Mahalanobis distances with 
# simulated envelope for the log-normal and 
# log-skew-normal linear regression models
#=========================================================

## logSN
#pdf("diag_logsn.pdf", height = 8.2, width = 8.2)
Mahadist.logsn <- sort((fit.logsn$resid.dp/fit.logsn$param$dp[4])^2)
teoquant.logsn <- qchisq((1:n/(n + 1)), df = 1)
e.logsn <- matrix(0, n, 100)
e1.logsn <- numeric(n)
e2.logsn <- numeric(n)
set.seed(4715)
for (i in 1:100){
  est.sam.logsn <- rsn(n, 0, 1, alpha=as.numeric(fit.logsn$param$dp[5]))
  e.logsn[, i] <- sort(est.sam.logsn^2)
}
for (i in 1:n){
  eo.logsn <- sort(e.logsn[i, ])
  e1.logsn[i] <- min(eo.logsn)
  e2.logsn[i] <- max(eo.logsn)
}
par(mfrow = c(1, 1), pty = "s")
plot(teoquant.logsn, Mahadist.logsn, xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "", xlim = c(0, 13), ylim = c(0, 16.5), 
     lty = 1, main = "", pch = 20, cex = 0.8)
at_1.lsn <- seq(from = 0, to = 12, by = 6)
at_2.lsn <- seq(from = 0, to = 16, by = 8)
axis(side = 1, at = at_1.lsn, cex.axis = 2.7, line = 0.1, tick = F)
axis(side = 2, at = at_2.lsn, cex.axis = 2.7, line = -0.8, tick = F)
title(xlab = "Theoretical quantiles", line = 3, cex.lab = 2.5)
title(ylab = "Observed quantiles", line = 2.4, cex.lab = 2.5)
abline(0,1,col = rgb(0.2, 0.2, 0.2), lty = 3, lwd = 1)
par(new = TRUE)
plot(teoquant.logsn, sort(e1.logsn), axes = F, xlab = "", 
     ylab = "", type = "l", xlim = c(0, 13), ylim = c(0, 16.5), 
     lwd = 1)
par(new = TRUE)
plot(teoquant.logsn, sort(e2.logsn), axes = F, xlab = "", 
     ylab = "", type = "l", xlim = c(0, 13), ylim = c(0, 16.5), 
     lwd = 1)
qqplot_logsn <- recordPlot()
#dev.off()

## logN
#pdf("diag_logn.pdf", height = 8.2, width = 8.2)
Mahadist.logn <- sort((fit.logn$resid.dp/fit.logn$param$dp[4])^2)
teoquant.logn <- qchisq((1:n/(n + 1)), df = 1)
e.logn <- matrix(0, n, 100)
e1.logn <- numeric(n)
e2.logn <- numeric(n)
set.seed(4715)
for (i in 1:100){
  est.sam.logn <- rsn(n, 0, 1, alpha=0)
  e.logn[, i] <- sort(est.sam.logn^2)
}
for (i in 1:n){
  eo.logn <- sort(e.logn[i, ])
  e1.logn[i] <- min(eo.logn)
  e2.logn[i] <- max(eo.logn)
}
par(mfrow = c(1, 1), pty = "s")
plot(teoquant.logn, Mahadist.logn, xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "", xlim = c(0, 13), ylim = c(0, 20), 
     lty = 1, main = "", pch = 20, cex = 0.8)
at_1.lsn <- seq(from = 0, to = 12, by = 6)
at_2.lsn <- seq(from = 0, to = 20, by = 10)
axis(side = 1, at = at_1.lsn, cex.axis = 2.7, line = 0.1, tick = F)
axis(side = 2, at = at_2.lsn, cex.axis = 2.7, line = -0.8, tick = F)
title(xlab = "Theoretical quantiles", line = 3, cex.lab = 2.5)
title(ylab = "Observed quantiles", line = 2.4, cex.lab = 2.5)
abline(0,1,col = rgb(0.2, 0.2, 0.2), lty = 3, lwd = 1)
par(new = TRUE)
plot(teoquant.logn, sort(e1.logn), axes = F, xlab = "", 
     ylab = "", type = "l", xlim = c(0, 13), ylim = c(0, 20), 
     lwd = 1)
par(new = TRUE)
plot(teoquant.logn, sort(e2.logn), axes = F, xlab = "", 
     ylab = "", type = "l", xlim = c(0, 13), ylim = c(0, 20), 
     lwd = 1)
qqplot_logn <- recordPlot()
#dev.off()

#==============================================
#Number of points outside the confidence bands
#==============================================
sum(Mahadist.logsn<sort(e1.logsn) | Mahadist.logsn>sort(e2.logsn))
sum(Mahadist.logn<sort(e1.logn) | Mahadist.logn>sort(e2.logn))

#=========================================
# Estimate, SE, 95% confidence interval 
# and Wald test for for the selected model 
#=========================================

# Regression coefficients
logsn.intercept_xi <- as.numeric(fit.logsn$param$dp[1])
logsn.gender_xi <- as.numeric(fit.logsn$param$dp[2])
logsn.age_xi <- as.numeric(fit.logsn$param$dp[3])
logsn.omega <- as.numeric(fit.logsn$param$dp[4]) 
logsn.alpha <- as.numeric(fit.logsn$param$dp[5])

logsn.var.intercept_xi <- fit.logsn$param.var$dp[1,1]
logsn.var.gender_xi <- fit.logsn$param.var$dp[2,2]
logsn.var.age_xi <- fit.logsn$param.var$dp[3,3]
logsn.var.omega <- fit.logsn$param.var$dp[4,4]
logsn.var.alpha <- fit.logsn$param.var$dp[5,5]


logsn.se.intercept_xi <- sqrt(logsn.var.intercept_xi)
logsn.se.gender_xi <- sqrt(logsn.var.gender_xi)
logsn.se.age_xi <- sqrt(logsn.var.age_xi)
logsn.se.omega <- sqrt(logsn.var.omega)
logsn.se.alpha <- sqrt(logsn.var.alpha)

# 95% confidence intervals
z0.975 <- qnorm(0.975)
low.lim.intercept_xi <- logsn.intercept_xi - z0.975 * logsn.se.intercept_xi
up.lim.intercept_xi <- logsn.intercept_xi + z0.975 * logsn.se.intercept_xi
low.lim.gender_xi <- logsn.gender_xi - z0.975 * logsn.se.gender_xi
up.lim.gender_xi <- logsn.gender_xi + z0.975 * logsn.se.gender_xi
low.lim.age_xi <- logsn.age_xi - z0.975 * logsn.se.age_xi
up.lim.age_xi <- logsn.age_xi + z0.975 * logsn.se.age_xi
low.lim.omega <- logsn.omega - z0.975 * logsn.se.omega
up.lim.omega <- logsn.omega + z0.975 * logsn.se.omega
low.lim.alpha <- logsn.alpha - z0.975 * logsn.se.alpha
up.lim.alpha <- logsn.alpha + z0.975 * logsn.se.alpha

# Wald tests
Wald.logsn.intercept_xi <- logsn.intercept_xi^2/logsn.var.intercept_xi
Wald.logsn.gender_xi <- logsn.gender_xi^2/logsn.var.gender_xi
Wald.logsn.age_xi <- logsn.age_xi^2/logsn.var.age_xi
Wald.logsn.omega <- logsn.omega^2/logsn.var.omega
Wald.logsn.alpha <- logsn.alpha^2/logsn.var.alpha


# p-values
pvalue.logsn.intercept_xi <- 1 - pchisq(Wald.logsn.intercept_xi, df = 1)
pvalue.logsn.gender_xi <- 1 - pchisq(Wald.logsn.gender_xi, df = 1)
pvalue.logsn.age_xi <- 1 - pchisq(Wald.logsn.age_xi, df = 1)
pvalue.logsn.omega <- 1 - pchisq(Wald.logsn.omega, df = 1)
pvalue.logsn.alpha <- 1 - pchisq(Wald.logsn.alpha, df = 1)

#=======================================================
# Table (Estimate, asymptotic standard error (SE),  
# lower and upper bounds of asymptotic 95% confidence  
# interval and p-value of the Wald test associated to 
# log-skew-normal linear regression model)        
#=======================================================

estim.weight <- c(logsn.intercept_xi,
                  logsn.gender_xi,
                  logsn.age_xi)
se.weight <- c(logsn.se.intercept_xi,
               logsn.se.gender_xi,
               logsn.se.age_xi)
low.weight <- c(low.lim.intercept_xi,
                low.lim.gender_xi,
                low.lim.age_xi)
up.weight <- c(up.lim.intercept_xi,
               up.lim.gender_xi,
               up.lim.age_xi)
pvalue.weight <- c(pvalue.logsn.intercept_xi,
                   pvalue.logsn.gender_xi,
                   pvalue.logsn.age_xi)
expl.vbles <- c("intercept", "gender", "age")
table_estim.weight <- data.frame(estim.weight,
                                 se.weight,
                                 low.weight,
                                 up.weight,
                                 pvalue.weight,
                                 row.names = expl.vbles)
table_estim.weight <- signif(round(table_estim.weight, 4), 6)
colnames(table_estim.weight) <- c("Estimate", "SE", "Lower", "Upper",
                                  "p-value")
table_estim.weight
#=======================================================
# Fitted quantile curves (for the percentiles 3, 5, 
# 25, 50, 75, 95, 97) for weight vs. age 
# (for girls and boys).
#=======================================================

q.stsn_03 <- qsn(0.03, xi=0, omega=1, alpha=as.numeric(fit.logsn$param$dp[5]), solver="RFB")
q.stsn_05 <- qsn(0.05, xi=0, omega=1, alpha=as.numeric(fit.logsn$param$dp[5]), solver="RFB")
q.stsn_25 <- qsn(0.25, xi=0, omega=1, alpha=as.numeric(fit.logsn$param$dp[5]), solver="RFB")
q.stsn_50 <- qsn(0.50, xi=0, omega=1, alpha=as.numeric(fit.logsn$param$dp[5]), solver="RFB")
q.stsn_75 <- qsn(0.75, xi=0, omega=1, alpha=as.numeric(fit.logsn$param$dp[5]), solver="RFB")
q.stsn_95 <- qsn(0.95, xi=0, omega=1, alpha=as.numeric(fit.logsn$param$dp[5]), solver="RFB")
q.stsn_97 <- qsn(0.97, xi=0, omega=1, alpha=as.numeric(fit.logsn$param$dp[5]), solver="RFB")

# Weight (female)
quant.weight.female.t1wf <- list()
quant.weight.female.t2wf <- list()
quant.weight.male.t1wf <- list()
quant.weight.male.t2wf <- list()
t1wf <- list(seq(2, 3.8, 0.05), seq(2, 3.8, 0.05), seq(2, 3.8, 0.05), seq(2, 3.8, 0.05),
             seq(2, 3.8, 0.05), seq(2, 3.8, 0.05), seq(2, 3.8, 0.05))
t2wf <- list(seq(4, 5, 0.05), seq(4, 5, 0.05),seq(4, 5, 0.05),seq(4, 5, 0.05),
             seq(4, 5, 0.05),seq(4, 5, 0.05),seq(4, 5, 0.05))
quant.weight.female.t1wf[[1]] <- exp(logsn.intercept_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_03)
quant.weight.female.t2wf[[1]] <- exp(logsn.intercept_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_03)
quant.weight.female.t1wf[[2]] <- exp(logsn.intercept_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_05)
quant.weight.female.t2wf[[2]] <- exp(logsn.intercept_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_05)
quant.weight.female.t1wf[[3]] <- exp(logsn.intercept_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_25)
quant.weight.female.t2wf[[3]] <- exp(logsn.intercept_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_25)
quant.weight.female.t1wf[[4]] <- exp(logsn.intercept_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_50)
quant.weight.female.t2wf[[4]] <- exp(logsn.intercept_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_50)
quant.weight.female.t1wf[[5]] <- exp(logsn.intercept_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_75)
quant.weight.female.t2wf[[5]] <- exp(logsn.intercept_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_75)
quant.weight.female.t1wf[[6]] <- exp(logsn.intercept_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_95)
quant.weight.female.t2wf[[6]] <- exp(logsn.intercept_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_95)
quant.weight.female.t1wf[[7]] <- exp(logsn.intercept_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_97)
quant.weight.female.t2wf[[7]] <- exp(logsn.intercept_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_97)

quant.weight.male.t1wf[[1]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_03)
quant.weight.male.t2wf[[1]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_03)
quant.weight.male.t1wf[[2]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_05)
quant.weight.male.t2wf[[2]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_05)
quant.weight.male.t1wf[[3]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_25)
quant.weight.male.t2wf[[3]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_25)
quant.weight.male.t1wf[[4]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_50)
quant.weight.male.t2wf[[4]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_50)
quant.weight.male.t1wf[[5]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_75)
quant.weight.male.t2wf[[5]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_75)
quant.weight.male.t1wf[[6]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_95)
quant.weight.male.t2wf[[6]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_95)
quant.weight.male.t1wf[[7]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t1wf[[1]] + logsn.omega*q.stsn_97)
quant.weight.male.t2wf[[7]] <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*t2wf[[1]] + logsn.omega*q.stsn_97)

#pdf("quant_weight_female.pdf", height = 8.2, width = 8.2)
par(mfrow = c(1, 1), pty = "s")
plot(age.female, weight.female, xaxt = "n", yaxt = "n", xlab = "", 
     ylab = "", xlim = c(2, 5), ylim = c(7,32), lty = 1, main = "", 
     cex = 1.4, col="gray64")
at1wf.weight.female <- seq(from = 2, to = 5, by = 1)
at2wf.weight.female <- seq(from = 7, to = 31, by = 8)
axis(side = 1, at = at1wf.weight.female, cex.axis = 2.7, line = 0.1, 
     tick = F)
axis(side = 2, at = at2wf.weight.female, cex.axis = 2.7, line = -0.8, 
     tick = F)
title(xlab = "Age", line = 3, cex.lab = 2.5)
title(ylab = "Weight (female)", line = 2.4, cex.lab = 2.5)
for (i in 1:7) {
  lines(t1wf[[i]], quant.weight.female.t1wf[[i]], lwd = 3)
  lines(t2wf[[i]], quant.weight.female.t2wf[[i]], lwd = 3)
}
est.weight.female03 <- exp(logsn.intercept_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_03)
text(3.9, est.weight.female03-0.2, expression(paste(3, th)), 
     srt = 14, cex = 1.2)
est.weight.female05 <- exp(logsn.intercept_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_05)
text(3.9, est.weight.female05+0.2, expression(paste(5, th)), 
     srt = 14, cex = 1.2)
est.weight.female25 <- exp(logsn.intercept_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_25)
text(3.9, est.weight.female25+0.08, expression(paste(25, th)), 
     srt = 14, cex = 1.2)
est.weight.female50 <- exp(logsn.intercept_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_50)
text(3.9, est.weight.female50+0.09, expression(paste(50, th)), 
     srt = 16, cex = 1.2)
est.weight.female75 <- exp(logsn.intercept_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_75)
text(3.9, est.weight.female75+0.11, expression(paste(75, th)), 
     srt = 16, cex = 1.2)
est.weight.female95 <- exp(logsn.intercept_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_95)
text(3.9, est.weight.female95+0.005, expression(paste(95, th)), 
     srt = 18, cex = 1.2)
est.weight.female97 <- exp(logsn.intercept_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_97)
text(3.9, est.weight.female97+0.1, expression(paste(97, th)), 
     srt = 18, cex = 1.2)
quantcurves_weightfemale <- recordPlot()
#dev.off()

#pdf("quant_weight_male.pdf", height = 8.2, width = 8.2)
par(mfrow = c(1, 1), pty = "s")
plot(age.male, weight.male, xaxt = "n", yaxt = "n", xlab = "", 
     ylab = "", xlim = c(2, 5), ylim = c(7,32), lty = 1, main = "", 
     cex = 1.4, col="gray64")
at1wf.weight.male <- seq(from = 2, to = 5, by = 1)
at2wf.weight.male <- seq(from = 7, to = 31, by = 8)
axis(side = 1, at = at1wf.weight.male, cex.axis = 2.7, line = 0.1, 
     tick = F)
axis(side = 2, at = at2wf.weight.male, cex.axis = 2.7, line = -0.8, 
     tick = F)
title(xlab = "Age", line = 3, cex.lab = 2.5)
title(ylab = "Weight (male)", line = 2.4, cex.lab = 2.5)
for (i in 1:7) {
  lines(t1wf[[i]], quant.weight.male.t1wf[[i]], lwd = 3)
  lines(t2wf[[i]], quant.weight.male.t2wf[[i]], lwd = 3)
}
est.weight.male03 <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_03)
text(3.9, est.weight.male03-0.2, expression(paste(3, th)), 
     srt = 14, cex = 1.2)
est.weight.male05 <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_05)
text(3.9, est.weight.male05+0.2, expression(paste(5, th)), 
     srt = 14, cex = 1.2)
est.weight.male25 <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_25)
text(3.9, est.weight.male25+0.08, expression(paste(25, th)), 
     srt = 14, cex = 1.2)
est.weight.male50 <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_50)
text(3.9, est.weight.male50+0.09, expression(paste(50, th)), 
     srt = 16, cex = 1.2)
est.weight.male75 <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_75)
text(3.9, est.weight.male75+0.11, expression(paste(75, th)), 
     srt = 16, cex = 1.2)
est.weight.male95 <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_95)
text(3.9, est.weight.male95+0.005, expression(paste(95, th)), 
     srt = 18, cex = 1.2)
est.weight.male97 <- exp(logsn.intercept_xi + logsn.gender_xi + logsn.age_xi*3.9 + logsn.omega*q.stsn_97)
text(3.9, est.weight.male97+0.1, expression(paste(97, th)), 
     srt = 18, cex = 1.2)
quantcurves_weightmale <- recordPlot()
#dev.off()
