#==========================================================
# DESCRIPTION: This file generates the PDF plot of the
# LSN distribution (Figure 1)
#==========================================================

density_lsn <- function(y,xi,omega,alpha,log = FALSE){
  # y: response
  # xi: scale parameter
  # omega: relative dispersion parameter
  # alpha: shape parameter
  if(y > 0){
   pdf_lsn <- (2/(y*omega))*dnorm((log(y/xi))/omega)*pnorm(alpha*(log(y/xi))/omega)
  
  if(log){
    return(log(pdf_lsn))
  }
  return(pdf_lsn)
  }
  else {
    if(log)
      return(-Inf)
    
    return(0)
  }
}
  
  
xi_A1 <- 1
omega_A1 <- 0.3
alpha_A1 <- 1
density_A1 <- Vectorize(function(y){
  density_lsn(y=y,xi=xi_A1,omega=omega_A1,alpha=alpha_A1)
})

xi_A2 <- 3
omega_A2 <- 0.3
alpha_A2 <- 1
density_A2 <- Vectorize(function(y){
  density_lsn(y=y,xi=xi_A2,omega=omega_A2,alpha=alpha_A2)
})

xi_A3 <- 5
omega_A3 <- 0.3
alpha_A3 <- 1
density_A3 <- Vectorize(function(y){
  density_lsn(y=y,xi=xi_A3,omega=omega_A3,alpha=alpha_A3)
})

xi_A4 <- 7
omega_A4 <- 0.3
alpha_A4 <- 1
density_A4 <- Vectorize(function(y){
  density_lsn(y=y,xi=xi_A4,omega=omega_A4,alpha=alpha_A4)
})

#pdf("dens_A_logSN.pdf", height = 5.9, width = 10.6)
par(mfrow=c(1,1),pty="s",mex=0.4)
par(mar=c(10,10,10,5))
lat_A <- seq(0,12,0.01)
plot(lat_A,density_A1(lat_A),type="l",lty = 1,xaxt="n",yaxt="n",main="",xlab="",ylab="",lwd = 2) #"solid"
lines(lat_A,density_A2(lat_A), lty = 5, lwd = 2) #"longdash"
lines(lat_A,density_A3(lat_A), lty = 2,lwd = 2) #"dashed"
lines(lat_A,density_A4(lat_A), lty = 3,lwd = 2) #"dotted"
axis_x_A1 <- seq(from = 0, to = 12, by = 6)
axis_y_A1 <- seq(from = 0, to = 1.2, by = 0.6)
axis(side = 1, at = axis_x_A1,cex.axis=2,line = 1,tick = FALSE)
axis(side = 2, at = axis_y_A1,cex.axis=2,line = -0.4,tick = FALSE)
title(xlab=expression(Y), line=7, cex.lab=2)
title(ylab="PDF", line=5, cex.lab=2)

legend("topright", legend = c(expression(paste(xi, " = ", 1)), expression(paste(xi, " = ", 3)), 
                              expression(paste(xi, " = ", 5)), expression(paste(xi, " = ", 7))), 
       lty = c(1,5,2,3), cex = 1.5,lwd = 1.2,bty="n")

#dev.off()
#dens_a <- recordPlot()

xi_B1 <- 3
omega_B1 <- 0.2
alpha_B1 <- 1
density_B1 <- Vectorize(function(y){
  density_lsn(y=y,xi=xi_B1,omega=omega_B1,alpha=alpha_B1)
})

xi_B2 <- 3
omega_B2 <- 0.3
alpha_B2 <- 1
density_B2 <- Vectorize(function(y){
  density_lsn(y=y,xi=xi_B2,omega=omega_B2,alpha=alpha_B2)
})

xi_B3 <- 3
omega_B3 <- 0.4
alpha_B3 <- 1
density_B3 <- Vectorize(function(y){
  density_lsn(y=y,xi=xi_B3,omega=omega_B3,alpha=alpha_B3)
})

xi_B4 <- 3
omega_B4 <- 0.5
alpha_B4 <- 1
density_B4 <- Vectorize(function(y){
  density_lsn(y=y,xi=xi_B4,omega=omega_B4,alpha=alpha_B4)
})

#pdf("dens_B_logSN.pdf", height = 5.9, width = 10.6)
par(mfrow=c(1,1),pty="s",mex=0.4)
par(mar=c(10,10,10,5))
lat_B <- seq(0,10,0.01)
plot(lat_B,density_B1(lat_B),type="l",lty = 1,xaxt="n",yaxt="n",main="",xlab="",ylab="",lwd = 2) #"solid"
lines(lat_B,density_B2(lat_B), lty = 5, lwd = 2) #"longdash"
lines(lat_B,density_B3(lat_B), lty = 2,lwd = 2) #"dashed"
lines(lat_B,density_B4(lat_B), lty = 3,lwd = 2) #"dotted"
axis_x_B1 <- seq(from = 0, to = 10, by = 5)
axis_y_B1 <- seq(from = 0, to = 0.6, by = 0.3)
axis(side = 1, at = axis_x_B1,cex.axis=2,line = 1,tick = FALSE)
axis(side = 2, at = axis_y_B1,cex.axis=2,line = -0.4,tick = FALSE)
title(xlab=expression(Y), line=7, cex.lab=2)
title(ylab="PDF", line=5, cex.lab=2)

legend("topright", legend = c(expression(paste(omega, " = ", 0.2)), expression(paste(omega, " = ", 0.3)), 
                              expression(paste(omega, " = ", 0.4)), expression(paste(omega, " = ", 0.5))), 
       lty = c(1,5,2,3), cex = 1.5,lwd = 1.2,bty="n")
#dev.off()

xi_C1 <- 4
omega_C1 <- 0.4
alpha_C1 <- -2
density_C1 <- Vectorize(function(y){
  density_lsn(y=y,xi=xi_C1,omega=omega_C1,alpha=alpha_C1)
})

xi_C2 <- 4
omega_C2 <- 0.4
alpha_C2 <- 0
density_C2 <- Vectorize(function(y){
  density_lsn(y=y,xi=xi_C2,omega=omega_C2,alpha=alpha_C2)
})

xi_C3 <- 4
omega_C3 <- 0.4
alpha_C3 <- 2
density_C3 <- Vectorize(function(y){
  density_lsn(y=y,xi=xi_C3,omega=omega_C3,alpha=alpha_C3)
})

xi_C4 <- 4
omega_C4 <- 0.4
alpha_C4 <- 5
density_C4 <- Vectorize(function(y){
  density_lsn(y=y,xi=xi_C4,omega=omega_C4,alpha=alpha_C4)
})

#pdf("dens_C_logSN.pdf", height = 5.9, width = 10.6)
par(mfrow=c(1,1),pty="s",mex=0.4)
par(mar=c(10,10,10,5))
lat_C <- seq(0,10,0.01)
plot(lat_C,density_C1(lat_C),type="l",lty = 1,xaxt="n",yaxt="n",main="",xlab="",ylab="",lwd = 2) #"solid"
lines(lat_C,density_C2(lat_C), lty = 5, lwd = 2) #"longdash"
lines(lat_C,density_C3(lat_C), lty = 2,lwd = 2) #"dashed"
lines(lat_C,density_C4(lat_C), lty = 3,lwd = 2) #"dotted"
axis_x_C1 <- seq(from = 0, to = 10, by = 5)
axis_y_C1 <- seq(from = 0, to = 0.4, by = 0.2)
axis(side = 1, at = axis_x_C1,cex.axis=2,line = 1,tick = FALSE)
axis(side = 2, at = axis_y_C1,cex.axis=2,line = -0.4,tick = FALSE)
title(xlab=expression(Y), line=7, cex.lab=2)
title(ylab="PDF", line=5, cex.lab=2)

legend("topright", legend = c(expression(paste(alpha, " = ", -2)), expression(paste(alpha, " = ", 0)), 
                              expression(paste(alpha, " = ", 2)), expression(paste(alpha, " = ", 5))), 
       lty = c(1,5,2,3), cex = 1.5,lwd = 1.2,bty="n")
#dev.off()