

## single run examples
id <- 0
rho <- 0.5
s2n <- 1
decay <- 100
sparse <- TRUE

source("code/simdata.R")
library(gamlr)

## gamma lasso
gl0 <- cv.gamlr(d$x, d$y.train)
gl1 <- cv.gamlr(d$x, d$y.train, gamma=1)
gl10 <- cv.gamlr(d$x, d$y.train, gamma=10)

pdf(file="sim_paths.pdf", width=7, height=2.5)
par(mfrow=c(1,3), 
    mai=c(.4,.4,.3,.2), 
    omi=c(.2,.2,0,0))
for(mod in list(gl0,gl1,gl10)){
    plot(mod$g, xlab="", ylab="", select=FALSE, 
        col=rgb(.5,.5,.5,.5), xlim=c(-3.15,max(log(gl0$gam$lam))), df=FALSE)
    dfi <- c(1,20,40,60,80)
    axis(3, at = log(mod$gamlr$lambda)[dfi], 
        labels = round(mod$gamlr$df[dfi]), tick = FALSE, line = -0.5)
     }
mtext(side=1, "log lambda", 
    font=3, outer=TRUE, cex=.7)
mtext(side=2, "coefficient", 
    font=3, outer=TRUE, cex=.7)
dev.off()

pdf(file="sim_cv.pdf", width=7, height=2.5)
par(mfrow=c(1,3), 
    mai=c(.4,.4,.3,.1), 
    omi=c(.2,.2,0,.2))
for(mod in list(gl0,gl1,gl10)){
    plot(mod, xlab="", ylab="", col="grey50",
        ylim=c(12,30), xlim=c(-3.15,max(log(gl0$gam$lam))), 
        select=FALSE, df=FALSE) 
    lines( log(mod$g$lam), exp(AICc(mod$g)/nrow(d$x)), lwd=2)
    dfi <- c(1,20,40,60,80)
    axis(3, at = log(mod$gamlr$lambda)[dfi], 
        labels = round(mod$gamlr$df[dfi]), tick = FALSE, line = -0.5)
}
mtext(side=1, "log lambda", 
    font=3, outer=TRUE, cex=.7)
mtext(side=2, "error", 
    font=3, outer=TRUE, cex=.7)
legend("topright",bty="n", lwd=3,
    col=c("grey75","black"), 
    legend=c("CV MSE","exp(AICc/n)"))
dev.off()
