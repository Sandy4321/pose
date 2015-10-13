## extract sim var
args <- commandArgs(TRUE)
id <- as.integer(args[1])
rho <- as.numeric(args[2])
s2n <- as.numeric(args[3])
decay <- as.numeric(args[4])

OUT=paste("sim-rho",rho,"-s2n",s2n,"-decay",decay,collapse='',sep='')
print(sessionInfo())

source("code/simdata.R")
library(gamlr)
source("code/simfit.R")
print(id)

fill <- function(v){ c(v, rep(tail(v,1),100-length(v)) )}

## basic properties
write(d$sigma, sprintf("results/%s-sigma.txt",OUT),append=TRUE)
times <- paste(round(c(tgl0,tgl2,tgl10,tmrgal,tsnet),1),collapse="|")
write(times, sprintf("results/%s-times.txt",OUT),append=TRUE)
write(paste(gl0$gamlr$lambda,collapse="|"),
        sprintf("results/%s-lambda.txt",OUT),append=TRUE)

## prediction
e0 <- predict(gl0$gamlr,d$x,select=0)-d$y.val
e2 <- predict(gl2$gamlr,d$x,select=0)-d$y.val
e10 <- predict(gl2$gamlr,d$x,select=0)-d$y.val
emrg <- predict(mrgal$gamlr,d$x,select=0)-d$y.val
ecp <- predict(cpbest,as.data.frame(as.matrix(d$x)))-d$y.val
esnet1se <- predict(snet,as.matrix(d$x),which="parms.1se")-d$y.val
esnetmin <- predict(snet,as.matrix(d$x),which="parms.min")-d$y.val

mse0 <- apply(e0^2,2,mean)
mse2 <- apply(e2^2,2,mean)
mse10 <- apply(e10^2,2,mean)
msemrg <- apply(emrg^2,2,mean)
msecp <- mean(ecp^2)
msesnet1se <- mean(esnet1se^2)
msesnetmin <- mean(esnetmin^2)

seg0 <- c(onese=gl0$seg.1se,min=gl0$seg.min,
          aicc=which.min(AICc(gl0$gamlr)),
          aic=which.min(AIC(gl0$gamlr)),
          bic=which.min(BIC(gl0$gamlr)))
seg2 <- c(onese=gl2$seg.1se,min=gl2$seg.min,
          aicc=which.min(AICc(gl2$gamlr)),
          aic=which.min(AIC(gl2$gamlr)),
          bic=which.min(BIC(gl2$gamlr)))
seg10 <- c(onese=gl10$seg.1se,min=gl10$seg.min,
           aicc=which.min(AICc(gl10$gamlr)),
           aic=which.min(AIC(gl10$gamlr)),
           bic=which.min(BIC(gl10$gamlr)))
segmrg <- c(onese=mrgal$seg.1se,min=mrgal$seg.min,
            aicc=which.min(AICc(mrgal$gamlr)),
            aic=which.min(AIC(mrgal$gamlr)),
            bic=which.min(BIC(mrgal$gamlr)))
write(paste(c(segmrg,seg0,seg2,seg10),collapse="|"),
      sprintf("results/%s-seg.txt",OUT),append=TRUE)

MSE <- c(cp=msecp,
         snet1se=msesnet1se,
         snetmin=msesnetmin,
         mrg=msemrg[segmrg],gl0=mse0[seg0],
         gl2=mse2[seg2],gl10=mse10[seg10])
write(paste(round(MSE,2),collapse="|"),
      sprintf("results/%s-mse.txt",OUT),append=TRUE)


R2 <- 1-MSE/mean( (d$y.val-mean(d$y.val))^2 )
write(paste(round(R2,2),collapse="|"),
      sprintf("results/%s-r2.txt",OUT),append=TRUE)

MSElong <- c( 
  mrg=fill(msemrg),gl0=fill(mse0),
  gl2=fill(mse2),gl10=fill(mse10))
write(paste(round(MSElong,2),collapse="|"),
      sprintf("results/%s-mselong.txt",OUT),append=TRUE)

R2long <- 1-MSElong/mean( (d$y.val-mean(d$y.val))^2 )
write(paste(round(R2long,2),collapse="|"),
      sprintf("results/%s-r2long.txt",OUT),append=TRUE)


## tracking the weights
wmrg <- matrix(wmrg,nrow=ncol(d$x),ncol=100)
w0 <- matrix(1,nrow=ncol(d$x),ncol=100)

getw <- function(fit){
  b <- coef(fit$gamlr,select=0)[-1,]
  gam <- fit$gamlr$gamma
  w <- matrix(1, nrow=nrow(b),ncol=ncol(b))
  if(gam!=0) w[,-1] <- 1/(1+gam*abs(as.matrix(b)[,-ncol(b)]))
  return(w) }
w2 <- getw(gl2)
w10 <- getw(gl10)

nu <- d$sigma^2/nrow(d$x)

S <- 1:sCp
writeE <- function(W,lam,f){  
  wsnorm <- apply(W[S,],2,
                  function(w) sqrt(sum(w^2)))/sqrt(sCp)
  wmin <- apply(W[-S,],2,min)
  L <- round(wsnorm/(wmin - sqrt(2*nu)/lam),2)

  write(paste(L,collapse="|"),
        sprintf("results/%s-L%s.txt",OUT,f),append=TRUE)
  write(paste(wsnorm,collapse="|"),
   	    sprintf("results/%s-wsnorm%s.txt",OUT,f),append=TRUE)
  write(paste(wmin,collapse="|"),
        sprintf("results/%s-wmin%s.txt",OUT,f),append=TRUE)
}

writeE(wmrg,mrgal$gamlr$lambda,"mrg")
writeE(w0,gl0$gamlr$lambda,"gl0")
writeE(w2,gl2$gamlr$lambda,"gl2")
writeE(w10,gl10$gamlr$lambda,"gl10")

source("code/simsens.R")








