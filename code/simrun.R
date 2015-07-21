## extract sim var
args <- commandArgs(TRUE)
id <- as.integer(args[1])
rho <- as.numeric(args[2])
s2n <- as.numeric(args[3])
OUT <- args[4]

print(sessionInfo())

source("code/simdata.R")
library(gamlr)
source("code/simfit.R")
print(id)

fill <- function(v){ c(v, rep(tail(v,1),100-length(v)) )}

## data properties
write(d$sigma, sprintf("results/%s/sigma.txt",OUT),append=TRUE)
times <- paste(round(c(tgl0,tgl2,tgl10,tmrgal,tsnet),1),collapse="|")
write(times, sprintf("results/%s/times.txt",OUT),append=TRUE)

## lambda grids
gllam <- paste(gl0$gamlr$lambda,collapse="|")
write(gllam, sprintf("results/%s/gllam.txt",OUT),append=TRUE)
mrglam <- paste(mrgal$gamlr$lambda,collapse="|")
write(mrglam, sprintf("results/%s/mrglam.txt",OUT),append=TRUE)

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
	sprintf("results/%s/seg.txt",OUT),append=TRUE)

MSE <- c(cp=msecp,
	snet1se=msesnet1se,
	snetmin=msesnetmin,
	mrg=msemrg[segmrg],gl0=mse0[seg0],
	gl2=mse2[seg2],gl10=mse10[seg10])
write(paste(round(MSE,2),collapse="|"),
	sprintf("results/%s/mse.txt",OUT),append=TRUE)


R2 <- 1-MSE/mean( (d$y.val-mean(d$y.val))^2 )
write(paste(round(R2,2),collapse="|"),
	sprintf("results/%s/r2.txt",OUT),append=TRUE)

MSElong <- c( 
	mrg=fill(msemrg),gl0=fill(mse0),
	gl2=fill(mse2),gl10=fill(mse10))
write(paste(round(MSElong,2),collapse="|"),
	sprintf("results/%s/mselong.txt",OUT),append=TRUE)

R2long <- 1-MSElong/mean( (d$y.val-mean(d$y.val))^2 )
write(paste(round(R2long,2),collapse="|"),
	sprintf("results/%s/r2long.txt",OUT),append=TRUE)

## support and sign recovery
getsupport <- function(fit){
	b <- coef(fit$gamlr,select=0)[-1,]
	apply(b,2,function(c) which(c!=0))
}

S0 <- getsupport(gl0)
S2 <- getsupport(gl2)
S10 <- getsupport(gl10)
Smrg <- getsupport(mrgal)
Ssnet1se <- which(coef(snet, which="parms.1se")[-1,]!=0)
Ssnetmin <- which(coef(snet, which="parms.min")[-1,]!=0)

s0 <- sapply(S0,length)
s2 <- sapply(S2,length)
s10 <- sapply(S10,length)
smrg <- sapply(Smrg,length)
ssnet1se <- length(Ssnet1se) 
ssnetmin <- length(Ssnetmin) 

s <- c(cp=sCp,snet1se=ssnet1se,
	snetmin=ssnetmin, mrg=smrg[segmrg],
	gl0=s0[seg0],gl2=s2[seg2],gl10=s10[seg10])
write(paste(round(s,2),collapse="|"),
	sprintf("results/%s/s.txt",OUT),append=TRUE)

slong <- c(mrg=fill(smrg),
	gl0=fill(s0),gl2=fill(s2),gl10=fill(s10))
write(paste(round(slong,2),collapse="|"),
	sprintf("results/%s/slong.txt",OUT),append=TRUE)

fp0 <- sapply(S0,function(set) sum(set>sCp))
fp2 <- sapply(S2,function(set) sum(set>sCp))
fp10 <- sapply(S10,function(set) sum(set>sCp))
fpmrg <- sapply(Smrg,function(set) sum(set>sCp))
fpsnet1se <- sum(Ssnet1se>sCp)
fpsnetmin <- sum(Ssnetmin>sCp)

fdr <- c(cp=0, 
	snet1se=fpsnet1se/ssnet1se,
	snetmin=fpsnetmin/ssnetmin,
	mrg=(fpmrg/smrg)[segmrg],gl0=(fp0/s0)[seg0],
	gl2=(fp2/s2)[seg2],gl10=(fp10/s10)[seg10])
fdr[is.nan(fdr)|is.infinite(fdr)] <- 0
write(paste(round(fdr*100,2),collapse="|"),
	sprintf("results/%s/fdr.txt",OUT),append=TRUE)

fdrlong <- c(mrg=fill(fpmrg/smrg),gl0=fill(fp0/s0),
	gl2=fill(fp2/s2),gl10=fill(fp10/s10))
fdrlong[is.nan(fdrlong)|is.infinite(fdrlong)] <- 0
write(paste(round(fdrlong*100,2),collapse="|"),
	sprintf("results/%s/fdrlong.txt",OUT),append=TRUE)

tp0 <- sapply(S0,function(set) sum(set<=sCp))
tp2 <- sapply(S2,function(set) sum(set<=sCp))
tp10 <- sapply(S10,function(set) sum(set<=sCp))
tpmrg <- sapply(Smrg,function(set) sum(set<=sCp))
tpsnet1se <- sum(Ssnet1se<=sCp)
tpsnetmin <- sum(Ssnetmin<=sCp)

sens <- c(cp=1,
	snet1se=tpsnet1se/sCp,
	snetmin=tpsnetmin/sCp,
	mrg=tpmrg[segmrg]/sCp,gl0=tp0[seg0]/sCp,
	gl2=tp2[seg2]/sCp,gl10=tp10[seg10]/sCp)
write(paste(round(sens*100,2),collapse="|"),
	sprintf("results/%s/sens.txt",OUT),append=TRUE)

senslong <- c(
	mrg=fill(tpmrg/sCp),gl0=fill(tp0/sCp),
	gl2=fill(tp2/sCp),gl10=fill(tp10/sCp))
write(paste(round(senslong*100,2),collapse="|"),
	sprintf("results/%s/senslong.txt",OUT),append=TRUE)

getsign <- function(b){
	b <- as.matrix(b)
	p <- nrow(b)
	apply(b,2,
		function(c){
			cnz <- which(c!=0)
			m <- mean(sign(c[cnz])==(-1)^cnz) 
			if(is.nan(m)) m <- 1
			m
			}
		)
	}
sgn0 <-  getsign( coef(gl0$gamlr,select=0)[-1,] )
sgn2 <-  getsign( coef(gl2$gamlr,select=0)[-1,] )
sgn10 <-  getsign( coef(gl10$gamlr,select=0)[-1,] )
sgnmrg <- getsign( coef(mrgal$gamlr,select=0)[-1,] )
sgncp <- getsign( coef(cpbest)[-1] )
sgnsnet1se <- getsign( coef(snet, which="parms.1se")[-1,] )
sgnsnetmin <- getsign( coef(snet, which="parms.min")[-1,] )

sgn <- c(cp=sgncp,
	snet1se=sgnsnet1se,
	snetmin=sgnsnetmin,
	mrg=sgnmrg[segmrg],gl0=sgn0[seg0],
	gl2=sgn2[seg2],gl10=sgn10[seg10])
write(paste(round(sgn,2),collapse="|"),
	sprintf("results/%s/sgn.txt",OUT),append=TRUE)

## tracking the weights
wmrg <- matrix(wmrg,nrow=ncol(d$x),ncol=100)

getw <- function(fit){
	b <- coef(fit$gamlr,select=0)[-1,]
	gam <- fit$gamlr$gamma
	w <- matrix(1, nrow=nrow(b),ncol=ncol(b))
	if(gam!=0) w[,-1] <- 1/(1+gam*abs(as.matrix(b)[,-ncol(b)]))
	return(w) }

w0 <- matrix(1,nrow=ncol(d$x),ncol=100)
w2 <- getw(gl2)
w10 <- getw(gl10)

nu <- d$sigma^2/nrow(d$x)

S <- 1:sCp
XXXXi <- t(d$x[,-S])%*%d$x[,S]%*%solve(t(d$x[,S])%*%d$x[,S])
getE <- function(W,lam,f){	
	wsnorm <- apply(W[S,],2,
		function(w) sqrt(sum(w^2)))/sqrt(sCp)
	wmin <- apply(W[-S,],2,min)
	cpineq <- as.integer(wmin > sqrt(2*nu)/lam)
	L <- round(wsnorm/(wmin - sqrt(2*nu)/lam),2)
	R <- XXXXi%*%W[S,] - 1 + sqrt(2*nu)/t(lam*t(W[-S,]))
	irrep <- round(apply(R,2,function(r) 100*mean(r<0)),2)

	write(paste(L*cpineq,collapse="|"),
		sprintf("results/%s/L%s.txt",OUT,f),append=TRUE)
	write(paste(cpineq,collapse="|"),
		sprintf("results/%s/cpineq%s.txt",OUT,f),append=TRUE)
	write(paste(irrep,collapse="|"),
		sprintf("results/%s/irrep%s.txt",OUT,f),append=TRUE)
	write(paste(wmin,collapse="|"),
		sprintf("results/%s/wmin%s.txt",OUT,f),append=TRUE)
	return(list(L=L,cpineq=cpineq,irrep=irrep))
}

Emrg <- getE(wmrg,mrgal$gamlr$lambda,"mrg")
E0 <- getE(w0,gl0$gamlr$lambda,"gl0")
E2 <- getE(w2,gl2$gamlr$lambda,"gl2")
E10 <- getE(w10,gl10$gamlr$lambda,"gl10")













