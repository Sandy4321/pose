## extract results
OUT <- "0"

## data properties
write(d$sigma, "results/sigma.txt",append=TRUE)
times <- paste(round(c(tgl0,tgl2,tgl10,tmrgal,tsnet,tscad),1),collapse="|")
write(times, "results/times.txt",append=TRUE)

## prediction
e0 <- predict(gl0$gamlr,d$x,select=0)-d$y.val
e2 <- predict(gl2$gamlr,d$x,select=0)-d$y.val
e10 <- predict(gl2$gamlr,d$x,select=0)-d$y.val
emrg <- predict(mrgal$gamlr,d$x,select=0)-d$y.val
ecp <- predict(cpbest,as.data.frame(as.matrix(d$x)))-d$y.val
esnet <- predict(snet,as.matrix(d$x),which="parms.min")-d$y.val
escad <- predict(mrgal,as.matrix(d$x))-d$y.val

mse0 <- apply(e0^2,2,mean)
mse2 <- apply(e2^2,2,mean)
mse10 <- apply(e10^2,2,mean)
msemrg <- apply(emrg^2,2,mean)
msecp <- mean(ecp^2)
msesnet <- mean(esnet^2)
msescad <- mean(escad^2) # 9.744

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

MSE <- c(cp=msecp,snet=msesnet,scad=msescad,
	mrg=msemrg[segmrg],gl0=mse0[seg0],
	gl2=mse2[seg2],gl10=mse10[seg10])
write(paste(round(MSE,2),collapse="|"),
	sprintf("results/%s/mse.txt",OUT),append=TRUE)

R2 <- 1-MSE/mean( (d$y.val-mean(d$y.val))^2 )
write(paste(round(R2,2),collapse="|"),
	sprintf("results/%s/r2.txt",OUT),append=TRUE)

## support and sign recovery
getsupport <- function(fit){
	b <- coef(fit$gamlr,select=0)[-1,]
	apply(b,2,function(c) which(c!=0))
}

S0 <- getsupport(gl0)
S2 <- getsupport(gl2)
S10 <- getsupport(gl10)
Smrg <- getsupport(mrgal)
Ssnet <- which(coef(snet)[-1,]!=0)
Sscad <- which(coef(scad)[-1]!=0)

s0 <- sapply(S0,length)
s2 <- sapply(S2,length)
s10 <- sapply(S10,length)
smrg <- sapply(Smrg,length)
ssnet <- length(Ssnet) 
sscad <- length(Sscad) # 137

s <- c(cp=sCp,snet=ssnet,scad=sscad,mrg=smrg[segmrg],
	gl0=s0[seg0],gl2=s2[seg2],gl10=s10[seg10])
write(paste(round(s,2),collapse="|"),
	sprintf("results/%s/s.txt",OUT),append=TRUE)

fp0 <- sapply(S0,function(set) sum(set>sCp))
fp2 <- sapply(S2,function(set) sum(set>sCp))
fp10 <- sapply(S10,function(set) sum(set>sCp))
fpmrg <- sapply(Smrg,function(set) sum(set>sCp))
fpsnet <- sum(Ssnet>sCp)
fpscad <- sum(Ssnet>sCp) # 96

fdr <- c(cp=0, snet=fpsnet/ssnet,scad=fpscad/sscad,
	mrg=fpmrg[segmrg]/smrg,gl0=fp0[seg0]/s0,
	gl2=fp2[seg2]/s2,gl10=fp10[seg10]/s10)
write(paste(round(fdr,2),collapse="|"),
	sprintf("results/%s/fdr.txt",OUT),append=TRUE)

tp0 <- sapply(S0,function(set) sum(set<=sCp))
tp2 <- sapply(S2,function(set) sum(set<=sCp))
tp10 <- sapply(S10,function(set) sum(set<=sCp))
tpmrg <- sapply(Smrg,function(set) sum(set<=sCp))
tpsnet <- sum(Ssnet<=sCp)
tpscad <- sum(Ssnet<=sCp) # 48

sens <- c(cp=1,snet=tpsnet/sCp,scad=tpscad/sCp,
	mrg=tpmrg[segmrg]/sCp,gl0=tp0[seg0]/sCp,
	gl2=tp2[seg2]/sCp,gl10=tp10[seg10]/sCp)
write(paste(round(sens,2),collapse="|"),
	sprintf("results/%s/sens.txt",OUT),append=TRUE)

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
sgnsnet <- getsign( coef(snet)[-1,] )
sgnscad <- getsign( coef(scad)[-1] ) # 0.6569343

sgn <- c(cp=sgncp,snet=sgnsnet,scad=sgnscad,
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
	irrep <- round(apply(R,2,function(r) mean(r<0)),2)

	write(paste(L,collapse="|"),
		sprintf("results/%s/L%s.txt",OUT,f),append=TRUE)
	write(paste(cpineq,collapse="|"),
		sprintf("results/%s/cpineq%s.txt",OUT,f),append=TRUE)
	write(paste(irrep,collapse="|"),
		sprintf("results/%s/irrep%s.txt",OUT,f),append=TRUE)
	return(list(L=L,cpineq=cpineq,irrep=irrep))
}

Emrg <- getE(wmrg,mrgal$gamlr$lambda,"mrg")
E0 <- getE(w0,gl0$gamlr$lambda,"gl0")
E2 <- getE(w2,gl2$gamlr$lambda,"gl2")
E10 <- getE(w10,gl10$gamlr$lambda,"gl10")













