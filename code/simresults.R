## extract results

## prediction
e0 <- predict(gl0$gamlr,d$x,select=0)-d$y.val
e2 <- predict(gl2$gamlr,d$x,select=0)-d$y.val
e10 <- predict(gl2$gamlr,d$x,select=0)-d$y.val
emrg <- predict(mrgal,d$x,s="min")-d$y.val
ecp <- predict(cpbest,as.data.frame(as.matrix(d$x)))-d$y.val
esnet <- predict(snet,as.matrix(d$x),which="parms.min")-d$y.val
escad <- predict(mrgal,as.matrix(d$x))-d$y.val

mse0 <- apply(e0^2,2,mean)
mse2 <- apply(e2^2,2,mean)
mse10 <- apply(e10^2,2,mean)
msemrg <- mean(emrg^2)
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

MSE <- c(cp=msecp,mrg=msemrg,snet=msesnet,gl0=mse0[seg0],
	gl2=mse2[seg2],gl10=mse10[seg10])
names(MSE)[-(1:3)] <- 
	paste(rep(c("gl0","gl2","gl10"),each=length(seg0)),
		names(seg0),sep=".")
#              cp             mrg            snet       gl0.onese   gl0.min.seg47 
#        7.906584        9.457747        8.913506       10.642369        8.978226 
#  gl0.aicc.seg41  gl0.aic.seg100   gl0.bic.seg11       gl2.onese   gl2.min.seg47 
#        9.292775       10.897815       13.529001        9.633760        9.075675 
#  gl2.aicc.seg41  gl2.aic.seg100   gl2.bic.seg11      gl10.onese  gl10.min.seg47 
#        8.890019       11.316384       13.521819       10.264028        9.230127 
# gl10.aicc.seg41 gl10.aic.seg100  gl10.bic.seg11 
#        9.150046       11.316384       13.601107 

	R2 <- 1-MSE/mean( (d$y.val-mean(d$y.val))^2 )

## support and sign recovery
getsupport <- function(fit){
	b <- coef(fit$gamlr,select=0)[-1,]
	apply(b,2,function(c) which(c!=0))
}

S0 <- getsupport(gl0)
S2 <- getsupport(gl2)
S10 <- getsupport(gl10)
Smrg <- which(coef(mrgal)[-1,]!=0)
Ssnet <- which(coef(snet)[-1,]!=0)
Sscad <- which(coef(scad)[-1]!=0)

s0 <- sapply(S0,length)
s2 <- sapply(S2,length)
s10 <- sapply(S10,length)
smrg <- length(Smrg)
ssnet <- length(Ssnet) 
sscad <- length(Sscad) # 137

fp0 <- sapply(S0,function(set) sum(set>sCp))
fp2 <- sapply(S2,function(set) sum(set>sCp))
fp10 <- sapply(S10,function(set) sum(set>sCp))
fpmrg <- sum(Smrg>sCp)
fpsnet <- sum(Ssnet>sCp)
fpscad <- sum(Ssnet>sCp) # 96

tp0 <- sapply(S0,function(set) sum(set<=sCp))
tp2 <- sapply(S2,function(set) sum(set<=sCp))
tp10 <- sapply(S10,function(set) sum(set<=sCp))
tpmrg <- sum(Smrg<=sCp)
tpsnet <- sum(Ssnet<=sCp)
tpscad <- sum(Ssnet<=sCp) # 48

# > fpscad/sscad
# [1] 0.7007299
# > tpscad/sCp
# [1] 0.6956522

getsign <- function(b){
	b <- as.matrix(b)
	p <- nrow(b)
	apply(b,2,
		function(c){
			cnz <- which(c!=0)
			mean(sign(c[cnz])==(-1)^cnz) 
			}
		)
	}
sgn0 <-  getsign( coef(gl0$gamlr,select=0)[-1,] )
sgn2 <-  getsign( coef(gl2$gamlr,select=0)[-1,] )
sgn10 <-  getsign( coef(gl10$gamlr,select=0)[-1,] )
sgnmrg <- getsign( coef(mrgal)[-1,] )
sgncp <- getsign( coef(cpbest)[-1] )
sgnsnet <- getsign( coef(snet)[-1,] )
sgnscad <- getsign( coef(scad)[-1] ) # 0.6569343

## tracking the weights

getw <- function(fit){
	b <- coef(fit$gamlr,select=0)[-1,]
	gam <- fit$gamlr$gamma
	w <- matrix(1, nrow=nrow(b),ncol=ncol(b))
	if(gam!=0) w[,-1] <- 1/(1+gam*abs(as.matrix(b)[,-ncol(b)]))
	return(w) }

w0 <- getw(gl0)
w2 <- getw(gl2)
w10 <- getw(gl10)
