## extract sim var
args <- commandArgs(TRUE)
id <- as.integer(args[1])
rho <- as.numeric(args[2])
s2n <- as.numeric(args[3])
decay <- as.numeric(args[4])
nobs <- 1000
binary <- TRUE

OUT=paste("sim-binary-rho",rho,"-s2n",s2n,"-decay",decay,collapse='',sep='')
print(sessionInfo())

source("code/simdata.R")
library(gamlr)
source("code/simfit.R")
print(id)

fill <- function(v){ c(v, rep(tail(v,1),100-length(v)) )}

## basic properties
write(d$sigma, sprintf("results/%s-sigma.txt",OUT),append=TRUE)
times <- paste(round(c(tgl0,tgl1,tgl10,tmrgal,tsnet),1),collapse="|")
write(times, sprintf("results/%s-times.txt",OUT),append=TRUE)
write(paste(gl0$gamlr$lambda,collapse="|"),
        sprintf("results/%s-lambda.txt",OUT),append=TRUE)

#model selection
seg0 <- c(onese=gl0$seg.1se,min=gl0$seg.min,
          aicc=which.min(AICc(gl0$gamlr)),
          aic=which.min(AIC(gl0$gamlr)),
          bic=which.min(BIC(gl0$gamlr)))
seg1 <- c(onese=gl1$seg.1se,min=gl1$seg.min,
          aicc=which.min(AICc(gl1$gamlr)),
          aic=which.min(AIC(gl1$gamlr)),
          bic=which.min(BIC(gl1$gamlr)))
seg10 <- c(onese=gl10$seg.1se,min=gl10$seg.min,
           aicc=which.min(AICc(gl10$gamlr)),
           aic=which.min(AIC(gl10$gamlr)),
           bic=which.min(BIC(gl10$gamlr)))
segmrg <- c(onese=mrgal$seg.1se,min=mrgal$seg.min,
            aicc=which.min(AICc(mrgal$gamlr)),
            aic=which.min(AIC(mrgal$gamlr)),
            bic=which.min(BIC(mrgal$gamlr)))
write(paste(c(segmrg,seg0,seg1,seg10),collapse="|"),
      sprintf("results/%s-seg.txt",OUT),append=TRUE)

# errors
print("prediction error")
source("code/simprederr.R")

print("sensitivity and specificity")
source("code/simsens.R")

print("sensitivity and specificity")
source("code/simestimerr.R")

## tracking the weights
print("weight stats")
source("code/simweights.R")



