## fits

source("code/simdata.R")
library(gamlr)

## draw the data
n <- 1e3
d <- dgp(id=0, n=n)

## gamma lasso
gl0 <- cv.gamlr(d$x, d$y.train, verb=1)
system.time(gl2 <- cv.gamlr(d$x, d$y.train, gamma=2))
#  user  system elapsed 
# 2.066   0.035   2.101 
# > sum(coef(gl2)!=0)
# [1] 65

gl10 <- cv.gamlr(d$x, d$y.train, gamma=10, verb=1)

## marginal adaptive lasso
wmrg <- 1/abs(cor(as.matrix(d$x),d$y.train))
mrgal <- cv.gamlr(d$x,d$y.train,varweight=wmrg)

## scad (way too slow!!!)
library(ncvreg)
system.time(scad <- cv.ncvreg(X=as.matrix(d$x),y=d$y.train,
				lambda.min=0.01,penalty="SCAD"))
#    user  system elapsed 
# 267.589   0.271 267.785 

## sparsenet 
library(sparsenet)
system.time(snet <- cv.sparsenet(as.matrix(d$x),d$y.train,nfolds=5))
#   user  system elapsed 
# 11.065   0.203  11.265 

## Cp optimal L0 pen
fitcp <- function(l){
	lm(d$y.train ~ as.matrix(d$x[,1:l]))}
cpfits <- lapply(1:200,fitcp)

pencp <- function(fit){
	0.5*sum(fit$resid^2) + d$sigma^2*length(coef(fit)) }
cpcosts <- sapply(cpfits,pencp)

sCp <- which.min(cpcosts)
cpbest <- cpfits[[sCp]]









