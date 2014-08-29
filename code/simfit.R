## fits
id <- 0 # debug
s2n <- 1
rho <- 0.5

source("code/simdata.R")
library(gamlr)

## draw the data
d <- dgp(id=id, s2n=s2n, rho=rho)

## gamma lasso
tgl0 <- system.time(gl0 <- cv.gamlr(d$x, d$y.train))[[3]]
tgl2 <- system.time(gl2 <- cv.gamlr(d$x, d$y.train, gamma=2))[[3]]
tgl10 <- system.time(gl10 <- cv.gamlr(d$x, d$y.train, gamma=10))[[3]]

## marginal adaptive lasso
tmrgal <- system.time({
	wmrg <- 1/abs(cor(as.matrix(d$x),d$y.train))
	mrgal <- cv.gamlr(d$x,d$y.train,varweight=wmrg) })[[3]]

## scad (way slow!!!)
library(ncvreg)
tscad <- system.time(scad <- cv.ncvreg(X=as.matrix(d$x),y=d$y.train,
				lambda.min=0.01,penalty="SCAD")[[3]]
#    user  system elapsed 
# 267.589   0.271 267.785 

## sparsenet 
library(sparsenet)
tsnet <- system.time(
	snet <- cv.sparsenet(as.matrix(d$x),d$y.train,nfolds=5))[[3]]
#   user  system elapsed 
# 11.065   0.203  11.265 

# same order as GL, 
# but GL is parallelizable over both gamma and fold

## Cp optimal L0 pen
fitcp <- function(l){
	lm(d$y.train ~ as.matrix(d$x[,1:l]))}
cpfits <- lapply(1:200,fitcp)

pencp <- function(fit){
	0.5*sum(fit$resid^2) + d$sigma^2*length(coef(fit)) }
cpcosts <- sapply(cpfits,pencp)

sCp <- which.min(cpcosts)
cpbest <- cpfits[[sCp]]









