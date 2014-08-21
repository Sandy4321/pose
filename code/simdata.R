## sim X

# use rho of 0, .5, .9, n=1000 and p=1000.  
# Gammas of 0, 2, 10 in the GL, with selection
# via BIC,AIC,AICc,CVmin,CV1se.

# Compare agains adaptive lasso with marginal weights,
# and scad via the ncvreg package in R.

# Track: Values of L, w.min, ||w.s||, along lambda.
# Track sigma and the Cp selected set (and its size),  
# whether the Cp weight inequality holds.

# For all, report prediction error on left-out true y,
# and prediction error on Cp rule, and selection error on Cp.

dgp <- function(id, n, p=n, 
		s2n=1, rho=0.9, decay=50){

	## fixed 
	xvar <- matrix(ncol=p,nrow=p)
	for(i in 1:p) 
		for(j in i:p) 
			xvar[i,j] <- 0.9^{abs(i-j)}
	C = chol(xvar)
	beta <- matrix( (-1)^(1:p)*exp(-(1:p)/decay) )

	## random
	require(Matrix)
	set.seed(id*2739)
	x <- rnorm(p*n)
	z <- rbinom(n*p, size=1, prob=.5)
	x <- matrix(x, nrow=n)%*%C
	x <- Matrix(x*z, sparse=TRUE)
	mu = as.vector(as.matrix(x%*%beta))
	sigma <- sd(mu)/s2n
	y <- c(mu,mu) + rnorm(2*n,sd=sigma)
	list(x=x,mu=mu,sigma=sigma,
		y.train=y[1:n],y.validate=y[n+1:n]) }













