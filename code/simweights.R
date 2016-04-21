wmrg <- matrix(wmrg,nrow=ncol(d$x),ncol=100)
w0 <- matrix(1,nrow=ncol(d$x),ncol=100)

getw <- function(fit){
  b <- coef(fit$gamlr,select=0)[-1,]
  gam <- fit$gamlr$gamma
  w <- matrix(1, nrow=nrow(b),ncol=ncol(b))
  if(gam!=0) w[,-1] <- 1/(1+gam*abs(as.matrix(b)[,-ncol(b)]))
  return(w) }
w1 <- getw(gl1)
w10 <- getw(gl10)

nu <- d$sigma^2/nrow(d$x)

S <- 1:sO
print(sO)
writeE <- function(W,lam,f){  
  #print(W)
  wsnorm <- lam*apply(W[S,],2,
                  function(w) sqrt(sum(w^2)))/sqrt(sO)
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
writeE(w1,gl1$gamlr$lambda,"gl1")
writeE(w10,gl10$gamlr$lambda,"gl10")


