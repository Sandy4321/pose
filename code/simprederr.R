
## prediction
e0 <- predict(gl0$gamlr,d$x,select=0)-d$y.val
e1 <- predict(gl1$gamlr,d$x,select=0)-d$y.val
e10 <- predict(gl10$gamlr,d$x,select=0)-d$y.val
emrg <- predict(mrgal$gamlr,d$x,select=0)-d$y.val
eoracle <- predict(oracle,as.data.frame(as.matrix(d$x)))-d$y.val
esnet1se <- predict(snet,as.matrix(d$x),which="parms.1se")-d$y.val
esnetmin <- predict(snet,as.matrix(d$x),which="parms.min")-d$y.val

mse0 <- apply(e0^2,2,mean)
mse1 <- apply(e1^2,2,mean)
mse10 <- apply(e10^2,2,mean)
msemrg <- apply(emrg^2,2,mean)
mseoracle <- mean(eoracle^2)
msesnet1se <- mean(esnet1se^2)
msesnetmin <- mean(esnetmin^2)

MSE <- c(oracle=mseoracle,
         snet1se=msesnet1se,
         snetmin=msesnetmin,
         mrg=msemrg[segmrg],gl0=mse0[seg0],
         gl1=mse1[seg1],gl10=mse10[seg10])
write(paste(round(MSE,2),collapse="|"),
      sprintf("results/%s-mse.txt",OUT),append=TRUE)

R2 <- 1-MSE/mean( (d$y.val-mean(d$y.val))^2 )
write(paste(round(R2,2),collapse="|"),
      sprintf("results/%s-r2.txt",OUT),append=TRUE)

MSElong <- c( 
  mrg=fill(msemrg),gl0=fill(mse0),
  gl1=fill(mse1),gl10=fill(mse10))
write(paste(round(MSElong,2),collapse="|"),
      sprintf("results/%s-mselong.txt",OUT),append=TRUE)

R2long <- 1-MSElong/mean( (d$y.val-mean(d$y.val))^2 )
write(paste(round(R2long,2),collapse="|"),
      sprintf("results/%s-r2long.txt",OUT),append=TRUE)

# best of the best
CVbest <- c(0,1,10)[which.min(c(gl0$cvm[gl0$seg.min],gl1$cvm[gl1$seg.min],gl10$cvm[gl10$seg.min]))]
ICbest <- c(0,1,10)[which.min(c(min(AICc(gl0$g)),min(AICc(gl1$g)),min(AICc(gl10$g))))]
MSECVbest <- get(sprintf("mse%d",CVbest))[get(sprintf("seg%d",CVbest))[2]]
MSEICbest <- get(sprintf("mse%d",ICbest))[get(sprintf("seg%d",ICbest))[3]]

bestests <- c(CVbest,MSECVbest,ICbest,MSEICbest)
write(paste(round(bestests,2),collapse="|"),
      sprintf("results/%s-bestgamma.txt",OUT),append=TRUE)

bestests[c(2,4)] <- bestests[c(2,4)]/mean( (d$y.val-mean(d$y.val))^2 )
write(paste(round(bestests,2),collapse="|"),
      sprintf("results/%s-R2bestgam.txt",OUT),append=TRUE)

