
## coefion
pad0 <- function(v) c(v, rep(0,1001-length(v)))

betaO <- pad0(coef(cpbest))
este0 <- coef(gl0$gamlr,d$x,select=0)-betaO
este1 <- coef(gl1$gamlr,d$x,select=0)-betaO
este10 <- coef(gl10$gamlr,d$x,select=0)-betaO
estemrg <- coef(mrgal$gamlr,d$x,select=0)-betaO
estecp <- pad0(coef(cpbest,as.data.frame(as.matrix(d$x))))-d$beta
estesnet1se <- coef(snet,as.matrix(d$x),which="parms.1se")-betaO
estesnetmin <- coef(snet,as.matrix(d$x),which="parms.min")-betaO

estmse0 <- apply(este0^2,2,mean)
estmse1 <- apply(este1^2,2,mean)
estmse10 <- apply(este10^2,2,mean)
estmsemrg <- apply(estemrg^2,2,mean)
estmsecp <- mean(estecp^2)
estmsesnet1se <- mean(estesnet1se^2)
estmsesnetmin <- mean(estesnetmin^2)

estMSE <- c(cp=estmsecp,
         snet1se=estmsesnet1se,
         snetmin=estmsesnetmin,
         mrg=msemrg[segmrg],gl0=estmse0[seg0],
         gl1=estmse1[seg1],gl10=estmse10[seg10])
write(paste(round(estMSE,2),collapse="|"),
      sprintf("results/%s-estmse.txt",OUT),append=TRUE)
