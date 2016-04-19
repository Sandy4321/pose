## support and sign recovery
getsupport <- function(fit){
  b <- coef(fit$gamlr,select=0)[-1,]
  apply(b,2,function(c) which(c!=0))
}

S0 <- getsupport(gl0)
S1 <- getsupport(gl1)
S10 <- getsupport(gl10)
Smrg <- getsupport(mrgal)
Ssnet1se <- which(coef(snet, which="parms.1se")[-1,]!=0)
Ssnetmin <- which(coef(snet, which="parms.min")[-1,]!=0)

s0 <- sapply(S0,length)
s1 <- sapply(S1,length)
s10 <- sapply(S10,length)
smrg <- sapply(Smrg,length)
ssnet1se <- length(Ssnet1se) 
ssnetmin <- length(Ssnetmin) 

s <- c(cp=sCp,snet1se=ssnet1se,
       snetmin=ssnetmin, mrg=smrg[segmrg],
       gl0=s0[seg0],gl1=s1[seg1],gl10=s10[seg10])
write(paste(round(s,2),collapse="|"),
      sprintf("results/%s-s.txt",OUT),append=TRUE)

slong <- c(mrg=fill(smrg),
           gl0=fill(s0),gl1=fill(s1),gl10=fill(s10))
write(paste(round(slong,2),collapse="|"),
      sprintf("results/%s-slong.txt",OUT),append=TRUE)

fp0 <- sapply(S0,function(set) sum(set>sCp))
fp1 <- sapply(S1,function(set) sum(set>sCp))
fp10 <- sapply(S10,function(set) sum(set>sCp))
fpmrg <- sapply(Smrg,function(set) sum(set>sCp))
fpsnet1se <- sum(Ssnet1se>sCp)
fpsnetmin <- sum(Ssnetmin>sCp)

fdr <- c(cp=0, 
         snet1se=fpsnet1se/ssnet1se,
         snetmin=fpsnetmin/ssnetmin,
         mrg=(fpmrg/smrg)[segmrg],gl0=(fp0/s0)[seg0],
         gl1=(fp1/s1)[seg1],gl10=(fp10/s10)[seg10])
fdr[is.nan(fdr)|is.infinite(fdr)] <- 0
write(paste(round(fdr*100,2),collapse="|"),
      sprintf("results/%s-fdr.txt",OUT),append=TRUE)

fdrlong <- c(mrg=fill(fpmrg/smrg),gl0=fill(fp0/s0),
             gl1=fill(fp1/s1),gl10=fill(fp10/s10))
fdrlong[is.nan(fdrlong)|is.infinite(fdrlong)] <- 0
write(paste(round(fdrlong*100,2),collapse="|"),
      sprintf("results/%s-fdrlong.txt",OUT),append=TRUE)

tp0 <- sapply(S0,function(set) sum(set<=sCp))
tp1 <- sapply(S1,function(set) sum(set<=sCp))
tp10 <- sapply(S10,function(set) sum(set<=sCp))
tpmrg <- sapply(Smrg,function(set) sum(set<=sCp))
tpsnet1se <- sum(Ssnet1se<=sCp)
tpsnetmin <- sum(Ssnetmin<=sCp)

sens <- c(cp=1,
          snet1se=tpsnet1se/sCp,
          snetmin=tpsnetmin/sCp,
          mrg=tpmrg[segmrg]/sCp,gl0=tp0[seg0]/sCp,
          gl1=tp1[seg1]/sCp,gl10=tp10[seg10]/sCp)
write(paste(round(sens*100,2),collapse="|"),
      sprintf("results/%s-sens.txt",OUT),append=TRUE)

senslong <- c(
  mrg=fill(tpmrg/sCp),gl0=fill(tp0/sCp),
  gl1=fill(tp1/sCp),gl10=fill(tp10/sCp))
write(paste(round(senslong*100,2),collapse="|"),
      sprintf("results/%s-senslong.txt",OUT),append=TRUE)

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
sgn1 <-  getsign( coef(gl1$gamlr,select=0)[-1,] )
sgn10 <-  getsign( coef(gl10$gamlr,select=0)[-1,] )
sgnmrg <- getsign( coef(mrgal$gamlr,select=0)[-1,] )
sgncp <- getsign( coef(cpbest)[-1] )
sgnsnet1se <- getsign( coef(snet, which="parms.1se")[-1,] )
sgnsnetmin <- getsign( coef(snet, which="parms.min")[-1,] )

sgn <- c(cp=sgncp,
         snet1se=sgnsnet1se,
         snetmin=sgnsnetmin,
         mrg=sgnmrg[segmrg],gl0=sgn0[seg0],
         gl1=sgn1[seg1],gl10=sgn10[seg10])
write(paste(round(sgn,2),collapse="|"),
      sprintf("results/%s-sgn.txt",OUT),append=TRUE)
