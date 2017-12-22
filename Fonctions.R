

permutation_nb <- function(x, nb, n_samples=99) {
  permutation_statistics0 <- matrix(0, n_samples)
  permutation_statistics1 <- matrix(0, n_samples)
  for (i in 1:n_samples) {
    t <-  x[sample(nrow(x), nb, replace=T),]
    m0 <- MASS::glm.nb(total~year*wipo_field, data=t)
    permutation_statistics0[i,] <- m0$aic
    m1 <- update(m0, . ~ . - year:wipo_field)
    permutation_statistics1[i,] <- m1$aic
  }
  result <- list(interaction=permutation_statistics0,
                 noninteraction=permutation_statistics1)
  return(result)
}


permutation_norm <- function(x, n_samples=9999) {
  permutation_statistics0 <- matrix(0, n_samples)
  permutation_statistics1 <- matrix(0, n_samples)
  for (i in 1:n_samples) {
    t <-  x[sample(nrow(x), nrow(x), replace=T),]
    m0 <- glm(total~year*wipo_field, data=t)
    permutation_statistics0[i,] <- m0$aic
    m1 <- update(m0, . ~ . - year:wipo_field)
    permutation_statistics1[i,] <- m1$aic
  }
  result <- list(interaction=permutation_statistics0,
                 noninteraction=permutation_statistics1)
  return(result)
}


empirical_pvalue_right <- function(interaction, noninteraction) {
  n_above <- sum(interaction >= noninteraction)
  n_samples <- length(interaction)
  return((n_above)/(n_samples))
}

empirical_pvalue_left <- function(interaction, noninteraction) {
  n_above <- sum(interaction <= noninteraction)
  n_samples <- length(interaction)
  return((n_above)/(n_samples))
}

empirical_pvalue_twosided <- function(interaction, noninteraction) {
  return(2*min(empirical_pvalue_left(interaction, noninteraction),
               empirical_pvalue_right(interaction, noninteraction)))  
}

permutation_pvalue_right <- function(p) {
  return(empirical_pvalue_right(p$interaction, p$noninteraction))
}

permutation_pvalue_left <- function(p) {
  return(empirical_pvalue_left(p$interaction, p$noninteraction))
}

permutation_pvalue_twosided <- function(p) {
  return(empirical_pvalue_twosided(p$interaction, p$noninteraction))
}




high.indep.test=function(n) {
  cat("High-dim Independence Test","\n")
  N=sum(n)
  r=length(n[,1,1]) #定义r×c×t矩阵大小
  c=length(n[1,,1])
  t=length(n[1,1,])
  ni__=apply(n,1,sum) #开始算各种小数据
  n_j_=apply(n,2,sum)
  n__k=apply(n,3,sum)
  ni=array(ni__,c(r,c,t))
  nj=aperm(array(n_j_,c(c,r,t)),c(2,1,3))
  nk=aperm(array(n__k,c(t,c,r)),c(3,2,1))
  n_jk=aperm(array(apply(n,c(2,3),sum),c(c,t,r)),c(3,1,2))
  ni_k=aperm(array(apply(n,c(1,3),sum),c(r,t,c)),c(1,3,2))
  nij_=array(apply(n,c(1,2),sum),c(r,c,t))
  G1=sum((-2)*n*log(ni*nj*nk/n/N^2)) #开始算各种检验统计量
  G2=sum((-2)*n*log(ni*n_jk/n/N))
  G3=sum((-2)*n*log(nj*ni_k/n/N))
  G4=sum((-2)*n*log(nk*nij_/n/N))
  G5=sum((-2)*n*log(ni_k*nij_/ni/n))
  G6=sum((-2)*n*log(nij_*n_jk/nj/n))
  G7=sum((-2)*n*log(ni_k*n_jk/nk/n))
  G=matrix(c(G1,G2,G3,G4,G5,G6,G7),nrow=1,ncol=7,dimnames=list(c("G^2"),
                                                               c("(A,B,C)","(A,BC)","(B,AC)","(C,AB)","(AB,AC)","(BA,BC)","(CA,CB)")))
  cat("检验统计量如下：","\n")
  print(G) #输出7个检验统计量
  P1=1-pchisq(G1,r*c*t-r-c-t+2)
  P2=1-pchisq(G2,(r-1)*(c*t-1))
  P3=1-pchisq(G3,(c-1)*(r*t-1))
  P4=1-pchisq(G4,(t-1)*(r*c-1))
  P5=1-pchisq(G5,r*(c-1)*(t-1))
  P6=1-pchisq(G6,c*(r-1)*(t-1))
  P7=1-pchisq(G7,t*(r-1)*(c-1))
  P=matrix(c(P1,P2,P3,P4,P5,P6,P7),nrow=1,ncol=7,dimnames=list(c("P-value"),
                                                               c("(A,B,C)","(A,BC)","(B,AC)","(C,AB)","(AB,AC)","(BA,BC)","(CA,CB)")))
  print(P) #输出检验统计量对应的7个p值
}