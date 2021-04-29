library(pROC)
library(softImpute)
source("SimonFunk.R")


LeiSim1 <- function(r,K,n,n1,n2=NULL,dg=3){



B0 <- matrix(1,K,K)
diag(B0) <- dg

if(is.null(n2)){
n2 <- (n-n1)/(K-1)

n.seq <- c(0,n1,rep(n2,K-1))
}else{
n3 <- (n-n1-n2)/(K-2)

n.seq <- c(0,n1,n2,rep(n3,K-2))
}
B <- r*B0

Theta <- matrix(0,nrow=n,ncol=K)
for(k in 2:(K+1)){
    Theta[((sum(n.seq[1:(k-1)])+1):(sum(n.seq[1:k]))),(k-1)] <- 1
}



P <- Theta%*%B%*%t(Theta)
A <- matrix(0,n,n)
for(i in 1:(n-1)){
    for(j in (i+1):n){
        A[i,j] <- rbinom(1,1,prob=P[i,j])
    }
}

A <- A + t(A)
return(A)

}





library(Matrix)

NCV <- function(A,cv,K){
n <- nrow(A)
sample.index <- sample.int(n)
max.fold.num <- ceiling(n/cv)
fold.index <- rep(1:cv,each=max.fold.num)[1:n]

cv.index <- fold.index[sample.index]
l2 <- log.like <- rep(0,cv)

for(k in 1:cv){
    holdout.index <- which(cv.index==k)
    train.index <- which(cv.index!=k)
    tmp.eval <- cv.evaluate(A,train.index,holdout.index,K)
    log.like[k] <- tmp.eval$loglike

    l2[k] <- tmp.eval$l2
}
return(list(log.like=log.like,l2=l2))
}



NCV.select <- function(A,max.K,cv=3,fast=FALSE,soft=TRUE,dc.est=1){
dc.avg.auc <- dc.avg.impute <- dc.avg.se <- dc.avg.log <- avg.auc <- avg.impute <- avg.se <- avg.log <- rep(0,max.K)
dc.avg.se[1] <- dc.avg.log[1] <- avg.impute[1] <- avg.se[1] <- avg.log[1] <- Inf
avg.auc[1] <- 0
dc.dev.mat <- dc.l2.mat <- sbm.dev.mat <- sbm.l2.mat <- matrix(0,cv,max.K)
n <- nrow(A)
sample.index <- sample.int(n)
max.fold.num <- ceiling(n/cv)
fold.index <- rep(1:cv,each=max.fold.num)[1:n]

cv.index <- fold.index[sample.index]
for(KK in 1:max.K){
    dc.impute.err.seq <- dc.auc.seq <- dc.no.edge.seq <- impute.err.seq <- auc.seq <- no.edge.seq <- dc.l2 <- l2 <- dc.log.like <- log.like <- rep(0,cv)
    print(paste("Start",KK))
    for(k in 1:cv){
        holdout.index <- which(cv.index==k)
        train.index <- which(cv.index!=k)
        if(fast){
            tmp.eval <- cv.evaluate.fast.all(A,train.index,holdout.index,KK)
            tmp.eval.dc <- cv.evaluate.DC.fast.all(A,train.index,holdout.index,KK,dc.est=dc.est)
        }else{
            tmp.eval <- cv.evaluate(A,train.index,holdout.index,KK)
            tmp.eval.dc <- cv.evaluate.DC(A,train.index,holdout.index,KK,dc.est=dc.est)
        }
        log.like[k] <- tmp.eval$loglike

        sbm.l2.mat[k,KK] <- l2[k] <- tmp.eval$l2
        no.edge.seq[k] <- tmp.eval$no.edge
        auc.seq[k] <- tmp.eval$AUC
        impute.err.seq[k] <- tmp.eval$impute.err
        sbm.dev.mat[k,KK] <- log.like[k] <- tmp.eval$loglike

        dc.l2.mat[k,KK] <- dc.l2[k] <- tmp.eval.dc$l2
        dc.dev.mat[k,KK] <- dc.log.like[k] <- tmp.eval.dc$loglike
    }

    avg.se[KK] <- mean(l2)
    avg.log[KK] <- mean(log.like)
    avg.auc[KK] <- mean(auc.seq)
    avg.impute[KK] <- mean(impute.err.seq)
    dc.avg.se[KK] <- mean(dc.l2)
    dc.avg.log[KK] <- mean(dc.log.like)


    print(paste("Finish ",KK,"....",sep=""))
}
   return(list(loglike=avg.log,SE=avg.se,impute.err=avg.impute,auc=avg.auc,no.edge.seq=no.edge.seq,dc.loglike=dc.avg.log,dc.SE=dc.avg.se,sbm.l2.mat=sbm.l2.mat,sbm.dev.mat=sbm.dev.mat,dc.l2.mat=dc.l2.mat,dc.dev.mat=dc.dev.mat))
}




cv.evaluate <- function(A,train.index,holdout.index,K){
    n <- nrow(A)
    A.new <- A[c(train.index,holdout.index),c(train.index,holdout.index)]
    n.holdout <- length(holdout.index)
    n.train <- n-n.holdout
    A1 <- A.new[1:n.train,]
    A1.svd <- irlba(A1+0.001,nu=K,nv=K)
    #A1.svd <- irlba(A1,nu=K,nv=K)

    if(K==1){
      A0 <- A1[1:n.train,1:n.train]
      pb <- sum(A0)/n.train^2
      if(pb < 1e-6) pb <- 1e-6
      if(pb > 1- 1e-6) pb <- 1-1e-6
      A.2 <- A.new[(n.train+1):n,(n.train+1):n]
      sum.index <- lower.tri(A.2)
      loglike <- -sum(A.2[sum.index]*log(pb)) - sum((1-A.2[sum.index])*log(1-pb))
      l2 <- sum((A.2[sum.index]-pb)^2)
      return(list(loglike=loglike,l2=l2,no.edge=NA,impute.err=NA,AUC=NA))
    } 

    V <- A1.svd$v
    km <- kmeans(V,centers=K,nstart=30,iter.max=30)

    degrees <- colSums(A1)
    no.edge <- sum(degrees==0)


    B <- matrix(0,K,K)
    for(i in 1:(K-1)){
        for(j in (i+1):K){
            N.1i <- intersect(1:n.train,which(km$cluster==i))
            N.2i <- intersect((n.train+1):n,which(km$cluster==i))
            N.1j <- intersect(1:n.train,which(km$cluster==j))
            N.2j <- intersect((n.train+1):n,which(km$cluster==j))
            B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j]) + sum(A.new[N.1j,N.2i])+1)/(length(N.1i)*length(N.1j)+length(N.1j)*length(N.2i)+length(N.1i)*length(N.2j)+1)
            #B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j])+1)/(length(N.1i)*length(N.1j)+length(N.1i)*length(N.2j)+1)
        }
    }
    B <- B+t(B)
    Theta <- matrix(0,n,K)
    for(i in 1:K){
            N.1i <- intersect(1:n.train,which(km$cluster==i))
            N.2i <- intersect((n.train+1):n,which(km$cluster==i))
            B[i,i] <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(length(N.1i)*(length(N.1i)-1)/2+length(N.1i)*length(N.2i)+1)
            Theta[which(km$cluster==i),i] <- 1

    }
    P.hat.holdout <- Theta[(n.train+1):n,]%*%B%*%t(Theta[(n.train+1):n,])
    P.hat.holdout[P.hat.holdout<1e-6] <- 1e-6
    P.hat.holdout[P.hat.holdout>(1-1e-6)] <- 1-1e-6
    A.2 <- A.new[(n.train+1):n,(n.train+1):n]
    sum.index <- lower.tri(A.2)
    loglike <- -sum(A.2[sum.index]*log(P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-P.hat.holdout[sum.index]))

    l2 <- sum((A.2[sum.index]-P.hat.holdout[sum.index])^2)
    return(list(loglike=loglike,l2=l2,no.edge=no.edge,impute.err=NA,AUC=NA))



}

cv.approx.evaluate <- function(A,train.index,holdout.index,K,soft=TRUE){
    n <- nrow(A)
    A.new <- A[c(train.index,holdout.index),c(train.index,holdout.index)]
    n.holdout <- length(holdout.index)
    n.train <- n-n.holdout
    A.miss <- A.new
    A.miss <- A.miss + 0.001
    A.miss[(n.train+1):n,(n.train+1):n] <- NA
    diag(A.miss) <- 0
    non.miss <- which(!is.na(A.miss))
    #A.miss[non.miss] <- A.miss[non.miss] + 0.001
    A1 <- A.new[1:n.train,]
    if(soft){
    SVD <- softImpute(x=as.matrix(A.miss),rank.max=K,maxit=150,type="svd")
    if(K==1) {
        A.new.approx <- list(A=matrix(SVD$u,ncol=K)%*%t(matrix(SVD$v,ncol=K))*as.numeric(SVD$d),SVD=SVD)
    }else{
        A.new.approx <- list(A=SVD$u%*%diag(SVD$d)%*%t(SVD$v),SVD=SVD)
        }
    }else{
    A.new.approx <- iter.SVD.core(A=A.miss,K=K,max.iter=300)#SF.SVD(A=A.miss,K=K,max.iter=500,learning.rate=0.02,verbose=FALSE)
    }
    V <- A.new.approx$SVD$v
    km <- kmeans(V,centers=K,nstart=5)
    degrees <- colSums(A1)
    no.edge <- sum(degrees==0)


    Omega <- which(is.na(A.miss))
    response <- A.new[Omega]
    predictors <- A.new.approx$A[Omega]
    impute.err <- sum(response-predictors)^2
    auc <- as.numeric(roc(response=response,predictor=predictors)$auc)


    if(K==1){
      A0 <- A1[1:n.train,1:n.train]
      pb <- sum(A0)/n.train^2
      if(pb < 1e-6) pb <- 1e-6
      if(pb > 1- 1e-6) pb <- 1-1e-6
      A.2 <- A.new[(n.train+1):n,(n.train+1):n]
      sum.index <- lower.tri(A.2)
      loglike <- -sum(A.2[sum.index]*log(pb)) - sum((1-A.2[sum.index])*log(1-pb))
      l2 <- sum((A.2[sum.index]-pb)^2)
      return(list(loglike=loglike,l2=l2,no.edge=NA,impute.err=impute.err,AUC=auc))
    } 

    B <- matrix(0,K,K)
    for(i in 1:(K-1)){
        for(j in (i+1):K){
            N.1i <- intersect(1:n.train,which(km$cluster==i))
            N.2i <- intersect((n.train+1):n,which(km$cluster==i))
            N.1j <- intersect(1:n.train,which(km$cluster==j))
            N.2j <- intersect((n.train+1):n,which(km$cluster==j))
            B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j]) + sum(A.new[N.1j,N.2i])+1)/(length(N.1i)*length(N.1j)+length(N.1j)*length(N.2i)+length(N.1i)*length(N.2j)+1)
            #B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j])+1)/(length(N.1i)*length(N.1j)+length(N.1i)*length(N.2j)+1)
        }
    }
    B <- B+t(B)
    Theta <- matrix(0,n,K)
    for(i in 1:K){
            N.1i <- intersect(1:n.train,which(km$cluster==i))
            N.2i <- intersect((n.train+1):n,which(km$cluster==i))
            B[i,i] <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(length(N.1i)*(length(N.1i)-1)/2+length(N.1i)*length(N.2i)+1)
            Theta[which(km$cluster==i),i] <- 1

    }
    P.hat.holdout <- Theta[(n.train+1):n,]%*%B%*%t(Theta[(n.train+1):n,])
    P.hat.holdout[P.hat.holdout<1e-6] <- 1e-6
    P.hat.holdout[P.hat.holdout>(1-1e-6)] <- 1-1e-6
    A.2 <- A.new[(n.train+1):n,(n.train+1):n]
    sum.index <- lower.tri(A.2)
    loglike <- -sum(A.2[sum.index]*log(P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-P.hat.holdout[sum.index]))

    l2 <- sum((A.2[sum.index]-P.hat.holdout[sum.index])^2)
    return(list(loglike=loglike,l2=l2,no.edge=no.edge,impute.err=impute.err,AUC=auc))


}



cv.impute.evaluate <- function(A,train.index,holdout.index,K){
    n <- nrow(A)
    A.new <- A[c(train.index,holdout.index),c(train.index,holdout.index)]
    n.holdout <- length(holdout.index)
    n.train <- n-n.holdout
    A.miss <- A.new
    A.miss[(n.train+1):n,(n.train+1):n] <- NA
    diag(A.miss) <- 0
    A1 <- A.new[1:n.train,]
    Omega <- which(is.na(A.miss))
    A.new.approx <- iter.SVD.core(A=A.miss,K=K)$A
    impute.err <- sum((A.new.approx[Omega]-A.new[Omega])^2)
    response <- A.new[Omega]
    predictors <- A.new.approx[Omega]
    tmp.roc <- roc(response=response,predictor=predictors)
    roc.auc <- as.numeric(tmp.roc$auc)


    degrees <- colSums(A1)
    no.edge <- sum(degrees==0)


    return(list(impute.err=impute.err,auc=roc.auc,no.edge=no.edge))
}


NCV.impute.select <- function(A,max.K,cv=3){
avg.se <- avg.log <- rep(0,max.K)
avg.se[1] <- avg.log[1] <- Inf
n <- nrow(A)
sample.index <- sample.int(n)
max.fold.num <- ceiling(n/cv)
fold.index <- rep(1:cv,each=max.fold.num)[1:n]

cv.index <- fold.index[sample.index]
for(KK in 2:max.K){
    no.edge.seq <- impute.seq <- auc.seq <- rep(0,cv)

    for(k in 1:cv){
        holdout.index <- which(cv.index==k)
        train.index <- which(cv.index!=k)
        tmp.eval <- cv.impute.evaluate(A,train.index,holdout.index,KK)

        impute.seq[k] <- tmp.eval$impute.err

        auc.seq[k] <- tmp.eval$auc
        no.edge.seq[k] <- tmp.eval$no.edge
    }

    print(paste("Finish ",KK,"....",sep=""))
}
   return(list(impute.select=which.min(impute.seq),auc.select=which.max(auc.seq[-1])+1,no.edge.seq=no.edge.seq))
}



cv.evaluate.DC <- function(A,train.index,holdout.index,K,dc.est=1){
    n <- nrow(A)
    A.new <- A[c(train.index,holdout.index),c(train.index,holdout.index)]
    n.holdout <- length(holdout.index)
    n.train <- n-n.holdout
    A1 <- A.new[1:n.train,]
    A1.svd <- irlba(A1,nu=K,nv=K)
    V <- A1.svd$v[,1:K]
    if(K==1) {V.norms <- abs(V)}else{
    V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
    }
    iso.index <- which(V.norms==0)
    Psi <- V.norms
    Psi <- Psi / max(V.norms)
    Psi.outer <- outer(Psi,Psi)
    inv.V.norms <- 1/V.norms
    inv.V.norms[iso.index] <- 1
    V.normalized <- diag(inv.V.norms)%*%V

    if(K==1){
            N.1i <- 1:n.train
            N.2i <- (n.train+1):n
            pb <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(sum(Psi.outer[N.1i,N.1i])/2 + sum(Psi.outer[N.1i,N.2i]) - sum(diag(Psi.outer))+1)


    P.hat.holdout <-  diag(Psi[(n.train+1):n])%*%matrix(1,ncol=(n-n.train),nrow=(n-n.train))%*%diag(Psi[(n.train+1):n])*pb
    P.hat.holdout[P.hat.holdout<1e-6] <- 1e-6
    P.hat.holdout[P.hat.holdout>(1-1e-6)] <- 1-1e-6
    A.2 <- A.new[(n.train+1):n,(n.train+1):n]
    sum.index <- lower.tri(A.2)
    loglike <- -sum(A.2[sum.index]*log(P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-P.hat.holdout[sum.index]))

    l2 <- sum((A.2[sum.index]-P.hat.holdout[sum.index])^2)
    return(list(loglike=loglike,l2=l2,no.edge=NA,impute.err=NA,AUC=NA))
    }


    km <- kmeans(V.normalized,centers=K,nstart=30,iter.max=30)

    degrees <- colSums(A1)
    no.edge <- sum(degrees==0)

    B <- matrix(0,K,K)
    
    for(i in 1:(K-1)){
        for(j in (i+1):K){
            N.1i <- intersect(1:n.train,which(km$cluster==i))
            N.2i <- intersect((n.train+1):n,which(km$cluster==i))
            N.1j <- intersect(1:n.train,which(km$cluster==j))
            N.2j <- intersect((n.train+1):n,which(km$cluster==j))
            B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j]) + sum(A.new[N.1j,N.2i])+1)/(sum(Psi.outer[N.1i,N.1j]) + sum(Psi.outer[N.1i,N.2j]) + sum(Psi.outer[N.1j,N.2i])+1)
            #B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j])+1)/(length(N.1i)*length(N.1j)+length(N.1i)*length(N.2j)+1)
        }
    }
    B <- B+t(B)
    Theta <- matrix(0,n,K)
    for(i in 1:K){
            N.1i <- intersect(1:n.train,which(km$cluster==i))
            N.2i <- intersect((n.train+1):n,which(km$cluster==i))
            B[i,i] <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(sum(Psi.outer[N.1i,N.1i])/2 + sum(Psi.outer[N.1i,N.2i]) - sum(diag(Psi.outer))+1)
            Theta[which(km$cluster==i),i] <- 1

    }
    tmp.imt.mat <- Theta[(n.train+1):n,]*Psi[(n.train+1):n]
    P.hat.holdout <-  tmp.imt.mat%*%B%*%t(tmp.imt.mat)
    #P.hat.holdout <-  diag(Psi[(n.train+1):n])%*%Theta[(n.train+1):n,]%*%B%*%t(Theta[(n.train+1):n,])%*%diag(Psi[(n.train+1):n])
    if(dc.est==2){
            B <- matrix(0,K,K)
            Theta <- matrix(0,n,K)
            A.new.na <- A.new
            A.new.na[(n.train+1):n,(n.train+1):n] <- NA
            for(i in 1:K){
                for(j in 1:K){
                N.i <- which(km$cluster==i)
                N.j <- which(km$cluster==j)
                B[i,j] <- sum(A.new.na[N.i,N.j],na.rm=TRUE)+0.01
                }
                Theta[N.i,i] <- 1
            }
            partial.d <- colSums(A.new.na,na.rm=TRUE)
            partial.gd <- colSums(B)
            phi <- rep(0,n)
            B.g <- Theta%*%partial.gd
            phi <- as.numeric(partial.d/B.g)
            P.hat <- diag(phi)%*%Theta%*%B%*%t(Theta)%*%diag(phi)
            diag(P.hat) <- 0
            P.hat.holdout <- P.hat[(n.train+1):n,(n.train+1):n]
    }
    P.hat.holdout[P.hat.holdout<1e-6] <- 1e-6
    P.hat.holdout[P.hat.holdout>(1-1e-6)] <- 1-1e-6
    A.2 <- A.new[(n.train+1):n,(n.train+1):n]
    sum.index <- lower.tri(A.2)
    loglike <- -sum(A.2[sum.index]*log(P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-P.hat.holdout[sum.index]))

    l2 <- sum((A.2[sum.index]-P.hat.holdout[sum.index])^2)
    return(list(loglike=loglike,l2=l2,no.edge=no.edge,impute.err=NA,AUC=NA))

}






#### Below are code provided by Kehui and Jing

Generate.B <- function(K, low.limit = 0.1, upp.limit = 0.9) {
  B <- matrix(low.limit+ (upp.limit-low.limit)*runif(K^2), ncol = K)
  tB <- t(B)
  B[lower.tri(B)] <- tB[lower.tri(tB)]
  return(B)
}

Generate.clust.size <- function(n, K, prob = NULL){
  if (is.null(prob)) {
    prob <- rep(1/K, K)
  }
  clust.size <- rmultinom(1, n, prob)
  return(clust.size)
}

Generate.psi <- function(n, K, clust.size, DCBM, low.limit = 0.2, upp.limit = 1){
  if (DCBM) {
    psi <- low.limit + (upp.limit-low.limit)* runif(n)
    for (k in 1:K){
      if (k==1){id1 <-1} else {id1 <- sum(clust.size[1:(k-1)])+1}
      id2 <- sum(clust.size[1:k])
      psi[id1:id2] <- psi[id1:id2] / max(psi[id1:id2])
    }
  } else {
    psi <- rep(1, n)
  }
  return(psi)
}

Generate.theta <- function(n, K, clust.size){
  theta <- matrix(0, n, K)
  for (k in 1:K){
      if (k==1){id1 <-1} else {id1 <- sum(clust.size[1:(k-1)])+1}
      id2 <- sum(clust.size[1:k])
      theta[id1:id2, k] <- 1
  }
  return(theta)
}



#### Below are the functions to generate DCSBM according to Chen & Lei's setting

Generate.A.DCSBM <- function(Theta,B,Psi){
    P <- diag(Psi)%*%Theta%*%B%*%t(Theta)%*%diag(Psi)
    n <- nrow(Theta)

    A <- matrix(0,n,n)
    upper.index <- which(upper.tri(P))
    upper.length <- length(upper.index)
    upper.edge.index <- which(runif(upper.length)<P[upper.index])
    A[upper.index[upper.edge.index]] <- 1
    A <- A+t(A)
    return(A)

}


Lei.DCSBM <- function(n,K){
    B <- matrix(0.1,K,K)
    diag(B) <- rep(0.25,K)
    cluster.size <- Generate.clust.size(n,K)
    Theta <- Generate.theta(n,K,cluster.size)
    Psi <- Generate.psi(n,K,cluster.size,DCBM=TRUE,low.limit=0.2,upp.limit=1)
    A <- Generate.A.DCSBM(Theta,B,Psi)
    return(A)
}



#Theta1 <- matrix(0,600,3)
#for(k in 1:600){
#    Theta1[k,clusters[k]] <- 1
#}


#P1 <- diag(psi.hat.final)%*%Theta1%*%B.Lei%*%t(Theta1)%*%diag(psi.hat.final)

#P2 <- diag(Psi)%*%Theta%*%B%*%t(Theta)%*%diag(Psi)






cv.evaluate.fast.all <- function(A,train.index,holdout.index,K){
    n <- nrow(A)
    A.new <- A[c(train.index,holdout.index),c(train.index,holdout.index)]
    n.holdout <- length(holdout.index)
    n.train <- n-n.holdout
    A1 <- A.new[1:n.train,]
    A1.svd <- irlba(A1+0.001,nu=K,nv=K)
    #A1.svd <- irlba(A1,nu=K,nv=K)

    if(K==1){
      A0 <- A1[1:n.train,1:n.train]
      pb <- sum(A0)/n.train^2
      if(pb < 1e-6) pb <- 1e-6
      if(pb > 1- 1e-6) pb <- 1-1e-6
      A.2 <- A.new[(n.train+1):n,(n.train+1):n]
      sum.index <- lower.tri(A.2)
      loglike <- -sum(A.2[sum.index]*log(pb)) - sum((1-A.2[sum.index])*log(1-pb))
      l2 <- sum((A.2[sum.index]-pb)^2)
      return(list(loglike=loglike,l2=l2,no.edge=NA,impute.err=NA,AUC=NA))
    } 

    V <- A1.svd$v
    km <- kmeans(V,centers=K,nstart=30,iter.max=30)

    degrees <- colSums(A1)
    no.edge <- sum(degrees==0)


    B <- matrix(0,K,K)
    for(i in 1:(K-1)){
        for(j in (i+1):K){
            N.1i <- intersect(1:n.train,which(km$cluster==i))
            N.2i <- intersect((n.train+1):n,which(km$cluster==i))
            N.1j <- intersect(1:n.train,which(km$cluster==j))
            N.2j <- intersect((n.train+1):n,which(km$cluster==j))
            B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j]) + sum(A.new[N.1j,N.2i])+1)/(length(N.1i)*length(N.1j)+length(N.1j)*length(N.2i)+length(N.1i)*length(N.2j)+1)
            #B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j])+1)/(length(N.1i)*length(N.1j)+length(N.1i)*length(N.2j)+1)
        }
    }
    B <- B+t(B)
    Theta <- matrix(0,n,K)
    for(i in 1:K){
            N.1i <- intersect(1:n.train,which(km$cluster==i))
            N.2i <- intersect((n.train+1):n,which(km$cluster==i))
            B[i,i] <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(length(N.1i)*(length(N.1i)-1)/2+length(N.1i)*length(N.2i)+1)
            Theta[which(km$cluster==i),i] <- 1

    }
    P.hat.holdout <- Theta[(n.train+1):n,]%*%B%*%t(Theta[(n.train+1):n,])
    P.hat.holdout[P.hat.holdout<1e-6] <- 1e-6
    P.hat.holdout[P.hat.holdout>(1-1e-6)] <- 1-1e-6
    A.2 <- A.new[(n.train+1):n,(n.train+1):n]
    sum.index <- lower.tri(A.2)
    loglike <- -sum(A.2[sum.index]*log(P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-P.hat.holdout[sum.index]))

    l2 <- sum((A.2[sum.index]-P.hat.holdout[sum.index])^2)
    return(list(loglike=loglike,l2=l2,no.edge=no.edge,impute.err=NA,AUC=NA))



}
