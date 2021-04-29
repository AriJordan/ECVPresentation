library('irlba')

SP <- function(A,K,regularize=TRUE,tau=1){
    n <- nrow(A)
    if(regularize){
        degrees <- colSums(A)
        degrees.reg <- degrees + mean(degrees)*tau
        L <- sqrt(diag(1/degrees.reg))%*%A%*%sqrt(diag(1/degrees.reg))
        SVD <- irlba(L,nu=K,nv=K)
    }else{
        SVD <- irlba(A,nu=K,nv=K)
    }
    km <- kmeans(SVD$v[,1:K],centers=K,nstart=20,iter.max=200)#,algorithm="Lloyd")
    return(list(cluster=km$cluster,loss=km$tot.withinss))
}


SSP <- function(A,K,regularize=TRUE,tau=1){
    n <- nrow(A)
    if(regularize){
        degrees <- colSums(A)
        degrees.reg <- degrees + mean(degrees)*tau
        L <- sqrt(diag(1/degrees.reg))%*%A%*%sqrt(diag(1/degrees.reg))
        #SVD <- irlba(L,nu=K,nv=K)
        SVD <- svd(L,nu=K,nv=K)
        V <- SVD$u[,1:K]
        V.norm <- apply(V,1,function(x)sqrt(sum(x^2)))+0.001
        V.normalized <- diag(1/V.norm)%*%V
    }else{
        SVD <- irlba(A,nu=K,nv=K)
        V <- SVD$v[,1:K]
        V.norm <- apply(V,1,function(x)sqrt(sum(x^2)))
        V.normalized <- diag(1/V.norm)%*%V
    }
    km <- kmeans(V.normalized,centers=K,nstart=20,iter.max=200,algorithm="Lloyd")
    return(list(cluster=km$cluster,loss=km$tot.withinss))
}


Arash.reg.SP <- function(A,K,tau=1,lap=FALSE){
    avg.d <- mean(colSums(A))
    A.tau <- A + tau*avg.d/nrow(A)
    if(!lap){SVD <- irlba(A.tau,nu=K,nv=K)}else{
         d.tau <- colSums(A.tau)
         L.tau <- diag(1/sqrt(d.tau))%*%A.tau%*%diag(1/sqrt(d.tau))
         #SVD <- svd(L.tau,nu=K,nv=K)
         SVD <- irlba(L.tau,nu=K,nv=K)
    }
    km <- kmeans(SVD$v[,1:K],centers=K,nstart=30,iter.max=30)#,algorithm="Lloyd")
    return(list(cluster=km$cluster,loss=km$tot.withinss))
}


Arash.reg.SSP <- function(A,K,tau=1,lap=FALSE){
    avg.d <- mean(colSums(A))
    A.tau <- A + tau*avg.d/nrow(A)
    if(!lap){SVD <- irlba(A.tau,nu=K,nv=K)
         V <- SVD$v[,1:K]
         V.norm <- apply(V,1,function(x)sqrt(sum(x^2)))
         V.normalized <- diag(1/V.norm)%*%V
         }else{
         d.tau <- colSums(A.tau)
         L.tau <- diag(1/sqrt(d.tau))%*%A.tau%*%diag(1/sqrt(d.tau))
         #SVD <- svd(L.tau,nu=K,nv=K)
         SVD <- irlba(L.tau,nu=K,nv=K)
         V <- SVD$v[,1:K]
         V.norm <- apply(V,1,function(x)sqrt(sum(x^2)))
         V.normalized <- diag(1/V.norm)%*%V
    }
    km <- kmeans(V.normalized,centers=K,nstart=50,iter.max=50)#,algorithm="Lloyd")
    return(list(cluster=km$cluster,loss=km$tot.withinss))
}


SBM.estimate <- function(A,m){
    n <- nrow(A)
    K <- length(unique(m))
    B <- matrix(0,K,K)
    M <- matrix(0,n,K)
    for(i in 1:K){
        for(j in i:K){
            if(i!=j){
                B[j,i] <- B[i,j] <- mean(A[which(m==i),which(m==j)])
            }else{
                n.i <- length(which(m==i))
                B[i,i] <- sum(A[which(m==i),which(m==i)])/(n.i^2 - n.i)
            }
        }
    }
    M[matrix(c(1:n,m),ncol=2)] <- 1
    P <- M%*%B%*%t(M)
    return(list(B=B,P=P))
}


SBM.eval <- function(A,B,m,m1){
    hat <- SBM.estimate(A,m1)
    n <- nrow(A)
    K <- length(unique(m))
    M <- matrix(0,n,K)
    M[matrix(c(1:n,m),ncol=2)] <- 1
    P <- M%*%B%*%t(M)
    err <- norm(P-hat$P,"F")^2/n
    return(err)
}


CASP <- function(A,K,regularize=FALSE,X,h,ac=TRUE){
    if(regularize){
        degrees <- colSums(A)
        degrees.reg <- degrees + mean(degrees)
        L <- sqrt(diag(1/degrees.reg))%*%A%*%sqrt(diag(1/degrees.reg))
        if(ac){
        SVD <- svd(L+h*X%*%t(X),nu=K,nv=K)
        }else{
        	SVD <- svd(L%*%L+h*X%*%t(X),nu=K,nv=K)
        }
    }else{
        SVD <- svd(A+h*X%*%t(X),nu=K,nv=K)
    }
    km <- kmeans(SVD$v[,1:K],centers=K,nstart=5,iter.max=30,algorithm="Lloyd")
    return(list(cluster=km$cluster,loss=km$tot.withinss))
}


PairMatch <- function(g1,g2){
    n <- length(g1)
    M1 <- M2 <- matrix(0,nrow=n,ncol=n)
    for(k in 1:n){
        M1[k,g1==g1[k]] <- 1
        M2[k,g2==g2[k]] <- 1
    }
    return(norm(M1-M2,"F")/n)
}



NMI <- function(g1,g2){
    K <- max(length(unique(g1)),length(unique(g2)))
    R <- matrix(0,K,K)
    for(i in 1:K){
        for(j in 1:K){
            R[i,j] <- length(intersect(which(g1==i),which(g2==j)))
        }
    }
    P <- R/length(g1)
    P.row <- rowSums(R+0.0001)/length(g1)
    P.col <- colSums(R+0.0001)/length(g2)
    P.normal <- diag(1/P.row)%*%P%*%diag(1/P.col)
    I <- -sum(P*log(P.normal),na.rm=TRUE)
    H.row <- -sum(P.row*log(P.row))
    H.col <- -sum(P.col*log(P.col))
    result <- -2*I/(H.row+H.col)
    return(result)
}
