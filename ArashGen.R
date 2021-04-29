## lambda: average degree
## n: network size
## beta: in-out ratio
## w:block specific degree
## Pi: block size proportion
## rho: 1-rho is hub-proportion. If rho >0, it indicates that we use degree corrected model.
## simple: if rho >0, then fix the degree of hubs and normal nodes by 1 and 0.2. If simple == FALSE and rho >0, go to more general degree distribution
## power: if TRUE, generate degrees according to exponential with (3*log(10)) as the rate, truncated between [0.2,1]. This is approximatly power law. If FALSE, generate uniform degrees over [0.2,1].

library(poweRlaw)


ArashSBM <- function(lambda,n,beta=0,K=3,w=rep(1,K),Pi=rep(1,K)/K,rho=0,simple=TRUE,power=TRUE,alpha=5){
    P0 <- diag(w)
    if(beta > 0){
        P0 <- matrix(1,K,K)
        diag(P0) <- w/beta
    }
    Pi.vec <- matrix(Pi,ncol=1)
    P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(rho*0.2+(1-rho))^2)
    if((rho >0) && (!simple) && (!power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(0.6)^2)   }
    if((rho >0) && (!simple) && (power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*((1.285)^2))   }

    M <- matrix(0,n,K)
    membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)
    for(i in 1:n){
        M[i,membership[i]] <- 1
    }
    A.bar <- M%*%P%*%t(M)
    node.degree <- rep(1,n)
    if(rho>0){
    if(simple){
        node.degree[runif(n)<rho] <- 0.2
    }else{
        if(power==FALSE){
            node.degree <- runif(n)*0.8 + 0.2
        }else{
            MM <- ceiling(n/300)
            degree.seed <- rplcon(300,1,alpha) ### sample fixed number of values from power law and then randomly assign to be degrees. Easier to control noises across different sample sizes.
            node.degree <- sample(rep(degree.seed,MM)[1:n],size=n,replace=FALSE)
        }
    }}
    DD <- diag(node.degree)
    A.bar <- DD%*%A.bar%*%DD
    diag(A.bar) <- 0
    #avg.d <- mean(colSums(A.bar))
    #A.bar <- A.bar*lambda/avg.d
    upper.index <- which(upper.tri(A.bar))
    upper.p <- A.bar[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    return(A)
}



ArashSBM.with.label <- function(lambda,n,beta=0,K=3,w=rep(1,K),Pi=rep(1,K)/K,rho=0,simple=TRUE,power=TRUE,alpha=5){
    P0 <- diag(w)
    if(beta > 0){
        P0 <- matrix(1,K,K)
        diag(P0) <- w/beta
    }
    Pi.vec <- matrix(Pi,ncol=1)
    P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(rho*0.2+(1-rho))^2)
    if((rho >0) && (!simple) && (!power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(0.6)^2)   }
    if((rho >0) && (!simple) && (power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*((1.285)^2))   }

    M <- matrix(0,n,K)
    membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)
    for(i in 1:n){
        M[i,membership[i]] <- 1
    }
    A.bar <- M%*%P%*%t(M)
    node.degree <- rep(1,n)
    if(rho>0){
    if(simple){
        node.degree[runif(n)<rho] <- 0.2
    }else{
        if(power==FALSE){
            node.degree <- runif(n)*0.8 + 0.2
        }else{
            MM <- ceiling(n/300)
            degree.seed <- rplcon(300,1,alpha)
            node.degree <- sample(rep(degree.seed,MM)[1:n],size=n,replace=FALSE)
        }
    }}
    DD <- diag(node.degree)
    A.bar <- DD%*%A.bar%*%DD
    diag(A.bar) <- 0
    #avg.d <- mean(colSums(A.bar))
    #A.bar <- A.bar*lambda/avg.d
    upper.index <- which(upper.tri(A.bar))
    upper.p <- A.bar[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    g <- membership
    mm <- degree.norm(node.degree,g,K)
    return(list(A=A,g=membership,P=A.bar,degree.norm=mm))
}


degree.norm <- function(d,g,K){

phi <- rep(0,K)
for(k in 1:K){
    index <- which(g==k)
    tmp <- d[index]
    tmp <- tmp/sum(tmp)
    phi[k] <- min(tmp)
}
return(min(phi))
}
ArashSBM.full <- function(lambda,n,beta=0,K=3,w=rep(1,K),Pi=rep(1,K)/K){
    P0 <- diag(w)
    if(beta > 0){
        P0 <- matrix(1,K,K)
        diag(P0) <- w/beta
    }
    Pi.vec <- matrix(Pi,ncol=1)
    P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec))
    M <- matrix(0,n,K)
    membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)
    for(i in 1:n){
        M[i,membership[i]] <- 1
    }
    A.bar <- M%*%P%*%t(M)
    A.bar.2 <- A.bar
    diag(A.bar) <- 0
    upper.index <- which(upper.tri(A.bar))
    upper.p <- A.bar[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    return(list(A=A,P=P,membership=membership,A.bar=A.bar.2))
}



BlurredSBM <- function(lambda,n,beta=0,K=3,w=rep(1,K),Pi=rep(1,K)/K,rho=0.2){
    P0 <- diag(w)
    if(beta > 0){
        P0 <- matrix(1,K,K)
        diag(P0) <- w/beta
    }
    theta <- diag(runif(K)*(1-rho)+rho)
    for(i in 1:(K-1)){
        theta[i,i+1] <- runif(1)*rho
        #theta[i,i+1] <- runif(1)*rho/2

    }
    P0 <- theta%*%P0%*%t(theta)
    Pi.vec <- matrix(Pi,ncol=1)
    P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec))
    M <- matrix(0,n,K)
    membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)
    for(i in 1:n){
        M[i,membership[i]] <- 1
    }
    A.bar <- M%*%P%*%t(M)
    A.bar.2 <- A.bar
    diag(A.bar) <- 0
    upper.index <- which(upper.tri(A.bar))
    upper.p <- A.bar[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    return(list(A=A,P=P,A.bar=A.bar.2))
}


GenericLowRank <- function(n,K){
    Z <- matrix(abs(runif(n*K)),nrow=n,ncol=K)
    S <- Z%*%t(Z)
    P <- S/max(S)

    A <- matrix(0,n,n)
    upper.tri.index <- which(upper.tri(P))
    prob.seq <- P[upper.tri.index]
    T <- length(prob.seq)
    binary.seq <- rep(0,T)
    for(t in 1:T){
        binary.seq[t] <- rbinom(1,1,prob=prob.seq[t])
    }
    A[upper.tri.index] <- binary.seq
    A <- A + t(A)
    return(A)
}


GenericDirectedLowRank <- function(n,K){
    Z1 <- matrix(abs(runif(n*K)),nrow=n,ncol=K)
    Z2 <- matrix(abs(runif(n*K)),nrow=n,ncol=K)
    S <- Z1%*%t(Z2)
    P <- S/max(S)

    A <- matrix(0,n,n)
    R <- matrix(runif(n^2),n,n)
    A[R<P] <- 1
    return(A)
}


GenericDirectedLowRank.full <- function(n,K){
    Z1 <- matrix(abs(runif(n*K)),nrow=n,ncol=K)
    Z2 <- matrix(abs(runif(n*K)),nrow=n,ncol=K)
    S <- Z1%*%t(Z2)
    P <- S/max(S)

    A <- matrix(0,n,n)
    R <- matrix(runif(n^2),n,n)
    A[R<P] <- 1
    return(list(A=A,P=P))
}





ArashSBMX <- function(lambda,n,beta=0,K=3,w=rep(1,K),Pi=rep(1,K)/K,rho=0,simple=TRUE,power=TRUE,p=1,sigma=1){
    P0 <- diag(w)
    if(beta > 0){
        P0 <- matrix(1,K,K)
        diag(P0) <- w/beta
    }
    Pi.vec <- matrix(Pi,ncol=1)
    P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(rho*0.2+(1-rho))^2)
    if((rho >0) && (!simple) && (!power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*(0.6)^2)   }
    if((rho >0) && (!simple) && (power)){ P <- lambda*P0/((n-1)*as.numeric(t(Pi.vec)%*%P0%*%Pi.vec)*((0.2364)^2))   }

    M <- matrix(0,n,K)
    membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)
    for(i in 1:n){
        M[i,membership[i]] <- 1
    }
    A.bar <- M%*%P%*%t(M)
    node.degree <- rep(1,n)
    if(rho>0){
    if(simple){
        node.degree[runif(n)<rho] <- 0.2
    }else{
        if(power==FALSE){
            node.degree <- runif(n)*0.8 + 0.2
        }else{
            node.degree <- rexp(n,rate=3*log(10))
            node.degree[node.degree<0.2] <- 0.2
            node.degree[node.degree>1] <- 1
        }
    }}
    DD <- diag(node.degree)
    A.bar <- DD%*%A.bar%*%DD
    diag(A.bar) <- 0
    upper.index <- which(upper.tri(A.bar))
    upper.p <- A.bar[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)


    X <- matrix(rnorm(n*p),nrow=n,ncol=p)*sigma
    MU <- seq(-1,1,length.out=K)
    for(i in 1:n){
        if(runif(1)>=0.2){
        X[i,] <- X[i,] + MU[membership[i]]
        }
    }

    return(list(A=A,X=X,g=membership))
}
