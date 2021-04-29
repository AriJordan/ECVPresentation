library(pROC)
library(softImpute)
source("SP.R")
source("evaluate_estimation.R")


source("SimonFunk.R")




holdout.evaluation.reg.DC <- function(holdout.index,A,h.seq,K,soft=TRUE,fast=FALSE,p.sample=1){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    A.new <- matrix(0,n,n)
    A.new[upper.tri(A.new)] <- A[edge.index]
    A.new[edge.index[holdout.index]] <- NA
    A.new <- A.new + t(A.new)
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    non.miss <- which(!is.na(A.new))
    #A.new[non.miss] <- A.new[non.miss] + 0.5
    partialacc <- matchness <- max.dis <- avg.dis <- sl2 <- l2 <- adj.l2 <- adj.loglike <- loglike <- sloglike <- rep(0,length(h.seq))

    if(K==1){
	if(soft){
	    SVD <- softImpute(x=as.matrix(A.new),rank.max=K,maxit=300,type="svd")
	    tmp.est <- list(A=matrix(SVD$u,ncol=K)%*%t(matrix(SVD$v,ncol=K))*as.numeric(SVD$d),SVD=SVD)

	}else{
		SVD <- softImpute(x=as.matrix(A.new),rank.max=K,maxit=500,type="svd")
	    tmp.est <- list(A=matrix(SVD$u,ncol=K)%*%t(matrix(SVD$v,ncol=K))*as.numeric(SVD$d),SVD=SVD)
	     tmp.est <- iter.SVD.core(A.new,K=K,init=tmp.est$SVD,fast=fast,p.sample=p.sample)
	}
    }else{
	if(soft){
	    SVD <- softImpute(x=as.matrix(A.new),rank.max=K,maxit=500,type="svd")
	    tmp.est <- list(A=SVD$u%*%diag(SVD$d)%*%t(SVD$v),SVD=SVD)

	}else{
	    tmp.est <- iter.SVD.core(A.new,K=K,fast=fast,p.sample=p.sample)#SF.SVD(A.new,K=k,verbose=FALSE,nstart=5)
	}
    }
    A.approx <- tmp.est$A
    #A.approx[A.approx<0] <- 0
    #A.approx[A.approx>1] <- 1
    A.approx[A.approx<(h*avg.d/n)] <- h*avg.d/n
    A.approx[A.approx>(1+h*avg.d/n)] <- 1+h*avg.d/n

    TMP.A <- A.approx
    TMP.A[non.miss] <- A[non.miss]
    #TMP.A[TMP.A <= 0.5] <- 0
    #TMP.A[TMP.A > 0.5] <- 1
    approx.degree <- rowSums(A)
    min.d <- 0
    D.tau <- mean(approx.degree)
    if(min(approx.degree)<0) min.d <- min(approx.degree)
    for(k in 1:length(h.seq)){
                print(k)
        h <- h.seq[k]

    approx.degree <- rowSums(TMP.A)+h*D.tau
    TMP.L <- diag(1/sqrt(approx.degree))%*%(TMP.A)%*%diag(1/sqrt(approx.degree))
    #TMP.L <- TMP.A + h*D.tau/n
        ## ACASC
        tilde.L <- TMP.L
        #svd.L <- svd((tilde.L+t(tilde.L))/2,nu=K,nv=K)
        svd.L <- irlba(tilde.L,nu=K,nv=K)
        L.U.approx <- matrix(svd.L$v[,1:K],ncol=K)

        V <- L.U.approx
        #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        if(K==1) {V.norms <- as.numeric(abs(V))}else{
            V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        }

        iso.index <- which(V.norms==0)
        Psi <- V.norms
        Psi <- Psi / max(V.norms)
        inv.V.norms <- 1/V.norms
        inv.V.norms[iso.index] <- 1

        V.normalized <- diag(as.numeric(inv.V.norms))%*%V

        #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        #Psi <- V.norms
        #Psi <- Psi / max(V.norms)
        #V.normalized <- diag(1/V.norms)%*%V
        Psi.outer <- outer(Psi,Psi)




        km <- kmeans(V.normalized,centers=K,nstart=30,iter.max=200)#,algorithm="Lloyd")

        group.dis <- rep(0,K)
        Ind.holdout <- A.holdout <- matrix(0,n,n)
        Ind.holdout[Omega] <- 1
        A.holdout[Omega] <- A[Omega]
        cluster.list <- list()
        for(i in 1:K){
              cluster.list[[i]] <- which(km$cluster==i)
        }
        count.mat <- sum.mat <- matrix(0,K,K)
        for(i in 1:K){
              sum.connect.same <- sum(A.holdout[cluster.list[[i]],cluster.list[[i]]])
              count.same <- sum(Ind.holdout[cluster.list[[i]],cluster.list[[i]]])
              sum.mat[i,i] <- sum.connect.same
              count.mat[i,i] <- count.same
              for(j in (i+1):K){
                  if(j <= K){
                       sum.mat[i,j] <- sum.ij <- sum(A.holdout[cluster.list[[i]],cluster.list[[j]]])
                       count.mat[i,j] <- count.ij <- sum(Ind.holdout[cluster.list[[i]],cluster.list[[j]]])

                  }
              }
              sum.mat <- sum.mat + t(sum.mat)
              diag(sum.mat) <- diag(sum.mat)/2
              count.mat <- count.mat + t(count.mat)
              diag(count.mat) <- diag(count.mat)/2
              count.mat <- count.mat + 0.01

        }
      sum.mat <- sum.mat + t(sum.mat)
      diag(sum.mat) <- diag(sum.mat)/2
      count.mat <- count.mat + t(count.mat)
      diag(count.mat) <- diag(count.mat)/2
      count.mat <- count.mat + 0.001
      for(i in 1:K){
          group.dis[i] <- sum.mat[i,i]/count.mat[i,i] - max(sum.mat[i,-i]/count.mat[i,-i])
      }
      all.reg.A <- A
      all.degree <- colSums(all.reg.A+h*D.tau)

      all.L <- diag(sqrt(1/all.degree))%*%(all.reg.A)%*%diag(sqrt(1/all.degree))
      all.L.v <- irlba(all.L,nu=K,nv=K)$v[,1:K]
      all.L.v.norm <- apply(all.L.v,1,function(x) sqrt(sum(x^2)))
      all.L.v.normalized <- diag(all.L.v.norm) %*% all.L.v

      full.km <- kmeans(all.L.v.normalized,centers=K,nstart=20,iter.max=200)
      partial.diff <- PairMatch(g1=full.km$cluster,g2=km$cluster)
      partial.Z <- full.Z <- matrix(0,n,K)
      partial.Z[cbind(1:n,km$cluster)] <- 1
      full.Z[cbind(1:n,full.km$cluster)] <- 1
      partial.acc <- accuracy(partial.Z,full.Z)

        B <- matrix(0,K,K)
        Psi.outer <- outer(Psi,Psi)
        Theta <- matrix(0,n,K)
        for(i in 1:K){
            for(j in i:K){
                N.i <- which(km$cluster==i)
                N.j <- which(km$cluster==j)
                if(i!=j){
                    B[i,j] <- B[j,i] <- (sum(A.new[N.i,N.j],na.rm=TRUE)+1)/(sum(Psi.outer[which(!is.na(A.new[N.i,N.j]))])+1)
                } else{
                   B[i,j] <- B[j,i] <- (1+sum(A.new[N.i,N.j],na.rm=TRUE))/(sum(Psi.outer[which(!is.na(A.new[N.i,N.j]))])-sum(diag(Psi.outer)[which(!is.na(diag(A.new[N.i,N.j])))])+1)#(sum(!is.na(A.new[N.i,N.j])) -sum(!is.na(diag(A.new[N.i,N.j]))))
                }

            }
            Theta[N.i,i] <- 1
        }
        P.hat <- diag(Psi)%*%Theta%*%B%*%t(Theta)%*%diag(Psi)
        diag(P.hat) <- 0
        #dc.block.sq.err[k] <- sum((P.hat[Omega]-A[Omega])^2)
        P.hat.Omega <- P.hat[Omega]
        A.Omega <- A[Omega]
        P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
        P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6

        l2[k] <- sqrt(sum((A.Omega - P.hat.Omega)^2))
        loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
        avg.dis[k] <- mean(group.dis)
        max.dis[k] <- max(group.dis)
        matchness[k] <- partial.diff
        partialacc[k] <- partial.acc
    }
    return(list(loglike=loglike,l2=l2,max.dis=max.dis,avg.dis=avg.dis,matchness=matchness,partialacc = partialacc))
}





holdout.evaluation.Arash.reg.DC.Weight <- function(holdout.index,A,h.seq,K,soft=FALSE,fast=FALSE,p.sample=1){
    n <- nrow(A)
    avg.d <- mean(colSums(A))
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    #tmp.est0 <- iter.SVD.core(A.new,K=K)


    gap <- partialacc <- matchness <- max.dis <- avg.dis <- sl2 <- l2 <- adj.l2 <- adj.loglike <- loglike <- sloglike <- rep(0,length(h.seq))
    A.new <- matrix(0,n,n)
    A.new[upper.tri(A.new)] <- A[edge.index]
    A.new[edge.index[holdout.index]] <- NA
    A.new <- A.new + t(A.new)
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    non.miss <- which(!is.na(A.new))

           edge.label.core <- matrix(c(1,2,3,2,4,5,3,5,6),3,3)
            if(K==2)  {edge.label.core <- matrix(c(1,2,2,3),2,2)}
        label.pair <- outer(as.character(1:K),as.character(1:K),function(x,y)return(paste(x,y,sep="+")))
            value.label.pair <- as.numeric(factor(label.pair))
            #mat.label.pair <- matrix(0,K,K)
            mat.label.pair <- matrix(value.label.pair,K,K)
            #mat.label.pair
            edge.label.core <- matrix(0,K,K)
            upper.index <- which(upper.tri(edge.label.core))
            edge.label.core[upper.index] <- mat.label.pair[upper.index]
            edge.label.core <- edge.label.core+t(edge.label.core)
            diag(edge.label.core) <- diag(mat.label.pair)

#tmp.est0 <- iter.SVD.core(A.new,K=K)
for(k in 1:length(h.seq)){
            print(k)
    h <- h.seq[k]

    if(K==1){
	if(soft){
	    SVD <- softImpute(x=as.matrix(A.new),rank.max=K,maxit=300,type="svd")
	    tmp.est <- list(A=matrix(SVD$u,ncol=K)%*%t(matrix(SVD$v,ncol=K))*as.numeric(SVD$d),SVD=SVD)

	}else{

	     tmp.est <- iter.SVD.core(A.new,K=K,fast=fast,p.sample=p.sample)
	}
    }else{
	if(soft){
	    SVD <- softImpute(x=as.matrix(A.new),rank.max=K,maxit=500,type="svd")
	    tmp.est <- list(A=SVD$u%*%diag(SVD$d)%*%t(SVD$v),SVD=SVD)

	}else{
	    tmp.est <- iter.SVD.core(A.new+h*avg.d/n,K=K,tau=h*avg.d/n,fast=fast,p.sample=p.sample)
	    	}
    }
    A.approx <- tmp.est$A
    A.approx[A.approx<(h*avg.d/n)] <- h*avg.d/n
    A.approx[A.approx>(1+h*avg.d/n)] <- 1+h*avg.d/n
    print("finish completion")
    TMP.A <- A.approx
    approx.degree <- rowSums(TMP.A)
    min.d <- 0
    D.tau <- mean(approx.degree)

    TMP.A.tau <- TMP.A
    approx.degree <- rowSums(TMP.A.tau)
            #print(min(approx.degree))
    TMP.L <- t(t(TMP.A.tau/sqrt(approx.degree))/sqrt(approx.degree))
        tilde.L <- TMP.L
        svd.L <- irlba(tilde.L,nu=K,nv=K)
        L.U.approx <- matrix(svd.L$v[,1:K],ncol=K)

        V <- L.U.approx
        if(K==1) {V.norms <- as.numeric(abs(V))}else{
            V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
        }

        iso.index <- which(V.norms==0)
        Psi <- V.norms
        Psi <- Psi / max(V.norms)
        inv.V.norms <- 1/V.norms
        inv.V.norms[iso.index] <- 1

        V.normalized <- V*as.numeric(inv.V.norms)

        km <- kmeans(V.normalized,centers=K,nstart=30,iter.max=200)
        print("finish training clustering")
      all.reg.A <- A + h*avg.d/n
      all.degree <- colSums(all.reg.A)

      all.L <- t(t(all.reg.A/sqrt(all.degree))/sqrt(all.degree))
      all.L.v <- irlba(all.L,nu=K,nv=K)$v[,1:K]
      all.L.v.norm <- apply(all.L.v,1,function(x) sqrt(sum(x^2)))
      all.L.v.normalized <- all.L.v/all.L.v.norm

      full.km <- kmeans(all.L.v.normalized,centers=K,nstart=100,iter.max=200)
      full.km.Theta <- matrix(0,n,K)
      print("finish full clustering")


        B <- matrix(0,K,K)
        Psi <- rep(0,n)
        #Psi.outer <- outer(Psi,Psi)
        Theta <- matrix(0,n,K)
        Theta[cbind(1:n,km$cluster)] <- 1
        full.km.Theta[cbind(1:n,full.km$cluster)] <- 1
             #edge.label.core
        full.km.edge.mat <- full.km.Theta%*%edge.label.core%*%t(full.km.Theta)
        km.edge.mat <- Theta%*%edge.label.core%*%t(Theta)
        full.km.edge.label <- full.km.edge.mat[edge.index[holdout.index]]
        km.edge.label <- km.edge.mat[edge.index[holdout.index]]
            print("finish assigning labels")
        #partial.acc <- NMI(full.km.edge.label,km.edge.label)
PS <- rep(0,K)
            for(i in 1:K){
                tmp.index <- which(full.km.edge.label==i)
                if(length(tmp.index)>0){
                predicted.label <- km.edge.label[tmp.index]
                pred.membership.mat <- Matrix(0,nrow=length(tmp.index),ncol=max(edge.label.core),sparse=TRUE)
                pred.membership.mat[cbind(1:length(tmp.index),predicted.label)] <- 1
                match.pair.count <- apply(pred.membership.mat,2,function(x){
                                              m <- sum(x)
                                              return(m*(m-1)/2)
                                          })
                ni <- length(tmp.index)
                PS[i] <- sum(match.pair.count)/(ni*(ni-1)/2)
                }else{
                  PS[i] <- 0
                }

            }
        partialacc[k] <- min(PS)
            print("finish calculatign PS!")
      #Z0 <- matrix(0,nrow=length(holdout.index),ncol=max(edge.label.core))
            Z0 <- Matrix(0,nrow=length(holdout.index),ncol=max(edge.label.core),sparse=TRUE)
      Z0[cbind(1:length(holdout.index),full.km.edge.label)] <- 1
      #G0 <- Z0%*%t(Z0)
      #Z1 <- matrix(0,nrow=length(holdout.index),ncol=max(edge.label.core))
            Z1 <- Matrix(0,nrow=length(holdout.index),ncol=max(edge.label.core),sparse=TRUE)
      Z1[cbind(1:length(holdout.index),km.edge.label)] <- 1
      #G1 <- Z1%*%t(Z1)
      #cocluster.err <- sum(abs(G0-G1))/2
            cocluster.err <- norm(t(Z0)%*%Z0,"F")^2+norm(t(Z1)%*%Z1,"F")^2 - 2*norm(t(Z1)%*%Z0,"F")^2
            gap[k] <- cocluster.err

    }
    print("finish calculating gap")
    return(list(partialacc = partialacc,gap=gap))
}






#### Warning: this function is written only for K<=3. For K>3, you will have to change the 'edge.label.core' matrix.

EdgeCV.REG.DC.Weight <- function(A,h.seq,K,cv=NULL,B=3,holdout.p=0.1,soft=FALSE,Arash=TRUE,fast=FALSE){
    n <- nrow(A)
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    holdout.index.list <- list()
    if(is.null(cv)){
        holdout.n <- floor(holdout.p*edge.n)

        for(j in 1:B){
            holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
        }
    }else{
        sample.index <- sample.int(edge.n)
        max.fold.num <- ceiling(edge.n/cv)
        fold.index <- rep(1:cv,each=max.fold.num)[edge.n]
        cv.index <- fold.index[sample.index]
        B <- cv
        for(j in 1:B){
            holdout.index.list[[j]] <- which(cv.index==j)
        }
    }
    if(Arash){
        result <- lapply(holdout.index.list,holdout.evaluation.Arash.reg.DC.Weight,A,h.seq,K,soft,fast=fast,p.sample=1-holdout.p)
    }else{
        result <- lapply(holdout.index.list,holdout.evaluation.reg.DC,A,h.seq,K,soft,fast=fast,p.sample=1-holdout.p)
    }

        gap.mat <- partial.mat  <- matrix(0,nrow=B,ncol=length(h.seq))
    no.edge.seq <- rep(0,B)
    for(b in 1:B){
        partial.mat[b,] <- result[[b]]$partialacc
        gap.mat[b,] <- result[[b]]$gap
    }

    gap.value <- colMeans(gap.mat)
    partial.value <- colMeans(partial.mat)
    gap.min <- min(gap.value)
    gap.which.min <- which.min(gap.value)
    gap.min.index <- apply(gap.mat,1,which.min)
    gap.min.avg <- round(mean(gap.min.index))
    tmp.table <- table(gap.min.index)
    gap.min.stable <- as.numeric(names(tmp.table)[which.max(tmp.table)])
    partial.max <- max(partial.value)
    partial.which.max <- which.max(partial.value)
    partial.max.index <- apply(partial.mat,1,which.max)
    partial.max.avg <- round(mean(partial.max.index))
    tmp.table <- table(partial.max.index)
    partial.max.stable <- as.numeric(names(tmp.table)[which.max(tmp.table)])

    return(list(partial=colMeans(partial.mat),gap.which.min=gap.which.min,partial.which.max=partial.which.max,gap.min.stable=gap.min.stable,gap.min.avg=gap.min.avg,partial.max.stable=partial.max.stable,partial.max.avg=partial.max.avg))
}
