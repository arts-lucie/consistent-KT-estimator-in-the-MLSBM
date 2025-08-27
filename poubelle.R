######################
# Variational M step #
######################
VMstep<-function(X, nX, N, Q, tau, pCts0, directed=FALSE){

  pCts <- list()
  pCts$n <- pCts0$n + colSums(tau)

  pCts$eta <- lapply(X, function(X_t) t(tau) %*% X_t %*% tau)
  if ( ! directed ){
    pCts$eta <- lapply(pCts$eta, function(eta_t) {
      diag(eta_t) <- diag(eta_t) / 2
      return(eta_t)
    })}
  T <- length(X)
  pCts$eta <- lapply(1:T, function(t) {pCts$eta[[t]] + pCts0$eta[[t]]})


  pCts$zeta <- lapply(nX, function(nX_t) t(tau) %*% nX_t %*% tau)
  if( ! directed ){
    pCts$zeta <- lapply(pCts$zeta, function(zeta_t) {
      diag(zeta_t) <- diag(zeta_t) / 2
      return(zeta_t)
    })}

  pCts$zeta <- lapply(1:T, function(t) {pCts$zeta[[t]] + pCts0$zeta[[t]]})

  return(pCts)
}

###############################
# Variational E step#
###############################
Estep<- function(N,Q,X,D,E,CstMat,mincut,maxcut,dmaxInner,fpnbiter,directed,Tau){
  T <- length(X)
  TauD <- replicate(T, matrix(0, nrow=N, ncol=Q), simplify=FALSE)
  TauE <- replicate(T, matrix(0, nrow=N, ncol=Q), simplify=FALSE)
  Tau_tE <- replicate(T, matrix(0, nrow=N, ncol=Q), simplify=FALSE)
  TauOld <- matrix(0, N, Q)

  dGamma <- N * Q
  dGammaOld <- 0.9
  eps <- .Machine$double.eps
  MaxDouble <- .Machine$double.xmax
  old_decrease <- TRUE
  decrease <- TRUE
  go_on <- TRUE

  for(i in 0:(fpnbiter-1)){
    if((dGamma < MaxDouble) && (dGamma > dmaxInner) && go_on){
      dGammaOld <- dGamma
      TauOld <- Tau

      TauD <- lapply(D, function(D_t) Tau %*% D_t)
      TauE <- lapply(E, function(E_t) Tau %*% E_t)

      if (directed) {
        Tau_tE <- lapply(E, function(E_t) Tau %*% t(E_t))
      }

      Tau <- CstMat + sum(sapply(1:T, function(t) { RestrictSum(TauD[[t]]) + RestrictProd(X[[t]], TauE[[t]])}))

      if (directed) {
        Tau <- Tau + sum(sapply(1:T, function(t) { RestrictProd(X[[t]], Tau_tE[[t]])}))
      }

      # Restrictions on Tau
      Tau <- pmax(Tau, mincut)
      Tau <- pmin(Tau, maxcut)
      Tau <- exp(Tau)
      Tau <- NormalizeRows(Tau)
      Tau <- pmax(Tau, eps)

      dGamma <- sum(abs(Tau - TauOld))

      if (i > 1) {
        old_decrease <- decrease
        decrease <- (dGamma - dGammaOld) < 0
        if (old_decrease && !decrease) {
          go_on <- FALSE
        }
      } else if (i == 1) {
        decrease <- (dGamma - dGammaOld) < 0
      }
    }
  }
  return(Tau)
}

###############################
# NormalizeRows #
###############################
NormalizeRows <- function(mat) {
  row_sums <- rowSums(mat)
  norm <- mat / row_sums
  norm[is.nan(norm)] <- 0
  return(norm)
}

###############################
# RestrictProd #
# Compute the product excluding the diagonal Cij = sum( Aik*Bkj, k = 1,.., i-1,i+1,..,col)
###############################
RestrictProd <- function(A, B) {
  diag(A) <- 0
  C <- A %*% B
  return(C)
}

###############################
# RestrictSum #
# A function that receives an n×q matrix and outputs an n×q matrix, where each element C(i,j) is the sum of all elements in column j except for row i.
###############################
RestrictSum <- function(A) {
  sum_col <- colSums(A)
  result <- matrix(rep(sum_col, each = nrow(A)), nrow = nrow(A), byrow = FALSE) - A
  return(result)
}


###############################
# Variational E step final #
###############################
VEstep<-function( X, N, Q, mincut, maxcut, tau, pCts, dmaxInner, fpnbiter,
                  directed=FALSE) {

  dGamma <- 1
  T <- length(X)

  A <- lapply(1:T, function(t) digamma(pCts$eta[[t]]))
  B <- lapply(1:T, function(t) digamma(pCts$zeta[[t]]))
  C <- lapply(1:T, function(t) digamma(pCts$eta[[t]] + pCts$zeta[[t]]))

  if (directed)
    D <- lapply(1:T,function(t) B[[t]] + t(B)[[t]] - C[[t]] - t(C)[[t]])
  else
    D <- lapply(1:T,function(t) B[[t]] - C[[t]] )

  E <- lapply(1:T,function(t) A[[t]] - B[[t]])

  digamma.n     <- digamma(pCts$n)
  digamma.sum.n <- digamma(sum(pCts$n))

  CstMat <- matrix( rep(digamma.n - digamma.sum.n, N), N, Q, byrow=TRUE)

  tau <- Estep(N,Q,X,D,E,CstMat,mincut,maxcut,dmaxInner,fpnbiter,directed,tau)

  return( tau )
}

###############
# Lower Bound #
###############
lowerBound_KT<-function(tau, pCts0, pCts, directed=FALSE){
  T <- length(pCts$eta)
  Q <- dim(tau)[2]

  if ( directed ) {
    KT <-lgamma(sum(pCts0$n)) - lgamma(sum(pCts$n)) - sum(lgamma(pCts0$n) - lgamma(pCts$n)) + sum(sapply(1:T, function(t) {
      sum(lbeta(pCts$eta[[t]], pCts$zeta[[t]]) - lbeta(pCts0$eta[[t]], pCts0$zeta[[t]]))
    }))
  } else {
    A <- lapply(1:T, function(t) {lbeta(pCts$eta[[t]], pCts$zeta[[t]])})
    B <- lapply(1:T, function(t) {lbeta(pCts0$eta[[t]], pCts0$zeta[[t]])})
    KT <-  lgamma(sum(pCts0$n)) - lgamma(sum(pCts$n)) - sum(lgamma(pCts0$n) - lgamma(pCts$n)) + sum(sapply (1:T , function(t) {
      (1/2)*sum(A[[t]][-seq(1, Q^2, by=Q+1)]) + sum(A[[t]][seq(1, Q^2, by=Q+1)]) - (1/2)*sum(B[[t]][-seq(1, Q^2, by=Q+1)]) - sum(B[[t]][seq(1, Q^2, by=Q+1)])
      }))
  }
  return(c(KT- sum(tau*log(tau)),KT))
}

#####################
# maj (majority vote) #
#####################
maj <- function(elem_i){
  unique <- unique(elem_i)
  unique[which.max(tabulate(match(elem_i,unique)))]
}


#####################
# Variational Bayes #
#####################
VariationalBayes <-
  function(N,m, qmin, qmax, nbiter, fpnbiter, emeps, fpeps, directed=FALSE, n0, eta0, zeta0) {

    ## Edge matrix to Adjacency matrix
    T <- length(m)
    X <- lapply(m, function(edges) {
      X_t <- matrix(0, N, N)
      X_t[ cbind(edges[1,], edges[2,]) ] <- 1
      if (!directed) {
        X_t[ cbind(edges[2,], edges[1,]) ] <- 1
      }
      return(X_t)
    })

    X <- lapply(X, function(X_t) {
                X_t <- sign(X_t)
                return(X_t)
              }) # check that the graph contains binary edges

    nX <- lapply(X, function(X_t) {
                nX_t <- matrix(1, N, N) - X_t
                diag(nX_t) <- 0
                return(nX_t)
               })

    vbOptions <- list(dmaxOuter = emeps, dmaxInner = fpeps)

    pCts0 <- list()

    # Result index
    i.res <- 1
    y <-  vector("list", qmax-qmin+1)

    clust <- lapply(1:T, function(t) hclust(dist(X[[t]]),method="ward.D2"))
    for (Q in qmin:qmax ) {

      maxcut <- log(.Machine$double.xmax) - log(Q)
      mincut <- log(.Machine$double.xmin)

      pCts0$n    <- rep(n0, Q)
      pCts0$eta  <- lapply(1:T, function(t) matrix(rep(eta0, Q^2), Q, Q))
      pCts0$zeta <- lapply(1:T, function(t) matrix(rep(zeta0, Q^2), Q, Q))

      # tau initialization
      zinit_liste <- lapply(1:T, function(t) cutree(clust[[t]],Q))
      zinit <- sapply(1:N, function(i) {
        elem_i <- sapply(zinit_liste, function(l) l[i])
        maj(elem_i)
      })
      Taus.init <- matrix(0,Q,N)
      Taus.init[cbind(zinit, 1:N)] <- 1

      # Require a N x Q matrix
      tau <- t(Taus.init)

      # Avoid overflow
      tau[tau<.Machine$double.xmin] <- .Machine$double.xmin

      delta <- 1
      i <- 0
      l <- list()
      KT <- list()


      while(  is.finite(delta) == TRUE
              & delta>vbOptions$dmaxOuter
              & (i < nbiter)
      ) {
        i <- i + 1

        pCts <- VMstep( X, nX, N, Q, tau, pCts0, directed )
        l[[i]] <- lowerBound_KT(tau, pCts0, pCts, directed)[1]
        KT[[i]] <- lowerBound_KT(tau, pCts0, pCts, directed)[2]
        tau <- VEstep( X, N, Q, mincut, maxcut, tau, pCts,
                       vbOptions$dmaxInner, fpnbiter, directed )

        if(i > 1) {
          delta <- abs((l[[i]] - l[[i-1]])/l[[i-1]])
        }
      }

      # Only intialization
      if ( nbiter == 0 ) {
        pCts <- pCts0
        l[[1]] <- 0
        KT[[1]] <- 0
      }

      # Store results
      # -------------
      #
      # ICL item is in fact the Bayesian criterion
      y[[i.res]]$criterion <- max( as.numeric(as.matrix(l)) )
      y[[i.res]]$KT <- max( as.numeric(as.matrix(KT)) )
      y[[i.res]]$alphas    <- pCts$n / sum( pCts$n )
      y[[i.res]]$Pis       <- lapply(1:T, function(t) pCts$eta[[t]] / (pCts$eta[[t]]+pCts$zeta[[t]]))
      y[[i.res]]$a <- pCts$n
      y[[i.res]]$eta <- pCts$eta
      y[[i.res]]$zeta <- pCts$zeta
      y[[i.res]]$Taus      <- t(tau)
      i.res <- i.res+1
    }
    return(y)

  }


###################################################################

#####################
# ML SBM #
#####################
library(igraph)

MLsbm <- function(n, pi, P, T){
  k <- length(pi)
  Z <- sample(1:k, n, replace=TRUE, prob=pi) # latent variables
  block_sizes <- tabulate(Z, nbins = k)
  G_ML <- list()
  for (t in 1:T) {
    G <-sample_sbm(n,P[[t]],block_sizes)
    V(G)$group <- rep(1:k,block_sizes)
    G_ML<-append(G_ML,list(G))
  }
  return(G_ML)
}

#####################
# Edge matrix from a ML graph #
#####################

matrice_aretes <- function(ml_graphe) {
  T <- length(ml_graphe)
  list_mat <- lapply(ml_graphe, function(graphe) {
    edgelist <- as_edgelist(graphe)
    matrice <- t(edgelist)
    return(matrice)
  })
  return(list_mat)
}

#####################
# penalization #
#####################
pen_ML <- function(k, n,T, eps = 10^-8) {
  if (k == 1) { return(0) }
  x <- 1:(k - 1)
  pen <- T*(x * (x + 1+1/T) + 1 + eps) / 2 * log(n)
  return(sum(pen))
}

#####################
# Computation of k-KT #
#####################
k_chap_ML <- function(ml_G, qmax=10, nbiter=10,fpnbiter=10,emeps=10^-8,fpeps=10^-8,n0=0.5, eta0=0.5, zeta0=0.5){
  n <- vcount(ml_G[[1]])
  m <- matrice_aretes(ml_G)
  res <- VariationalBayes(N=n,m = m, qmin =1, qmax=min(qmax,n), nbiter=nbiter, fpnbiter=fpnbiter, emeps=emeps, fpeps=fpeps, directed=FALSE, n0=n0, eta0=eta0, zeta0=zeta0)
  log_KT <- sapply(1:min(qmax,n), function(l) res[[l]]$KT)
  pen <- sapply(1:min(qmax,n), function(l) pen_ML(l,n,T))
  test <- log_KT-pen
  k <- which.max(test)
  return(k)
}

