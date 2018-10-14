## regmix_em
regmix_em <- function(y, xmat, pi.init, beta.init , sigma.init, 
                     control = list(maxit = 500, tol = 1e-5)) {
  n <- nrow(xmat)
  k <- ncol(beta.init)
  p <- matrix(0, nrow = nrow(xmat), ncol = ncol(beta.init))
  p_nume<-p
  beta1 <- beta.init
  pi1<-pi.init
for (r in 1:control$maxit) {
  for (j in 1:k) {
    p_nume[,j] <- as.matrix((pi.init[j]*(2*3.14159265*sigma.init^2)^(-0.5)*
                 exp(-(y-as.matrix(xmat)%*%as.matrix(beta.init[,j]))^2/2/sigma.init^2)))
  }
  p_deno <- rowSums(p_nume)
  p <- p_nume/p_deno
  pi1 <- colSums(p)/n
  for (j in 1:k) {
    beta1_2nd <- matrix(0, nrow = ncol(xmat), ncol=k)
    for (i in 1:n) {
      beta1_2nd[,j] <- beta1_2nd[,j] + t(xmat[i,]*p[i,j]*y[i])
    }
    beta1[,j] <- (sum(diag(as.matrix(xmat)%*%t(as.matrix(xmat)))
                      *p[,j]))^(-1)*beta1_2nd[,j]
  }
  sigma_2 <- sum( (y-as.matrix(xmat)%*%as.matrix(beta1))^2 * p)
  sigma_1 <- sqrt(sigma_2)
  if ((max(abs(pi1-pi.init)) <= control$tol) &(max(abs(beta1-beta.init)) <= control$tol) &
      (max(sigma_1-sigma.init) <= control$tol )) break
  pi.init<-pi1
  beta.init<-beta1
  sigma.init<-sigma_1
}
  return(list(pi=pi1, beta=beta1, sigma=sigma_1, iteration=r))
}

## regmix_sim
regmix_sim <- function(n, pi, beta, sigma) {
  K <- ncol(beta)
  p <- NROW(beta)
  xmat <- matrix(rnorm(n * p), n, p) # normal covaraites
  error <- matrix(rnorm(n * K, sd = sigma), n, K)
  ymat <- xmat %*% beta + error # n by K matrix
  ind <- t(rmultinom(n, size = 1, prob = pi))
  y <- rowSums(ymat * ind)
  data.frame(y, xmat)
}

## simulation
n <- 400
pi <- c(.3, .4, .3)
bet <- matrix(c( 1,  1,  1, 
                 -1, -1, -1), 2, 3)
sig <- 1
set.seed(1205)
dat <- regmix_sim(n, pi, bet, sig)
regmix_em(y = dat[,1], xmat = dat[,-1], 
           pi.init = pi/pi/length(pi),
           beta.init = bet*1,
           sigma.init = sig / sig, 
           control = list(maxit = 500, tol = 1e-5))

