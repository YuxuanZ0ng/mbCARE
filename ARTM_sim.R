# R code for generating time series data from an ARTM model with two states
library(MASS)
library(dplyr)


ARTM_sim <- function(n, K, L, beta_1k = 1, pi_DA, pi_trans, FC = 2, confounder = "FALSE", covariate = c("TRUE", "FALSE")){
  
  # The number of visits is the same
  t <- rep(1:L, n)
  L <- rep(L, n)
  N <- length(t)
  ID <- 1:N 
  
  samp_info <- data.frame(
    ID = ID,
    sub = rep(1:n, times = L),
    t = t
  )
  
  min_time_indices <- samp_info %>%
    group_by(sub) %>%
    slice(which.min(t)) %>%
    ungroup() %>%
    pull(ID) %>%
    match(samp_info$ID)
  
  # Generate variables
  X_it <- matrix(rbinom(N, size=1, prob=0.5), N, 1)
  beta_1k <- matrix(rep(beta_1k, K), ncol = 1)
  phi_k <- 1  
  Y_k0 <- rnorm(K, 3, 1)
  
  if(confounder == "TRUE"){
    V_it <- cbind(rep(1, N), rnorm(N))
    beta_k <- matrix(c(runif(K, 3, 10), rnorm(K, 0, 1)), ncol = 2)
  }else{
    V_it <- matrix(rep(1, N), ncol = 1)
    beta_k <- matrix(runif(K, 3, 10), ncol = 1)
  }
  
  
  if(covariate == "TRUE"){
    U_it <- cbind(rep(1, N), matrix(rbinom(N, size=1, prob=0.5), N, 1))
    # Generate parameters related to identity distribution 
    alpha0_k <- matrix(c(rep(0, K), rep(0, K)), nrow = K, ncol = 2)  
    alpha01_k <- matrix(c(rep(0, K), rep(0, K - 0.4 * K), runif(0.3 * K, 1, 2), rep(0, 0.1 * K)), nrow = K, ncol = 2)  
    alpha10_k <- matrix(c(rep(0, K), rep(0, K - 0.3 * K), runif(0.3 * K, 1, 2)), nrow = K, ncol = 2)
  }else{
    U_it <- matrix(rep(1, N), ncol = 1)
    # Generate parameters related to identity distribution 
    alpha0_k <- matrix(rnorm(K), nrow = K, ncol = 1)  
    alpha01_k <- matrix(rnorm(K), nrow = K, ncol = 1)
    alpha10_k <- matrix(rnorm(K), nrow = K, ncol = 1)
  }
  
  
  # generate H
  if(pi_DA == 0){
    H <- rep(0, K)
    Delta_k <- rep(0, K)
  }else{
    H <- rep(0, K)
    Delta_k <- rep(0, K)
    C1 <- 1 : floor(pi_DA * K)
    H[C1] <- 1
    Delta_k[C1] <- FC
  }
  
  
  d_it <- matrix(rnorm(N, X_it * FC / 2, 1), ncol = 1)
  
  # Initialize the Y, Z matrices
  Y <- matrix(0, nrow = N, ncol = K)
  Z <- matrix(0, nrow = N, ncol = K)
  
  
  tau <- samp_info$t - c(NA, samp_info$t[-N])
  tau[min_time_indices] <- NA
  
  
  # Simulate the data
  for (i in 1:n) {
    for (k in 1:K) {
      s <- samp_info$ID[samp_info$sub==i][1]
      # Simulate the initial hidden state Z using the first set of alphas
      Z[s, k] <- rbinom(1, 1, 1-plogis(t(U_it[s, ]) %*% alpha0_k[k, ]))
      mu <-  Y_k0[k] * exp(-1 * beta_1k[k]) + X_it[s, ] * Delta_k[k] + t(V_it[s, ]) %*% beta_k[k, ] + d_it[s, ]
      Y[s, k] <- ifelse(Z[s, k] == 1, rnorm(1, mu, phi_k), 0)
      
      for (l in 2:L[i]) {
        s <- samp_info$ID[samp_info$sub==i][l]
        # Transition probabilities for Z
        p01 <- plogis(t(U_it[s, ]) %*% alpha01_k[k, ])
        p10 <- plogis(t(U_it[s, ]) %*% alpha10_k[k, ])
        # Update Z based on transition probabilities
        Z[s, k] <- ifelse(Z[s-1, k] == 0, rbinom(1, 1, p01), rbinom(1, 1, 1-p10))
        
        # Calculate mu_ikt for the observation equation
        mu <- Y[s-1, k] * exp(-beta_1k[k] * tau[s]) + X_it[s, ] * Delta_k[k] + t(V_it[s, ]) %*% beta_k[k, ] + d_it[s, ]
        # Simulate Y based on Z
        Y[s, k] <- ifelse(Z[s, k] == 1, rnorm(1, mu, phi_k), 0)
      }
    }
  }
  
  
  count <- ceiling(exp(Y))
  count[Z == 0] <- 0
  count[count == Inf] <- max(count[count != Inf])
  eta = sum(Z == 0) / (N * K)
  
  dimnames(U_it) <- list(
    sample = paste("samp", 1:N),
    dim = paste("p", 1:dim(U_it)[2])        
  )
  
  dimnames(V_it) <- list(
    sample = paste("samp", 1:N),
    dim = paste("p", 1:dim(V_it)[2])        
  )
  
  dimnames(beta_k) <- list(
    taxon = paste("taxon", 1:K),  
    dim = paste("p", 1:dim(V_it)[2])        
  )
  
  dimnames(Y) <- list(
    sample = paste("samp", 1:N),
    taxon = paste("taxon", 1:K)      
  )
  
  dimnames(count) <- list(
    sample = paste("samp", 1:N),
    taxon = paste("taxon", 1:K)        
  )
  
  dimnames(Z) <- list(
    sample = paste("samp", 1:N),
    taxon = paste("taxon", 1:K)        
  )
 
  dimnames(alpha0_k) <- list(
    taxon = paste("taxon", 1:K),  
    dim = paste("p", 1:dim(U_it)[2])        
  )
  
  dimnames(alpha01_k) <- list(
    taxon = paste("taxon", 1:K),  
    dim = paste("p", 1:dim(U_it)[2])        
  )
  
  dimnames(alpha10_k) <- list(
    taxon = paste("taxon", 1:K),  
    dim = paste("p", 1:dim(U_it)[2])        
  )

  names(H) <- paste0("taxon", seq(K))
  
  names(Delta_k) <- paste0("taxon", seq(K))
  
  names(beta_1k) <- paste0("taxon", seq(K))
  
  names(d_it) <- paste("samp", 1:N)
  
  names(X_it) <- paste("samp", 1:N)
  
  eta = sum(Z == 0) / (N * K)
  test.data=list(samp_info, count, Y, Z, U_it, V_it, X_it, d_it, alpha0_k, alpha01_k, alpha10_k, beta_1k, beta_k, Delta_k, H, eta)
  names(test.data)=c("samp_info", "count", "Y", "Z", "U", "V", "X", "d", "alpha0", "alpha01", "alpha10", "beta1", "beta", "Delta", "H", "eta")
  
  return(test.data)
}

