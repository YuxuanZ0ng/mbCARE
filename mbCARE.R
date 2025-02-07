library(foreach)
library(dplyr)
library(stats)

mbCARE <- function(count, X, U, V, samp_info, interval.prob = 0.5, adjust = "BH", signif = 0.05, prev.cut = 0.05) {
  library(foreach)
  library(dplyr)
  library(stats)
  
  # Sort samp_info by sub and t.
  sorted_samp_info <- samp_info %>%
    arrange(sub, t)
  sorted_indices <- match(sorted_samp_info$ID, samp_info$ID)
  X <- X[sorted_indices, , drop = FALSE]
  U <- U[sorted_indices, , drop = FALSE]
  V <- V[sorted_indices, , drop = FALSE]
  count <- count[sorted_indices, , drop = FALSE]
  samp_info <- sorted_samp_info
  
  min_time_indices <- samp_info %>%
    group_by(sub) %>%
    filter(t == min(t)) %>%
    ungroup() %>%
    pull(ID) %>%
    match(samp_info$ID)
  
  n <- length(unique(samp_info$sub))
  L <- sapply(unique(samp_info$sub), function(i){sum(samp_info$sub == i)})
  N <- nrow(count)
  Y <- log(count + 1)
  Y_AR1 <- Y
  Y_AR1[2:N,] <- Y[1:(N-1), ]
  Y_AR1[min_time_indices, ] <- NA 
  tax_keep <- which(apply(Y_AR1, 2, function(x) {mean(x != 0, na.rm = TRUE) > prev.cut})) #filtering rare taxa
  Y <- Y[,tax_keep]
  Y_AR1 <- Y_AR1[,tax_keep]
  K <- dim(Y)[2]  # Number of taxons
  tau <- samp_info$t - c(NA, samp_info$t[-N])
  tau[min_time_indices] <- NA
  tau <- tau / mean(tau, na.rm = T)
  
  p1 <- dim(U)[2] # Number of covariates in state transition process
  p2 <- dim(V)[2] # Number of covariates in Norm
  colnames(V) <- paste0("V", 1:p2)
  
  # taxon with structure zero
  if(length(unique(X)) == 2){
    tax_struc <- c(which(colSums(Y[X == unique(X)[1],], na.rm = TRUE) == 0), which(colSums(Y[X == unique(X)[2],], na.rm = TRUE)==0))
    structure_zero <- (colSums(Y[X == unique(X)[1],], na.rm = TRUE) == 0) - (colSums(Y[X == unique(X)[2],], na.rm = TRUE)==0)
    structure_zero[structure_zero == 0] <- "non-struc-zero"
    structure_zero[structure_zero == 1] <- unique(X)[1]
    structure_zero[structure_zero == -1] <- unique(X)[2]
  }else{
    tax_struc <- integer(0)
    structure_zero <- FALSE
  }
  
  # Initialization of parameters
  Theta <- list(d = rep(0, N),
                alpha = matrix(runif(K*p1), ncol = p1),
                alpha0 = matrix(runif(K*p1), ncol = p1),
                alpha01 = matrix(runif(K*p1), ncol = p1),
                alpha10 = matrix(runif(K*p1), ncol = p1),
                alpha_p = matrix(rep(1, K*p1), ncol = p1),
                alpha0_p = matrix(rep(1, K*p1), ncol = p1),
                alpha01_p = matrix(rep(1, K*p1), ncol = p1),
                alpha10_p = matrix(rep(1, K*p1), ncol = p1),
                alpha_q = matrix(rep(1, K*p1), ncol = p1),
                alpha0_q = matrix(rep(1, K*p1), ncol = p1),
                alpha01_q = matrix(rep(1, K*p1), ncol = p1),
                alpha10_q = matrix(rep(1, K*p1), ncol = p1),
                beta1 = matrix(runif(K), ncol = 1),
                beta = matrix(runif(K*p2), ncol = p2),
                Delta = matrix(rep(0, K), ncol = 1),
                Delta_sd = matrix(rep(1, K), ncol = 1),
                na_keep = matrix(rep(0, K), ncol = 1))
  
  matrices <- c("alpha0", "alpha", "alpha01", "alpha10", "alpha0_p", "alpha_p", "alpha01_p", "alpha10_p", "alpha0_q", "alpha_q", "alpha01_q", "alpha10_q")
  Theta[matrices] <- lapply(Theta[matrices], function(x) {
    dimnames(x)[1] <- dimnames(Y)[2]
    dimnames(x)[2] <- dimnames(U)[2]
    return(x)
  })
  matrices1 <- c("beta1", "beta", "Delta", "Delta_sd", "na_keep")
  Theta[matrices1] <- lapply(Theta[matrices1], function(x) {
    dimnames(x)[1] <- dimnames(Y)[2]
    return(x)
  })
  
  out_gamma_xi <- compute_gamma_xi(Y, Y_AR1, N, K)
  gamma_ikt <- out_gamma_xi$gamma_ikt# Compute gamma
  xi_ikt <- out_gamma_xi$xi_ikt # Compute xi
  
  
  Theta <- update_parameters(gamma_ikt, xi_ikt, Y, Y_AR1, X, U, V, p1, p2, n, K, tau, min_time_indices, Theta)
  
  
  #MCI approach for bias correction
  BC <- list(bias = NA,
             Delta = matrix(rep(NA, K), ncol = 1),
             Delta_p = matrix(rep(NA, K), ncol = 1),
             Delta_q = matrix(rep(NA, K), ncol = 1),
             reject = matrix(rep(NA, K), ncol = 1),
             structure_zero = structure_zero)
  
  matrices2 <- c("Delta_q", "Delta_p", "Delta", "reject")
  BC[matrices2] <- lapply(BC[matrices2], function(x) {
    dimnames(x)[1] <- dimnames(Y)[2]
    return(x)
  })
  
  p <- sort(Theta$Delta[Theta$na_keep == 0])
  K1 <- length(p)
  pp <- function(p,i){p[i + floor(K1 * interval.prob)] - p[i]}
  a <- which.min(sapply(1:(K1 - floor(K1 * interval.prob)), pp, p=p))[1]
  BC$bias <- mean(p[a:(a + floor(K1 * interval.prob))])
  BC$Delta <- Theta$Delta - BC$bias
  BC$Delta_p[, 1] <- sapply(BC$Delta/Theta$Delta_sd, function(x) 2*pnorm(abs(x), mean = 0, sd = 1, lower.tail = F))
  BC$Delta_q[, 1] <- p.adjust(BC$Delta_p, method = adjust)
  BC$Delta[Theta$na_keep == 1, 1] <- NA # failed converge in non-linear regression
  BC$Delta_p[Theta$na_keep == 1, 1] <- 1
  BC$Delta_q[Theta$na_keep == 1, 1] <- 1 
  BC$Delta_p[tax_struc, 1] <- 1  # taxon with structure zero
  BC$Delta_q[tax_struc, 1] <- 1  # taxon with structure zero
  BC$reject[, 1] <- c(BC$Delta_q < signif)
  
  Theta$alpha_q <- apply(Theta$alpha_p, 2, function(x) {p.adjust(x, method = adjust)})
  Theta$alpha0_q <- apply(Theta$alpha0_p, 2, function(x) {p.adjust(x, method = adjust)})
  Theta$alpha01_q <- apply(Theta$alpha01_p, 2, function(x) {p.adjust(x, method = adjust)})
  Theta$alpha10_q <- apply(Theta$alpha10_p, 2, function(x) {p.adjust(x, method = adjust)})
  Theta$beta1[Theta$na_keep == 1, 1] <- NA # failed converge in non-linear regression
  Theta$beta[Theta$na_keep == 1,] <- NA # failed converge in non-linear regression
  Theta$na_keep[tax_struc] <- "Structure Zero"  # taxon with structure zero
  
  # Return the final parameters and estimated hidden states
  out <- list(Theta = Theta,
              bias_correction = BC,
              tax_keep = tax_keep)
  return(out)
  
}


# Compute gamma and xi
compute_gamma_xi <- function(Y, Y_AR1, N, K) {
  gamma_ikt <- array(0, c(N, K))
  xi_ikt <- array(0, c(N, K, 2, 2))
  
  # Compute gamma
  gamma_ikt[,] = ifelse(Y[,] == 0, 1, 0)
  
  # Compute xi for all t except the last one
  xi_ikt[, , 1, 1] <- ifelse(is.na(Y_AR1), NA, ifelse((Y_AR1 == 0 & Y == 0), 1, 0))
  xi_ikt[, , 1, 2] <- ifelse(is.na(Y_AR1), NA, ifelse((Y_AR1 == 0 & Y > 0), 1, 0))
  xi_ikt[, , 2, 1] <- ifelse(is.na(Y_AR1), NA, ifelse((Y_AR1 > 0 & Y == 0), 1, 0))
  xi_ikt[, , 2, 2] <- ifelse(is.na(Y_AR1), NA, ifelse((Y_AR1 > 0 & Y > 0), 1, 0))
  
  out=list(gamma = gamma_ikt, xi = xi_ikt)
  names(out)=c("gamma_ikt", "xi_ikt")
  return(out)
}


# Update Parameters
update_parameters <- function(gamma_ikt, xi_ikt, Y, Y_AR1, X, U, V, p1, p2, n, K, tau, min_time_indices, Theta) {
  
  updated_Theta <- lapply(1:K, function(k) {

    # Update alpha0, alpha01, and alpha10 for taxon k
    logis <- data.frame(gamma_ikt = gamma_ikt[,k], U = U)
    fit_alpha <- glm(gamma_ikt ~ . -1, family=stats::quasibinomial(link = "logit"), control=list(maxit=100), data=logis)
    
    # Update alpha0, alpha01, and alpha10 for taxon k
    logis0 <- data.frame(gamma_ikt = gamma_ikt[min_time_indices, k], U = U[min_time_indices, ])
    fit_alpha0 <- glm(gamma_ikt ~ . -1, family=stats::quasibinomial(link = "logit"), control=list(maxit=100), data=logis0)
    
    # Similar for alpha01 and alpha10
    U_flat <- U[-min_time_indices, ]
    weight <- xi_ikt[-min_time_indices, k, 1, 1]+xi_ikt[-min_time_indices, k, 1, 2]
    if(sum(weight) == 0){
      alpha01 = rep(NA, p1)
      alpha01_p = rep(1, p1)
    }else{
      logis01 <- data.frame(xi_ikt_01 = as.vector(xi_ikt[-min_time_indices,k,1,2]/weight), weight = weight, U_flat = U_flat)
      fit_alpha01 <- glm(xi_ikt_01 ~ U_flat -1, family=stats::quasibinomial(link = "logit"), weights=weight, control=list(maxit=100), data=logis01)
      alpha01 = coef(fit_alpha01)
      alpha01_p = summary(fit_alpha01)[["coefficients"]][,4]
    }
    
    weight <- xi_ikt[-min_time_indices, k, 2, 1]+xi_ikt[-min_time_indices, k, 2, 2]
    if(sum(weight) == 0){
      alpha10 = rep(NA, p1)
      alpha10_p = rep(1, p1)
    }else{
    logis10 <- data.frame(xi_ikt_10 = as.vector(xi_ikt[-min_time_indices,k,2,1]/weight), weight = weight, U_flat = U_flat)
    fit_alpha10 <- glm(xi_ikt_10 ~ U_flat -1, family=stats::quasibinomial(link = "logit"), weights=weight, control=list(maxit=100), data=logis10)
    alpha10 = coef(fit_alpha10)
    alpha10_p = summary(fit_alpha10)[["coefficients"]][,4]
    }
    
    # Update Delta, beta1, and beta for taxon k
    Norm <- data.frame(gamma_ikt = as.vector(1-gamma_ikt[,k]), Y = Y[,k], 
                       V , Y_AR = Y_AR1[,k], tau = tau, X = X)
    Norm_1 <- na.omit(Norm) %>%
      filter(gamma_ikt == 1)
    
    formula_string <- "Y ~"
    for (i in 1:p2) {
      formula_string <- paste0(formula_string, " a", i, " * V", i, " +")
    }
    formula_string <- paste0(formula_string, " exp(- beta_1k * tau) * Y_AR + delta * X")
    
    start_values <- list(beta_1k = 1, delta = 1)
    for (i in 1:p2) {
      start_values[[paste0("a", i)]] <- 1
    }
    model_formula <- as.formula(formula_string)
    
    suppressWarnings(glm <- try(nls(model_formula, 
                                    data = Norm_1,
                                    start = start_values,
                                    na.action = na.exclude), silent = TRUE))
    if (inherits(glm, "try-error")) {
      suppressWarnings(glm <- try(nls(model_formula, 
                                      data = Norm_1,
                                      start = start_values,
                                      algorithm = "plinear",
                                      na.action = na.exclude), silent = TRUE))
      if (inherits(glm, "try-error")) {
        Theta$Delta[k] <- 0
        Theta$Delta_sd[k] <- 0
        Theta$na_keep[k] <- 1
        Theta$beta[k,] <- rep(0, p2)
        Theta$beta1[k,] <- 0
      }else{ 
        Theta$Delta[k] <- ifelse(is.na(summary(glm)[["coefficients"]]["delta",1]), 0, summary(glm)[["coefficients"]]["delta",1])
        Theta$Delta_sd[k] <- ifelse(is.na(summary(glm)[["coefficients"]]["delta",2]), 0, summary(glm)[["coefficients"]]["delta",2])
        Theta$na_keep[k] <- is.na(summary(glm)[["coefficients"]]["delta",1])
        Theta$beta[k,] <- summary(glm)[["coefficients"]][3:(2+p2),1]
        Theta$beta1[k,] <- ifelse(is.na(summary(glm)[["coefficients"]]["beta_1k",1]), 0, summary(glm)[["coefficients"]]["beta_1k",1])
      }
    }else{ 
      Theta$Delta[k] <- ifelse(is.na(summary(glm)[["coefficients"]]["delta",1]), 0, summary(glm)[["coefficients"]]["delta",1])
      Theta$Delta_sd[k] <- ifelse(is.na(summary(glm)[["coefficients"]]["delta",2]), 0, summary(glm)[["coefficients"]]["delta",2])
      Theta$na_keep[k] <- is.na(summary(glm)[["coefficients"]]["delta",1])
      Theta$beta[k,] <- summary(glm)[["coefficients"]][3:(2+p2),1]
      Theta$beta1[k,] <- ifelse(is.na(summary(glm)[["coefficients"]]["beta_1k",1]), 0, summary(glm)[["coefficients"]]["beta_1k",1])
    }
    
    
    return(list(alpha = coef(fit_alpha),
                alpha0 = coef(fit_alpha0),
                alpha01 = alpha01,
                alpha10 = alpha10,
                alpha0_p = summary(fit_alpha0)[["coefficients"]][,4],
                alpha_p = summary(fit_alpha)[["coefficients"]][,4],
                alpha01_p = alpha01_p,
                alpha10_p = alpha10_p,
                beta = Theta$beta[k,],
                beta1 = Theta$beta1[k,],
                Delta = Theta$Delta[k],
                Delta_sd = Theta$Delta_sd[k],
                na_keep = Theta$na_keep[k]))
  })
  
  for (k in 1:K) {
    Theta$alpha[k, ] <- updated_Theta[[k]]$alpha
    Theta$alpha0[k, ] <- updated_Theta[[k]]$alpha0
    Theta$alpha01[k, ] <- updated_Theta[[k]]$alpha01
    Theta$alpha10[k, ] <- updated_Theta[[k]]$alpha10
    Theta$alpha_p[k, ] <- updated_Theta[[k]]$alpha_p
    Theta$alpha0_p[k, ] <- updated_Theta[[k]]$alpha0_p
    Theta$alpha01_p[k, ] <- updated_Theta[[k]]$alpha01_p
    Theta$alpha10_p[k, ] <- updated_Theta[[k]]$alpha10_p
    Theta$beta[k, ] <- updated_Theta[[k]]$beta
    Theta$beta1[k, ] <- updated_Theta[[k]]$beta1
    Theta$Delta[k] <- updated_Theta[[k]]$Delta
    Theta$Delta_sd[k] <- updated_Theta[[k]]$Delta_sd
    Theta$na_keep[k] <- updated_Theta[[k]]$na_keep
  }
  
  
  #Update d_it
  gamma_2 <- 1 - gamma_ikt
  res <- Y - sapply(1:K, function(k) {Y_AR1[,k] * exp(- Theta$beta1[k,] * tau) + X * Theta$Delta[k] + V %*% Theta$beta[k,]})
  Theta$d <- rowSums(gamma_2 * res) / rowSums(gamma_2)
  return(Theta)
}

