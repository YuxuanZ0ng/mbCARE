source("mbCARE.R")
source("ANCOMBC2.R")
library(doSNOW)
library(foreach)
library(ggplot2)
library(ggsci) 
library(MASS)
library(dplyr)
library(NBZIMM)
library(ZIBR)


ZIMM_sim <- function(n, K, L, pi_DA, pi_trans, FC = 2, covariate = c("TRUE", "FALSE")){
  
  t <- rep(1:L, n)
  L <- rep(L, n)
  N <- length(t)
  ID <- 1:N 
  
  samp_info <- data.frame(
    ID = ID,
    sub = rep(1:n, times = L),
    t = t
  )
  
  # Generate variables
  X_it <- matrix(rbinom(N, size=1, prob=0.5), N, 1)
  phi_k <- 1  
  V_it <- matrix(rep(1, N), ncol = 1)
  beta_k <- matrix(runif(K, 3, 10), ncol = 1)
  
  
  
  if(covariate == "TRUE"){
    U_it <- cbind(rep(1, N), matrix(rbinom(N, size=1, prob=0.5), N, 1))
    # Generate parameters related to identity distribution 
    alpha_k <- cbind(rep(0, K), c(rep(0, K-pi_trans * K), runif(pi_trans * K, 1, 2)))
  }else{
    U_it <- matrix(rep(1, N), ncol = 1)
    # Generate parameters related to identity distribution 
    alpha_k <- matrix(rnorm(K), nrow = K, ncol = 1)  
  }
  
  # Generate Z 
  Z <- array(0, c(N, K))
  
  # Simulate the data
  for (i in 1:N) {
    for (k in 1:K) {
      Z[i, k] <- rbinom(1, 1, 1-plogis(t(U_it[i, ]) %*% alpha_k[k, ]))
    }
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
  R_i <- rnorm(n, 0, 1)
  R <- rep(R_i, time = L)
  
  # Initialize the Y matrices
  Y <- matrix(0, nrow = N, ncol = K)
  
  
  # Simulate the data
  for (s in 1:N) {
    for (k in 1:K) {
      mu <- X_it[s, ] * Delta_k[k] + t(V_it[s, ]) %*% beta_k[k, ] + R[s] + d_it[s, ]
      Y[s, k] <- ifelse(Z[s, k] == 1, rnorm(1, mu, phi_k), 0)
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
  
  dimnames(alpha_k) <- list(
    taxon = paste("taxon", 1:K),  
    dim = paste("p", 1:dim(U_it)[2])        
  )
  
  names(H) <- paste0("taxon", seq(K))
  
  names(Delta_k) <- paste0("taxon", seq(K))
  
  names(d_it) <- paste("samp", 1:N)
  
  names(X_it) <- paste("samp", 1:N)
  
  eta = sum(Z == 0) / (N * K)
  test.data=list(samp_info, count, Y, Z, U_it, V_it, X_it, d_it, alpha_k, beta_k, Delta_k, H, eta)
  names(test.data)=c("samp_info", "count", "Y", "Z", "U", "V", "X", "d", "alpha", "beta", "Delta", "H", "eta")
  
  return(test.data)
}


# Define the parameters
ncores <- 64
n.sim <- 100
n <- 50 # Number of subjects
K <- 50  # Number of taxons
L <- c(5, 10, 20, 50) # Number of time points
sim_total = n.sim * 4
sim.seed = matrix(1:sim_total, n.sim, 4)

cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
DATA <- list()
DATA <- foreach(i = 1:4, .combine='c') %do% {
  L.i <- L[i]
  foreach(j = sim.seed[,i], .packages = c("MASS", "dplyr"), .combine='c') %dopar% {
    set.seed(j)
    data <- ZIMM_sim(n = n, K = K, L = L.i, pi_DA = 0.2, pi_trans = 0.4, FC = 1, covariate = "TRUE")
    list(data)}
}
stopCluster(cl)


cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
res <- foreach(i = 1:sim_total, .combine = 'cbind', .errorhandling = "pass") %dopar% {
  library(foreach)
  library(dplyr)
  library(doRNG)
  library(MicrobiomeStat)
  library(NBZIMM)
  library(ZIBR)
  
  data=DATA[[i]]
  id <- which(data[["H"]] == 1)
  id0 <- which(data[["alpha"]][,2]!=0)
  
  ############ mbCARE ############
  count <- data[["count"]]
  X <- data[["X"]]
  U <- data[["U"]]
  V <- data[["V"]]
  samp_info <- data[["samp_info"]]
  N <- dim(count)[1]
  K <- dim(count)[2]
  
  
  suppressWarnings(out_mbCARE <- try(mbCARE(count, X, U, V, samp_info, prev.cut = 0), silent = TRUE))
  if (inherits(out_mbCARE, "try-error")) {
    power_mbCARE_q <- 0; fdr_mbCARE_q <- 0
    power_mbCARE_0 <- 0; fdr_mbCARE_0 <- 0 
  }else{
    id_mbCARE_q <- which(out_mbCARE[["bias_correction"]][["Delta_q"]] <= 0.05)
    power_mbCARE_q <- length(intersect(id_mbCARE_q,id))/length(id)
    fdr_mbCARE_q <- 1-length(intersect(id_mbCARE_q,id))/length(id_mbCARE_q)
    id_mbCARE_0 <- which(out_mbCARE[["Theta"]][["alpha_q"]][,2] <= 0.05)
    power_mbCARE_0 <- length(intersect(id_mbCARE_0,id0))/length(id0)
    fdr_mbCARE_0 <- 1-length(intersect(id_mbCARE_0,id0))/length(id_mbCARE_0)
  }
  
  
  ############ ANCOMBC ############
  counts <- t(count)
  colnames(counts) <- paste("sample", 1:N)
  rownames(counts) <- paste("taxon", 1:K)
  ID <- data.frame(
    sample = rownames(count),
    sub = samp_info$sub,  
    time = samp_info$t, 
    X = X,
    U = as.numeric(U[,2])
  )
  assays = S4Vectors::SimpleList(counts = counts)
  smd = S4Vectors::DataFrame(ID)
  tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)
  
  suppressWarnings(out_ANCOMBC <- try(ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                                               fix_formula = "X", 
                                               rand_formula = "(1 | sub)",
                                               p_adj_method = "BH", 
                                               prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                               group = "X", struc_zero = FALSE, neg_lb = FALSE,
                                               alpha = 0.05, n_cl = 1, verbose = FALSE,
                                               global = FALSE, pairwise = FALSE, 
                                               dunnet = FALSE, trend = FALSE,
                                               iter_control = list(tol = 1e-2, max_iter = 10, 
                                                                   verbose = FALSE),
                                               em_control = list(tol = 1e-5, max_iter = 100),
                                               lme_control = lme4::lmerControl(), 
                                               mdfdr_control = list(fwer_ctrl_method = "BH", B = 100)), 
                                      silent = TRUE))
  if (inherits(out_ANCOMBC, "try-error")) {
    power_ANCOMBC_q <- 0
    fdr_ANCOMBC_q <- 0
  } else {
    id_ANCOMBC_q <- which(out_ANCOMBC[["res"]][["q_X1"]] <= 0.05)
    power_ANCOMBC_q <- length(intersect(id_ANCOMBC_q, id)) / length(id)
    fdr_ANCOMBC_q <- 1 - length(intersect(id_ANCOMBC_q, id)) / length(id_ANCOMBC_q)
  }
  
  
  
  ############ LinDA ############
  suppressWarnings(out_LinDA <- try(linda(feature.dat = counts, meta.dat = smd,
                                          formula = "~ X + (1 | sub)",
                                          alpha = 0.05, 
                                          prev.filter = 0, 
                                          mean.abund.filter = 0,
                                          adaptive = TRUE,
                                          max.abund.filter = 0,
                                          p.adj.method = "BH",
                                          n.cores = 1, 
                                          verbose = FALSE), 
                                    silent = TRUE))
  if (inherits(out_LinDA, "try-error")) {
    power_LinDA_p <- 0
    fdr_LinDA_p <- 0
    power_LinDA_q <- 0
    fdr_LinDA_q <- 0
  } else {
    id_LinDA_q <- which( out_LinDA[["output"]][["X"]][["padj"]] <= 0.05)
    power_LinDA_q <- length(intersect(id_LinDA_q,id))/length(id)
    fdr_LinDA_q <- 1-length(intersect(id_LinDA_q,id))/length(id_LinDA_q)
  }
  
  
  
  ############ ZIGMM ############ correlation = corCAR1(form = ~ time),
  otu_log = log2(count+1)
  suppressWarnings(f <- try(mms(y = otu_log, fixed = ~ X, data = ID,
                                random = ~ 1 | sub, zi_fixed = ~ U, zi_random = NULL,
                                min.p = 0, method = "zig"), 
                            silent = TRUE))
  if (inherits(f, "try-error")) {
    power_ZIGMM_p <- 0
    fdr_ZIGMM_p <- 0
    power_ZIGMM_0 <- 0
    fdr_ZIGMM_0 <- 0
  } else {
    out = fixed(f)$dist
    out_ZIGMM = out[out[,2]=="X", ]
    id_ZIGMM_q <- which(p.adjust(out_ZIGMM$pvalue, "BH") <= 0.05)
    power_ZIGMM_q <- length(intersect(id_ZIGMM_q,id))/length(id)
    fdr_ZIGMM_q <- 1-length(intersect(id_ZIGMM_q,id))/length(id_ZIGMM_q)
    
    out_zi = fixed(f)$zi
    out_ZIGMM_zi = out_zi[out_zi[,2]=="U", ]
    ZIGMM_zi_q <- p.adjust(out_ZIGMM_zi$pvalue, "BH")
    id_ZIGMM_0 <- which(ZIGMM_zi_q<=0.05)
    power_ZIGMM_0 <- length(intersect(id_ZIGMM_0,id0))/length(id0)
    fdr_ZIGMM_0 <- 1-length(intersect(id_ZIGMM_0,id0))/length(id_ZIGMM_0)
  }
  
  
  ############ ZIBR ############
  RA <- count / rowSums(count)
  logistic_cov_matrix = as.matrix(ID$U)
  subject_ind_matrix = as.matrix(ID$sub)
  time_ind_matrix = as.matrix(ID$time)
  beta_cov_matrix = as.matrix(ID$X)
  
  suppressWarnings(out_ZIBR <- try(sapply(1:K, function(k) {
    zibr_fit <- zibr(logistic_cov=logistic_cov_matrix, beta_cov=beta_cov_matrix,
                     Y=as.matrix(RA[,k]), subject_ind=subject_ind_matrix, time_ind=time_ind_matrix)
    c(zibr_fit[["beta_est_table"]][2, 2], zibr_fit[["logistic_est_table"]][2,2])}
  ), silent = TRUE))
  
  if (inherits(out_ZIBR, "try-error")) {
    power_ZIBR_q <- 0
    fdr_ZIBR_q <- 0
    power_ZIBR_0 <- 0
    fdr_ZIBR_0 <- 0
  } else {
    id_ZIBR_q <- which(p.adjust(out_ZIBR[1,], "BH") <= 0.05)
    power_ZIBR_q <- length(intersect(id_ZIBR_q,id))/length(id)
    fdr_ZIBR_q <- 1-length(intersect(id_ZIBR_q,id))/length(id_ZIBR_q)
    id_ZIBR_0 <- which(p.adjust(out_ZIBR[2,], "BH") <= 0.05)
    power_ZIBR_0 <- length(intersect(id_ZIBR_0,id0))/length(id0)
    fdr_ZIBR_0 <- 1-length(intersect(id_ZIBR_0,id0))/length(id_ZIBR_0)
  }
  
  c(power_mbCARE_q, power_ANCOMBC_q, power_LinDA_q, power_ZIGMM_q, power_ZIBR_q,
    fdr_mbCARE_q, fdr_ANCOMBC_q, fdr_LinDA_q, fdr_ZIGMM_q, fdr_ZIBR_q)
}
stopCluster(cl)


res[is.na(res)]=0
res=res*100
rownames(res)=c("power.mbCARE.q", "power.ANCOMBC", "power.LinDA", "power.ZIGMM", "power.ZIBR",
                "fdr.mbCARE.q", "fdr.ANCOMBC", "fdr.LinDA", "fdr.ZIGMM", "fdr.ZIBR")


##########bar plot#########
f1 = function(i){
  index = sim.seed[,i]
  data.frame(metric = rep(c("power", "FDP"), each=5),
             method = rep(c("mbCARE", "ANCOM-BC2", "LinDA", "ZIGMM", "ZIBR"), 2),
             value = rowMeans(res[, index]))
}

power.fdr = foreach(i = 1:4, .combine=rbind) %do% f1(i)

out = data.frame(L = rep(c("5", "10", "20", "50"), each = 10),
                 power.fdr)
colnames(out) = c("L", "metric", "method", "value")
out$L = factor(out$L, levels = c("5", "10", "20", "50"))
out$metric = factor(out$metric, levels =c("power", "FDP"))
out$method = factor(out$method, levels =c("mbCARE", "ANCOM-BC2", "LinDA", "ZIGMM", "ZIBR"))
line = data.frame(metric=factor(c("power", "FDP"), levels =c("power", "FDP")),y=c(NA,5))

p1 = ggplot(out, aes(x = L, y = value, col = method)) + theme_bw()
p1 = p1+ geom_line(aes(group = method, col = method)) +
  geom_point(size = 3, shape = 21, aes(fill = method)) + 
  facet_grid(vars(metric), labeller = labeller(.cols = label_both))+
  geom_hline(data = line, aes(yintercept=y), linetype = "dashed")+ theme(legend.position = "top")+
  labs( y = 'FDP and power (%)', x = expression(italic(L)))+ 
  scale_fill_npg()+scale_color_npg()
p1

