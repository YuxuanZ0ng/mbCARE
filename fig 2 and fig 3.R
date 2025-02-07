source("ARTM_sim.R")
source("mbCARE.R")
source("ANCOMBC2.R")
library(doSNOW)
library(foreach)
library(ggplot2)
library(ggsci) 


# Define the parameters
ncores <- 64
n.sim <- 100
n <- 50 # Number of subjects
K <- 50  # Number of taxons
L <- c(5, 10, 20, 50) # Number of time points
pi_DA <- 0.2
pi_trans <- 0.2
covariate <- "TRUE"
confounder = "TRUE"
FC <- 1
sim_total = n.sim * 4
sim.seed = matrix(1:sim_total, n.sim, 4)



cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
DATA <- list()
DATA <- foreach(i = 1:4, .combine='c') %do% {
  L.i <- L[i]
  foreach(j = sim.seed[,i], .packages = 'MASS', .combine='c') %dopar% {
    library(dplyr)
    set.seed(j)
    data <- ARTM_sim(n = n, K = K, L = L.i, pi_DA = pi_DA, 
                     pi_trans = pi_trans, FC = FC, 
                     confounder = confounder, covariate = covariate)
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
  id01 <- which(data[["alpha01"]][,2] != 0)
  id10 <- which(data[["alpha10"]][,2] != 0)
  
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
    power_mbCARE_01 <- 0; fdr_mbCARE_01 <- 0 
    power_mbCARE_10 <- 0; fdr_mbCARE_10 <- 0
  }else{
    id_mbCARE_q <- which(out_mbCARE[["bias_correction"]][["Delta_q"]] <= 0.05)
    power_mbCARE_q <- length(intersect(id_mbCARE_q,id))/length(id)
    fdr_mbCARE_q <- 1-length(intersect(id_mbCARE_q,id))/length(id_mbCARE_q)
    id_mbCARE_01 <- which(out_mbCARE[["Theta"]][["alpha01_q"]][,2] <= 0.05)
    id_mbCARE_10 <- which(out_mbCARE[["Theta"]][["alpha10_q"]][,2] <= 0.05)
    power_mbCARE_01 <- length(intersect(id_mbCARE_01,id01))/length(id01)
    fdr_mbCARE_01 <- 1-length(intersect(id_mbCARE_01,id01))/length(id_mbCARE_01)
    power_mbCARE_10 <- length(intersect(id_mbCARE_10,id10))/length(id10)
    fdr_mbCARE_10 <- 1-length(intersect(id_mbCARE_10,id10))/length(id_mbCARE_10)
  }
  
  
  counts <- t(count)
  colnames(counts) <- paste("sample", 1:N)
  rownames(counts) <- paste("taxon", 1:K)
  ID <- data.frame(
    sample = rownames(count),
    sub = samp_info$sub,  
    time = samp_info$t, 
    X = X,
    V = as.numeric(V[,2]),
    U = as.numeric(U[,2])
  )
  
  
  ############ ZIGMM ############
  otu_log = log2(count+1)
  suppressWarnings(f <- try(mms(y = otu_log, fixed = ~ X + V, data = ID,
                                random = ~ 1 | sub, zi_fixed = ~ U, zi_random = NULL,
                                min.p = 0, method = "zig"), 
                            silent = TRUE))
  if (inherits(f, "try-error")) {
    power_ZIGMM_01 <- 0
    fdr_ZIGMM_01 <- 0
    power_ZIGMM_10 <- 0
    fdr_ZIGMM_10 <- 0
  } else {
    out_zi = fixed(f)$zi
    out_zi_X = as.data.frame(out_zi[out_zi[,2]=="U", c(1,5)])
    id_ZIGMM_q <- which(p.adjust(out_zi_X$pvalue, "BH") <= 0.05)
    power_ZIGMM_01 <- length(intersect(id_ZIGMM_q,id01))/length(id01)
    fdr_ZIGMM_01 <- 1-length(intersect(id_ZIGMM_q,id01))/length(id_ZIGMM_q)
    power_ZIGMM_10 <- length(intersect(id_ZIGMM_q,id10))/length(id10)
    fdr_ZIGMM_10 <- 1-length(intersect(id_ZIGMM_q,id10))/length(id_ZIGMM_q)
  }
  
  ############ ZIBR ############
  RA <- count / rowSums(count)
  logistic_cov_matrix = as.matrix(U[,2])
  subject_ind_matrix = as.matrix(ID$sub)
  time_ind_matrix = as.matrix(ID$time)
  beta_cov_matrix = as.matrix(cbind(ID$X, ID$V))
  
  suppressWarnings(out_ZIBR <- try(sapply(1:K, function(k) {
    zibr_fit <- zibr(logistic_cov=logistic_cov_matrix, beta_cov=beta_cov_matrix,
                     Y=as.matrix(RA[,k]), subject_ind=subject_ind_matrix, time_ind=time_ind_matrix)
    zibr_fit[["logistic_est_table"]][2, 2]}), silent = TRUE))
  
  if (inherits(out_ZIBR, "try-error")) {
    power_ZIBR_01 <- 0
    fdr_ZIBR_01 <- 0
    power_ZIBR_10 <- 0
    fdr_ZIBR_10 <- 0
  } else {
    id_ZIBR_q <- which(p.adjust(out_ZIBR, "BH") <= 0.05)
    power_ZIBR_01 <- length(intersect(id_ZIBR_q,id01))/length(id01)
    fdr_ZIBR_01 <- 1-length(intersect(id_ZIBR_q,id01))/length(id_ZIBR_q)
    power_ZIBR_10 <- length(intersect(id_ZIBR_q,id10))/length(id10)
    fdr_ZIBR_10 <- 1-length(intersect(id_ZIBR_q,id10))/length(id_ZIBR_q)
  }
  
  
  c(power_mbCARE_01, power_mbCARE_10, power_ZIGMM_01, power_ZIGMM_10, power_ZIBR_01, power_ZIBR_10,
    fdr_mbCARE_01, fdr_mbCARE_10, fdr_ZIGMM_01, fdr_ZIGMM_10, fdr_ZIBR_01, fdr_ZIBR_10)
}
stopCluster(cl)


res[is.na(res)]=0
res=res*100
rownames(res)=c("power.mbCARE.01", "power.mbCARE.10", "power.ZIGMM.01", "power.ZIGMM.10", "power.ZIBR.01", "power.ZIBR.10",
                "fdr.mbCARE.01", "fdr.mbCARE.10", "fdr.ZIGMM.01", "fdr.ZIGMM.10", "fdr.ZIBR.01", "fdr.ZIBR.10")


##########line plot#########
f = function(i){
  index = sim.seed[,i]
  data.frame(metric = rep(c("power", "FDP"), each = 6),
             state = rep(c("0-1", "1-0"), 6),
             method = rep(rep(c("mbCARE", "ZIGMM", "ZIBR"), each = 2), 2),
             value = rowMeans(res[, index]))
}

power.fdr = foreach(i = 1:4, .combine=rbind) %do% f(i)

out = data.frame(L = rep(c("5", "10", "20", "50"), each = 12),
                 power.fdr)
colnames(out) = c("L", "metric", "state", "method", "value")
out$L = factor(out$L, levels = c("5", "10", "20", "50"))
out$metric = factor(out$metric, levels =c("power", "FDP"))
out$state = factor(out$state, levels =c("0-1", "1-0"))
out$method = factor(out$method, levels =c("mbCARE", "ZIGMM", "ZIBR"))
line = data.frame(metric=factor(c("power", "FDP"), levels =c("power", "FDP")), y=c(NA,5))

p = ggplot(out, aes(x = L, y = value, col = method)) + theme_bw()
p = p + geom_line(aes(group = method, col = method)) +
  geom_point(size = 3, shape = 21, aes(fill = method)) + 
  facet_grid(metric ~ state)+
  geom_hline(data = line, aes(yintercept=y), linetype = "dashed")+ theme(legend.position = "top")+
  labs( y = 'FDP and power (%)', x = expression(italic(L))) + 
  scale_color_manual(values = c("mbCARE" = "#E64B35FF", "ZIGMM" = "#3C5488FF", "ZIBR" = "#F39B7FFF")) + 
  scale_fill_manual(values = c("mbCARE" = "#E64B35FF", "ZIGMM" = "#3C5488FF", "ZIBR" = "#F39B7FFF"))
p



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
  id01 <- which(data[["alpha01"]][,2] != 0)
  id10 <- which(data[["alpha10"]][,2] != 0)
  
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
    power_mbCARE_01 <- 0; fdr_mbCARE_01 <- 0 
    power_mbCARE_10 <- 0; fdr_mbCARE_10 <- 0
  }else{
    id_mbCARE_q <- which(out_mbCARE[["bias_correction"]][["Delta_q"]] <= 0.05)
    power_mbCARE_q <- length(intersect(id_mbCARE_q,id))/length(id)
    fdr_mbCARE_q <- 1-length(intersect(id_mbCARE_q,id))/length(id_mbCARE_q)
    id_mbCARE_01 <- which(out_mbCARE[["Theta"]][["alpha01_q"]][,2] <= 0.05)
    id_mbCARE_10 <- which(out_mbCARE[["Theta"]][["alpha10_q"]][,2] <= 0.05)
    power_mbCARE_01 <- length(intersect(id_mbCARE_01,id01))/length(id01)
    fdr_mbCARE_01 <- 1-length(intersect(id_mbCARE_01,id01))/length(id_mbCARE_01)
    power_mbCARE_10 <- length(intersect(id_mbCARE_10,id10))/length(id10)
    fdr_mbCARE_10 <- 1-length(intersect(id_mbCARE_10,id10))/length(id_mbCARE_10)
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
    V = as.numeric(V[,2]),
    U = as.numeric(U[,2])
  )
  assays = S4Vectors::SimpleList(counts = counts)
  smd = S4Vectors::DataFrame(ID)
  tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)
  
  suppressWarnings(out_ANCOMBC <- try(ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                                               fix_formula = "X + V", 
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
                                          formula = "~ X + V + (1 | sub)",
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
    id_LinDA_p <- which( out_LinDA[["output"]][["X"]][["pvalue"]] <= 0.05)
    id_LinDA_q <- which( out_LinDA[["output"]][["X"]][["padj"]] <= 0.05)
    power_LinDA_p <- length(intersect(id_LinDA_p,id))/length(id)
    fdr_LinDA_p <- 1-length(intersect(id_LinDA_p,id))/length(id_LinDA_p)
    power_LinDA_q <- length(intersect(id_LinDA_q,id))/length(id)
    fdr_LinDA_q <- 1-length(intersect(id_LinDA_q,id))/length(id_LinDA_q)
  }
  
  
  
  ############ ZIGMM ############
  otu_log = log2(count+1)
  suppressWarnings(f <- try(mms(y = otu_log, fixed = ~ X + V, data = ID,
                                random = ~ 1 | sub, zi_fixed = ~ U, zi_random = NULL,
                                min.p = 0, method = "zig"), 
                            silent = TRUE))
  if (inherits(f, "try-error")) {
    power_ZIGMM_p <- 0
    fdr_ZIGMM_p <- 0
  } else {
    out = fixed(f)$dist
    out_ZIGMM = out[out[,2]=="X", ]
    id_ZIGMM_q <- which(p.adjust(out_ZIGMM$pvalue, "BH") <= 0.05)
    power_ZIGMM_q <- length(intersect(id_ZIGMM_q,id))/length(id)
    fdr_ZIGMM_q <- 1-length(intersect(id_ZIGMM_q,id))/length(id_ZIGMM_q)
  }
  
  ############ ZIBR ############
  RA <- count / rowSums(count)
  logistic_cov_matrix = as.matrix(U[,2])
  subject_ind_matrix = as.matrix(ID$sub)
  time_ind_matrix = as.matrix(ID$time)
  beta_cov_matrix = as.matrix(cbind(ID$X, ID$V))
  
  suppressWarnings(out_ZIBR <- try(sapply(1:K, function(k) {
    zibr_fit <- zibr(logistic_cov=logistic_cov_matrix, beta_cov=beta_cov_matrix,
                     Y=as.matrix(RA[,k]), subject_ind=subject_ind_matrix, time_ind=time_ind_matrix)
    zibr_fit[["beta_est_table"]][2, 2]}), silent = TRUE))
  
  if (inherits(out_ZIBR, "try-error")) {
    power_ZIBR_q <- 0
    fdr_ZIBR_q <- 0
  } else {
    id_ZIBR_q <- which(p.adjust(out_ZIBR, "BH") <= 0.05)
    power_ZIBR_q <- length(intersect(id_ZIBR_q,id))/length(id)
    fdr_ZIBR_q <- 1-length(intersect(id_ZIBR_q,id))/length(id_ZIBR_q)
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
  labs( y = 'FDP and power (%)', x =  expression(italic(L)))+ 
  scale_fill_npg()+scale_color_npg()
p1
