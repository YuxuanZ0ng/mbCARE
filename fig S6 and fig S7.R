library(MASS)
source("mbCARE.R")
source("ANCOMBC2.R")
library(MicrobiomeStat)
library(NBZIMM)
library(ggplot2)
library(ggsci)
library(dplyr)
library(reshape2)
library(ggVennDiagram)
library(UpSetR)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(patchwork)


metadata <- read.csv("BabyGut/BabyGut16S_metadata.csv")
otu <- read.csv("BabyGut/BabyGut16S_otu.csv")
tax <- read.csv("BabyGut/BabyGut16S_tax.csv")
otu1 <- otu[,-1]
rownames(otu1) <- otu[,1]
otu1 <- as.data.frame(t(otu1))
otu1$X <- rownames(otu1)

num_cols <- sapply(otu1, is.numeric)
num_cols <- names(num_cols[num_cols])

otu_tax_merged <- merge(otu1, tax, by = "X")

genus_aggregated <- otu_tax_merged %>%
  group_by(Genus) %>%
  dplyr::summarise(across(all_of(num_cols), ~sum(.x, na.rm = TRUE)))

otu2 <- genus_aggregated[,-1]
rownames(otu2) <- genus_aggregated$Genus
otu2 <- as.data.frame(t(otu2))
otu2$ID <- rownames(otu2)


data_subset <- metadata %>% 
  dplyr::select(ID, InfantID, X.days, FOOD) %>%
  mutate(Diet = ifelse(FOOD == "BreastOnly", "breast only", "not breast only"))


data_subset$X.days <- as.numeric(data_subset$X.days)



data_subset <- data_subset %>%
  arrange(InfantID, X.days)


otu_subset <- otu2 %>%
  filter(ID %in% data_subset$ID)

otu_subset <- otu_subset[match(data_subset$ID, otu_subset$ID), ]

rownames(otu_subset) <- otu_subset$ID
count <- otu_subset[,1:61]
tax_keep <- which(apply(count, 2, function(x) {mean(x != 0, na.rm = TRUE) > 0.05})) #filtering rare taxa
count <- count[, tax_keep]


samp_info <- data.frame(
  ID = data_subset$ID,
  sub = data_subset$InfantID,
  t = data_subset$X.days
)

data.gen=function(count){
  n <- length(unique(samp_info$sub))
  L <- sapply(unique(samp_info$sub), function(i){sum(samp_info$sub == i)})
  N <- nrow(count)
  X <- matrix(rbinom(N, size=1, prob=0.5), N, 1)
  U <- cbind(rep(1, N), matrix(rbinom(N, size=1, prob=0.5), N, 1))
  V <- matrix(rep(1, N), ncol = 1)
  test.data=list(count, samp_info, X, U, V)
  names(test.data)=c("count", "samp_info", "X", "U", "V")
  return(test.data)
}


# Define the parameters
ncores <- 64
n.sim <- 100

cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
DATA <- list()
DATA <- foreach(j = 1:n.sim, .packages = 'MASS', .combine='c') %dopar% {
  set.seed(j)
  data <- data.gen(count)
  list(data)}
stopCluster(cl)



cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
p_mbCARE <- foreach(i = 1:n.sim, .combine = 'cbind', .errorhandling = "pass") %dopar% {
  library(foreach)
  library(dplyr)
  library(doRNG)
  library(MASS)
  
  data=DATA[[i]]
  count <- data[["count"]]
  X <- data[["X"]]
  U <- data[["U"]]
  V <- data[["V"]]
  samp_info <- data[["samp_info"]]
  
  suppressWarnings(out_mbCARE <- try(mbCARE(count, X, U, V, samp_info, prev.cut = 0), silent = FALSE))
  if (inherits(out_mbCARE, "try-error")) {
    p = NA; q = NA; p0 = NA;  q0 = NA;  p01 = NA;  q01 = NA;  p10 = NA; q10 = NA
  }else{
    p = out_mbCARE[["bias_correction"]][["Delta_p"]]
    q = out_mbCARE[["bias_correction"]][["Delta_q"]]
    p0 = out_mbCARE[["Theta"]][["alpha_p"]][,2]
    q0 = out_mbCARE[["Theta"]][["alpha_q"]][,2]
    p01 = out_mbCARE[["Theta"]][["alpha01_p"]][,2]
    q01 = out_mbCARE[["Theta"]][["alpha01_q"]][,2]
    p10 = out_mbCARE[["Theta"]][["alpha10_p"]][,2]
    q10 = out_mbCARE[["Theta"]][["alpha10_q"]][,2]
  }
  
  c(p, q, p0, q0, p01, q01, p10, q10)
}
stopCluster(cl)

View(apply(p_mbCARE < 0.05, 1, mean))


cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
p_ANCOMBC <- foreach(i = 1:n.sim, .combine = 'cbind', .errorhandling = "pass") %dopar% {
  library(doRNG)
  
  data=DATA[[i]]
  samp_info <- data[["samp_info"]]
  count <- data[["count"]]
  X <- data[["X"]]
  N <- dim(count)[1]
  K <- dim(count)[2]
  counts <- t(count)
  colnames(counts) <- paste("sample", 1:N)
  rownames(counts) <- paste("taxon", 1:K)
  ID <- data.frame(
    sample = rownames(count),
    sub = samp_info$sub,  
    time = samp_info$t, 
    X = X
  )
  assays = S4Vectors::SimpleList(counts = counts)
  smd = S4Vectors::DataFrame(ID)
  tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)
  
  ############ ANCOMBC ############
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
    p <- NA
    q <- NA
  } else {
    p <- out_ANCOMBC[["res"]][["p_X1"]]
    q <- out_ANCOMBC[["res"]][["q_X1"]]
  }
  
  c(p, q)
}
stopCluster(cl)



cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
p_LinDA <- foreach(i = 1:n.sim, .combine = 'cbind', .errorhandling = "pass") %dopar% {
  library(doRNG)
  library(MicrobiomeStat)
  
  data=DATA[[i]]
  samp_info <- data[["samp_info"]]
  count <- data[["count"]]
  X <- data[["X"]]
  N <- dim(count)[1]
  K <- dim(count)[2]
  counts <- t(count)
  colnames(counts) <- paste("sample", 1:N)
  rownames(counts) <- paste("taxon", 1:K)
  ID <- data.frame(
    sample = rownames(count),
    sub = samp_info$sub,  
    time = samp_info$t, 
    X = X
  )
  assays = S4Vectors::SimpleList(counts = counts)
  smd = S4Vectors::DataFrame(ID)
  tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)
  
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
    p <- NA
    q <- NA
  } else {
    p <- out_LinDA[["output"]][["X"]][["pvalue"]]
    q <- out_LinDA[["output"]][["X"]][["padj"]]
  }
  
  c(p, q)
}
stopCluster(cl)



cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
p_ZIGMM <- foreach(i = 1:n.sim, .combine = 'cbind', .errorhandling = "pass") %dopar% {
  library(doRNG)
  library(NBZIMM)
  
  data=DATA[[i]]
  samp_info <- data[["samp_info"]]
  count <- data[["count"]]
  X <- data[["X"]]
  U <- data[["U"]]
  N <- dim(count)[1]
  K <- dim(count)[2]
  ID <- data.frame(
    sample = rownames(count),
    sub = samp_info$sub,  
    time = samp_info$t, 
    X = X,
    U = as.numeric(U[,2])
  )
  ############ ZIGMM ############
  otu_log = log2(count+1)
  
  suppressWarnings(f <- try(mms(y = otu_log, fixed = ~ X, data = ID,
                                random = ~ 1 | sub, zi_fixed = ~ U, zi_random = NULL,
                                min.p = 0, method = "zig"), 
                            silent = TRUE))
  if (inherits(f, "try-error")) {
    p <- NA
    q <- NA
    p0 <- NA
    q0 <- NA
  } else {
    out = fixed(f)$dist
    out_ZIGMM = out[out[,2]=="X", ]
    p <- out_ZIGMM$pvalue
    q <- p.adjust(out_ZIGMM$pvalue, "BH")
    out_zi = fixed(f)$zi
    p0 = out_zi[out_zi[,2]=="U", ]$pvalue
    q0 <- p.adjust(p0, "BH")
  }
  
  c(p, q, p0, q0)
  
}
stopCluster(cl)


K <- ncol(count)

Negcontr <- data.frame(rowSums(p_mbCARE[1:(2*K),] < 0.05)/100, rowSums(p_ANCOMBC<0.05)/100,  rowSums(p_LinDA<0.05)/100,
                       rowSums(p_ZIGMM[1:(2*K),] < 0.05)/100)
Negcontr[is.na(Negcontr)] <- 0
colnames(Negcontr) <- c("mbCARE","ANCOM-BC2","LinDA","ZIGMM")
Negcontr_p <- Negcontr[1:K,]
Negcontr_p$taxon <- 1:K
rownames(Negcontr_p) <- colnames(count)
Negcontr_q <- Negcontr[(K+1):(2*K),]
Negcontr_q$taxon <- 1:K
rownames(Negcontr_q) <- colnames(count)


Negcontr_long <- pivot_longer(Negcontr_p, cols = c("mbCARE","ANCOM-BC2","LinDA","ZIGMM"),
                              names_to = "Method", values_to = "Value")


p <- ggplot(Negcontr_long, aes(x = factor(taxon), y = Value)) + 
  geom_point() + 
  facet_wrap(~ Method, scales = "fixed") + ylim(c(0,1))+
  labs(x = "Taxon", y = "type I error") + theme_bw()+ 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())+ 
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") 
p


Negcontr0 <- data.frame(rowSums(p_mbCARE[(4*K+1):(6*K),] < 0.05)/100,  
                        rowSums(p_mbCARE[(6*K+1):(8*K),] < 0.05)/100, rowSums(p_ZIGMM[(2*K+1):(4*K),] < 0.05)/100)
Negcontr0[is.na(Negcontr0)] <- 0
colnames(Negcontr0) <- c("mbCARE01","mbCARE10","ZIGMM0")
Negcontr0_p <- Negcontr0[1:K,]
Negcontr0_p$taxon <- 1:K
rownames(Negcontr0_p) <- colnames(count)
Negcontr0_q <- Negcontr0[(K+1):(2*K),]
Negcontr0_q$taxon <- 1:K
rownames(Negcontr0_q) <- colnames(count)

custom_labels <- c(
  "mbCARE01"  = "mbCARE(alpha[k]^\"01\")",
  "mbCARE10"  = "mbCARE(alpha[k]^\"10\")",
  "ZIGMM0"    = "ZIGMM(alpha[k]^0)"
)


Negcontr0_long <- pivot_longer(Negcontr0_p, cols = c("mbCARE01","mbCARE10","ZIGMM0"), 
                               names_to = "Method", values_to = "Value")


p0 <- ggplot(Negcontr0_long, aes(x = factor(taxon), y = Value)) + 
  geom_point() + 
  facet_wrap(~ Method, scales = "fixed", labeller = labeller(Method = custom_labels, .default = label_parsed), nrow = 2) +  # Specify two rows
  ylim(c(0, 1)) +
  labs(x = "Taxon", y = "Type I Error") + 
  theme_bw() + 
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()
  ) + 
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  theme(
    strip.text = element_text(size = 8)  # Optional: Adjust facet label text size
  )

p0

