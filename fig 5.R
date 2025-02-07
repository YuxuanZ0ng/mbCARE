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


n <- length(unique(samp_info$sub))
L <- sapply(unique(samp_info$sub), function(i){sum(samp_info$sub == i)})
N <- nrow(count)
X <- matrix(as.numeric(data_subset$FOOD != "BreastOnly"), N, 1)
U <- cbind(rep(1, N), X)
V <- matrix(rep(1, N), ncol = 1)

############ mbCARE ############
out_mbCARE <- mbCARE(count, X, U, V, interval.prob = 0.8, samp_info = samp_info, prev.cut = 0.05)
q_mbCARE_01 <- out_mbCARE[["Theta"]][["alpha01_q"]][,2]
mbCARE01 <- data.frame(class = "Colonization", taxon = names(q_mbCARE_01), "adjusted p-value" = q_mbCARE_01, "effect size" = out_mbCARE[["Theta"]][["alpha01"]][,2])

q_mbCARE_10 <- out_mbCARE[["Theta"]][["alpha10_q"]][,2]
mbCARE10 <- data.frame(class = "Extinction", taxon = names(q_mbCARE_10), "adjusted p-value" = q_mbCARE_10, "effect size" = out_mbCARE[["Theta"]][["alpha10"]][,2])

q_mbCARE <- out_mbCARE[["bias_correction"]][["Delta_q"]]
Delta <- out_mbCARE[["bias_correction"]][["Delta"]]
mbCAREq <- data.frame(class = "Abundance", taxon = names(q_mbCARE_10), "adjusted p-value" = q_mbCARE, "effect size" = Delta)


mbCARE01 <- mbCARE01 %>%
  mutate(log_adj_pval = 1/(adjusted.p.value)) %>%
  mutate(significant = adjusted.p.value < 0.05)

p2 <- ggplot(mbCARE01, aes(x = effect.size, y = log_adj_pval)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("grey", "#b10026")) + # grey for non-significant, red for significant
  geom_text_repel(data = subset(mbCARE01, adjusted.p.value <= 0.05),
                  aes(label = taxon),
                  size = 4, max.overlaps = 10) + # Adjust max.overlaps to limit labels
  labs(
    x = expression(alpha[k]^{"01"}),
    y = expression(frac(1, "Adjusted P-Value"))
  ) +
  theme_classic() +
  xlim(c(-4, 4)) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1/(0.05), linetype = "dashed", color = "grey")  +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")
p2


mbCAREq <- mbCAREq %>%
  mutate(log_adj_pval = -log(adjusted.p.value)) %>%
  mutate(significant = adjusted.p.value < 0.05)

p3 <- ggplot(mbCAREq, aes(x = effect.size, y = log_adj_pval)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("grey", "#b10026")) + # grey for non-significant, red for significant
  geom_text_repel(data = subset(mbCAREq, adjusted.p.value <= 0.05),
                  aes(label = taxon),
                  size = 4, max.overlaps = 10) + # Adjust max.overlaps to limit labels
  labs(
    x = expression(Delta[k]),
    y = expression(-log("Adjusted P-Value"))
  ) +
  theme_classic() +
  xlim(c(-4, 7)) +
  theme(legend.position = "none") +
  geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "grey")  +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")
p3

# no extinct taxa were founded
which(out_mbCARE[["Theta"]][["alpha10_q"]][,2] <= 0.05)


id_mbCARE_01 <- which(out_mbCARE[["Theta"]][["alpha01_q"]][,2] <= 0.05)
names(id_mbCARE_01)

p4 <- p2+p3 + 
  plot_layout(widths = c(1.5, 1))
p4


################# DA ###############################
id_mbCARE <- which(q_mbCARE[,1] <= 0.05)
names(id_mbCARE)


counts <- t(count[,out_mbCARE[["tax_keep"]]])
colnames(counts) <- paste("sample", 1:N)
ID <- data.frame(
  sample = colnames(counts),
  sub = samp_info$sub,  
  time = samp_info$t, 
  X = X
)

## ANCOMBC will report error when test Actinobacillus, so we removed it.
counts <- counts[!(rownames(counts) %in% c("Actinobacillus")),]
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
q_ANCOMBC <- out_ANCOMBC[["res"]][["q_X1"]]
names(q_ANCOMBC) <- rownames(out_ANCOMBC[["feature_table"]])
id_ANCOMBC_q <- which(out_ANCOMBC[["res"]][["q_X1"]] < 0.05)


############ LinDA ############
counts <- t(count[,out_mbCARE[["tax_keep"]]])
colnames(counts) <- paste("sample", 1:N)
ID <- data.frame(
  sample = colnames(counts),
  sub = samp_info$sub,  
  time = samp_info$t, 
  X = X
)
suppressWarnings(out_LinDA <- try(linda(feature.dat = counts, meta.dat = smd,
                                        formula = "~ X + (1 | sub)",
                                        alpha = 0.05, 
                                        prev.filter = 0.05, 
                                        mean.abund.filter = 0,
                                        adaptive = TRUE,
                                        max.abund.filter = 0,
                                        p.adj.method = "BH",
                                        n.cores = 1, 
                                        verbose = FALSE), 
                                  silent = TRUE))
q_LinDA <- out_LinDA[["output"]][["X"]][["padj"]]
names(q_LinDA) <- rownames(out_LinDA[["feature.dat.use"]])
id_LinDA <- which(out_LinDA[["output"]][["X"]][["padj"]] <= 0.05)



############ ZIGMM ############
otu_log <- log2(count[,out_mbCARE[["tax_keep"]]]+1)
suppressWarnings(f <- try(mms(y = otu_log, fixed = ~ X, data = ID,
                              random = ~ 1 | sub, zi_fixed = ~ X, zi_random = NULL,
                              min.p = 0.05, method = "zig"), 
                          silent = TRUE))
out <- fixed(f)$dist
out_ZIGMM <- out[out[,2]=="X", ]
q_ZIGMM <- p.adjust(out_ZIGMM$pvalue, "BH")
names(q_ZIGMM) <- out_ZIGMM$responses
id_ZIGMM_q <- which(p.adjust(out_ZIGMM$pvalue, "BH") <= 0.05)

df_mbCARE <- data.frame(name = rownames(q_mbCARE), q_mbCARE = q_mbCARE)
df_LinDA <- data.frame(name = names(q_LinDA), q_LinDA = q_LinDA)
df_ZIGMM <- data.frame(name = names(q_ZIGMM), q_ZIGMM = q_ZIGMM)
df_ANCOMBC <- data.frame(name = names(q_ANCOMBC), q_ANCOMBC = q_ANCOMBC)

q_sum <- df_mbCARE %>%
  full_join(df_LinDA, by = "name") %>%
  full_join(df_ZIGMM, by = "name") %>%
  full_join(df_ANCOMBC, by = "name")

rownames(q_sum) <- q_sum$name
q_sum <- q_sum[,-1]
q_sum[is.na(q_sum)] <- 1
findings <- apply((q_sum <= 0.05), 2 , as.numeric)
findings <- as.data.frame(findings)
colnames(findings) <- c("mbCARE", "LinDA", "ZIGMM", "ANCOM-BC2")
rownames(findings) <- rownames(q_sum)
UpSetR::upset(findings,
              point.size = 1.8,
              line.size = 1,
              main.bar.color = "#2a83a2",
              sets.bar.color = c("#3b7960","#3b7960","#3b7960","#D55E00"))

