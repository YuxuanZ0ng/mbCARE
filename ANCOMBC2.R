library(doRNG)

################# sim_data  #########################
sim_plnm = function(abn_table, taxa_are_rows = TRUE,
                    prv_cut = 0.1, n, lib_mean, disp) {
  abn_table = as.matrix(abn_table)
  if (!taxa_are_rows) {
    abn_table = t(abn_table)
  }
  prevalence = apply(abn_table, 1, function(x)
    sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
  tax_keep = which(prevalence >= prv_cut)
  txt = paste0("The number of taxa after filtering is: ", length(tax_keep))
  message(txt)
  
  if (length(tax_keep) > 0) {
    abn_table = abn_table[tax_keep, , drop = FALSE]
  } else {
    stop("No taxa remain under the current cutoff", call. = FALSE)
  }
  rel_table = t(t(abn_table)/colSums(abn_table))
  
  n_taxa = nrow(rel_table)
  mean_rel = rowMeans(rel_table)
  cov_rel = cov(t(rel_table))
  overdisp = cov_rel - diag(mean_rel) + outer(mean_rel, mean_rel)
  exp_cov = matrix(1, nrow = n_taxa, ncol = n_taxa) +
    diag(1/mean_rel) %*% overdisp %*% diag(1/mean_rel)
  
  ev = eigen(exp_cov)
  Lambda = ev$values
  if (all(Lambda <= 0)) {
    top_txt = paste0("All eigenvalues are nonpositive \n",
                     "Please use a different dataset")
    stop(top_txt, call. = FALSE)
  }
  Q = ev$vectors
  pos_idx = which(Lambda > 0)
  if (length(pos_idx) == 1) {
    log_Lambda = log(Lambda[pos_idx])
  } else {
    log_Lambda = diag(log(Lambda[pos_idx]))
  }
  cov_est = Q[, pos_idx, drop = FALSE] %*% log_Lambda %*%
    t(Q[, pos_idx, drop = FALSE])
  
  if (!all(eigen(cov_est)$values > 0)) {
    cov_est = as.matrix(Matrix::nearPD(cov_est)$mat)
  }
  
  mean_est = log(mean_rel) - 0.5 * diag(cov_est)
  
  N = rnbinom(n = n, mu = lib_mean, size = disp)
  sim_data = .rplnm(mu = mean_est, sigma = sqrt(cov_est),
                    n = n, N = N)
  return(t(sim_data))
}






############## ANCOMBC2 #########################
ancombc2 = function(data, assay.type = assay_name, assay_name = "counts",
                    rank = tax_level, tax_level = NULL,
                    fix_formula , rand_formula = NULL,
                    p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                    prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                    group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                    alpha = 0.05, n_cl = 1, verbose = FALSE,
                    global = FALSE, pairwise = FALSE,
                    dunnet = FALSE, trend = FALSE,
                    iter_control = list(tol = 1e-2,
                                        max_iter = 20,
                                        verbose = FALSE),
                    em_control = list(tol = 1e-5, max_iter = 100),
                    lme_control = lme4::lmerControl(),
                    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                    trend_control = list(contrast = NULL,
                                         node = NULL,
                                         solver = "ECOS",
                                         B = 100)){
  if (n_cl > 1) {
    cl = parallel::makeCluster(n_cl)
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }
  
  # 1. Data pre-processing
  # Check for aliases
  if (!is.null(assay.type)) {
    assay_name = assay.type
  }
  
  if (!is.null(rank)) {
    tax_level = rank
  }
  
  # TSE data construction
  tse_obj = .tse_construct(data = data, assay_name = assay_name,
                           tax_level = tax_level, phyloseq = NULL)
  tse = tse_obj$tse
  assay_name = tse_obj$assay_name
  tax_level = tse_obj$tax_level
  tse_alt = tse_obj$tse_alt
  
  # Identify taxa with structural zeros
  if (struc_zero) {
    zero_ind = .get_struc_zero(tse = tse, tax_level = tax_level,
                               assay_name = assay_name,
                               alt = TRUE, group = group, neg_lb = neg_lb)
    # Taxa with structural zeros will be removed from ANCOM-BC2 analyses
    tax_idx = apply(zero_ind[, -1], 1, function(x) all(x == FALSE))
    tax_keep = which(tax_idx)
  }else{
    zero_ind = NULL
    # Taxa with structural zeros will be removed from ANCOM-BC2 analyses
    tax_keep = seq(nrow(tse_alt))}
  
  # Filter data by prevalence and library size
  core1 = .data_core(tse = tse, tax_level = tax_level,
                     assay_name = assay_name, alt = FALSE,
                     prv_cut = prv_cut, lib_cut = lib_cut,
                     tax_keep = NULL, samp_keep = NULL)
  O1 = core1$feature_table
  samp_keep = colnames(O1)
  
  core2 = .data_core(tse = tse, tax_level = tax_level,
                     assay_name = assay_name, alt = TRUE,
                     prv_cut = prv_cut, lib_cut = lib_cut,
                     tax_keep = tax_keep, samp_keep = samp_keep)
  O2 = core2$feature_table
  n_tax = nrow(O2)
  tax_name = rownames(O2)
  
  # Metadata and arguments check
  qc = .data_qc(meta_data = core2$meta_data,
                formula = fix_formula, group = group,
                struc_zero = struc_zero, global = global,
                pairwise = pairwise, dunnet = dunnet,
                mdfdr_control = mdfdr_control,
                trend = trend, trend_control = trend_control)
  meta_data = qc$meta_data
  global = qc$global
  pairwise = qc$pairwise
  dunnet = qc$dunnet
  trend = qc$trend
  trend_control = qc$trend_control
  
  # 2. Estimation of the sample-specific biases
  options(na.action = "na.pass") # Keep NA's in rows of x
  x = model.matrix(formula(paste0("~", fix_formula)), data = meta_data)
  options(na.action = "na.omit") # Switch it back
  fix_eff = colnames(x)
  n_fix_eff = length(fix_eff)
  
  if (nrow(O1) < 50) {
    warn_txt = sprintf(paste0("The number of taxa used for estimating ",
                              "sample-specific biases is: ",
                              nrow(O1),
                              "\nA large number of taxa (>50) is required ",
                              "for the consistent estimation of biases"))
    warning(warn_txt, call. = FALSE)
  }
  o1 = log(O1)
  o1[is.infinite(o1)] = NA
  y1 = o1 - rowMeans(o1, na.rm = TRUE)
  
  # Obtain initial estimates
  if (verbose) {
    message("Obtaining initial estimates ...")
  }
  
  if (is.null(rand_formula)) {
    para1 = .iter_mle(x = x, y = y1, meta_data = meta_data,
                      formula = fix_formula,
                      theta = NULL, tol = iter_control$tol,
                      max_iter = iter_control$max_iter,
                      verbose = iter_control$verbose)
  } else {
    para1 = .iter_remle(x = x, y = y1, meta_data = meta_data,
                        fix_formula = fix_formula,
                        rand_formula = rand_formula,
                        lme_control = lme_control,
                        theta = NULL, tol = iter_control$tol,
                        max_iter = iter_control$max_iter,
                        verbose = iter_control$verbose)
  }
  beta1 = para1$beta
  var_hat1 = para1$var_hat
  
  # Apply E-M algorithm
  if (verbose) {
    message("Estimating sample-specific biases ...")
  }
  fun_list = list(.bias_em)
  
  bias1 = foreach(i = seq_len(ncol(beta1)), .combine = rbind) %dorng% {
    output = fun_list[[1]](beta = beta1[, i],
                           var_hat = var_hat1[, i],
                           tol = em_control$tol,
                           max_iter = em_control$max_iter)
  }
  bias1 = data.frame(bias1, row.names = fix_eff, check.names = FALSE)
  colnames(bias1) = c("delta_em", "delta_wls", "var_delta")
  
  delta_em = bias1$delta_em
  delta_wls = bias1$delta_wls
  var_delta = bias1$var_delta
  
  # Obtain the final estimates for sample-specific biases
  beta1 = t(t(beta1) - delta_em)
  theta_hat = matrix(NA, nrow = nrow(y1), ncol = ncol(y1))
  for (i in seq_len(nrow(y1))) {
    theta_hat[i, ] = y1[i, ] - base::rowSums(t(t(x) * beta1[i, ]), na.rm = TRUE)
  }
  theta_hat = colMeans(theta_hat, na.rm = TRUE)
  names(theta_hat) = colnames(y1)
  
  if (any(is.na(theta_hat))) {
    warn_txt = sprintf(paste("Estimation of sampling fractions failed for the following samples:",
                             paste(names(which(is.na(theta_hat))), collapse = ", "),
                             "These samples may have an excessive number of zero values",
                             sep = "\n"))
    warning(warn_txt, call. = FALSE)
  }
  
  # 3. Obtain unbiased estimates
  o2 = log(O2)
  o2[is.infinite(o2)] = NA
  y2 = o2 - rowMeans(o2, na.rm = TRUE)
  y_bias_crt = data.frame(t(t(y2) - theta_hat), check.names = FALSE)
  if (is.null(rand_formula)) {
    para2 = .iter_mle(x = x, y = y2, meta_data = meta_data,
                      formula = fix_formula,
                      theta = theta_hat, tol = iter_control$tol,
                      max_iter = iter_control$max_iter,
                      verbose = iter_control$verbose)
  } else {
    para2 = .iter_remle(x = x, y = y2, meta_data = meta_data,
                        fix_formula = fix_formula,
                        rand_formula = rand_formula,
                        lme_control = lme_control,
                        theta = theta_hat, tol = iter_control$tol,
                        max_iter = iter_control$max_iter,
                        verbose = iter_control$verbose)
  }
  beta_hat = para2$beta
  var_hat = para2$var_hat
  dof = para2$dof
  
  # Account for the variance of delta
  var_hat = sweep(var_hat, 2, var_delta, "+") +
    2 * sqrt(sweep(var_hat, 2, var_delta, "*"))
  
  # Add a small positive constant to stabilize the variance
  if (is.null(s0_perc)) {
    s02 = 0
  } else {
    s02 = apply(var_hat, 2, function(x)
      stats::quantile(x, s0_perc, na.rm = TRUE))
  }
  var_hat = t(t(var_hat) + s02)
  var_hat[is.na(beta_hat)] = NA
  se_hat = sqrt(var_hat)
  vcov_hat = lapply(seq_len(n_tax), function(i) {
    diag(para2$vcov_hat[[i]]) = var_hat[i, ]
    return(para2$vcov_hat[[i]])
  })
  
  # 4. Sensitivity analysis for pseudo-count addition to 0s
  if (pseudo_sens) {
    if (verbose) {
      message("Sensitivity analysis for pseudo-count addition to 0s: ...")
    }
    pseudo_list = seq(0.01, 0.5, 0.01)
    fun_list = list(.get_p)
    
    if (!is.null(group)) {
      all_levels = as.character(unique(meta_data[, group]))
      n_levels = length(all_levels)
    } else {
      n_levels = NULL
    }
    
    # The sensitivity score is calculated as the standard deviation of the
    # negative log p-values derived from bias-corrected abundance regression
    # analyses across different pseudo-count additions.
    pseudo = NULL
    ss_list = foreach(pseudo = pseudo_list) %dorng% {
      count1 = core2$feature_table
      count2 = count1
      count2[count2 == 0] = pseudo
      log_count = log(count2)
      log_resid = log_count - rowMeans(log_count, na.rm = TRUE)
      log_resid_crt = t(t(log_resid) - theta_hat)
      Y = data.frame(t(log_resid_crt), check.names = FALSE)
      
      p_list = lapply(Y, fun_list[[1]], data = meta_data,
                      formula = fix_formula, group = group,
                      n_levels = n_levels, pairwise = pairwise,
                      global = global, trend = trend)
      p = do.call(rbind, p_list)
      data.frame(p, check.names = FALSE)
    }
    
    ss_3d = array(unlist(ss_list), c(dim(ss_list[[1]]), length(ss_list)))
    ss_tab = apply(ss_3d, c(1, 2), function(x) {
      sum(x > alpha)/length(pseudo_list)
    })
    ss_tab = data.frame(taxon = tax_name, ss_tab,  check.names = FALSE)
    rownames(ss_tab) = NULL
    colnames(ss_tab) = c("taxon", colnames(ss_list[[1]]))
    ss_flag = ss_tab
  } else {
    ss_tab = NULL
  }
  
  # 5. Primary results
  if (verbose) {
    message("ANCOM-BC2 primary results ...")
  }
  W = beta_hat/se_hat
  p_hat = 2 * pt(abs(W), df = dof, lower.tail = FALSE)
  p_hat[is.na(p_hat)] = 1
  q_hat = apply(p_hat, 2, function(x) p.adjust(x, method = p_adj_method))
  diff_abn = q_hat <= alpha & !is.na(q_hat)
  
  beta_prim = data.frame(beta_hat, check.names = FALSE)
  se_prim = data.frame(se_hat, check.names = FALSE)
  W_prim = data.frame(W, check.names = FALSE)
  p_prim = data.frame(p_hat, check.names = FALSE)
  q_prim = data.frame(q_hat, check.names = FALSE)
  diff_prim = data.frame(diff_abn, check.names = FALSE)
  colnames(beta_prim) = paste0("lfc_", colnames(beta_hat))
  colnames(se_prim) = paste0("se_", colnames(se_hat))
  colnames(W_prim) = paste0("W_", colnames(W))
  colnames(p_prim) = paste0("p_", colnames(p_hat))
  colnames(q_prim) = paste0("q_", colnames(q_hat))
  colnames(diff_prim) = paste0("diff_", colnames(diff_abn))
  res = do.call("cbind", list(data.frame(taxon = tax_name),
                              beta_prim, se_prim, W_prim,
                              p_prim, q_prim, diff_prim))
  if (pseudo_sens) {
    ss_flag_prim = ss_flag[fix_eff]
    for (col in fix_eff) {
      ss_flag_prim[[col]] = with(ss_flag_prim,
                                 (ss_flag_prim[[col]] == 0 & res[[paste0("diff_", col)]] == TRUE) |
                                   (ss_flag_prim[[col]] == 1 & res[[paste0("diff_", col)]] == FALSE))
    }
    colnames(ss_flag_prim) = paste0("passed_ss_", colnames(ss_flag_prim))
    res = cbind(res, ss_flag_prim)
  }
  rownames(res) = NULL
  
  # 6. Results of global test
  if (global) {
    if (verbose) {
      message("ANCOM-BC2 global test ...")
    }
    if (is.null(rand_formula)) {
      res_global = .ancombc_global_F(x = x, group = group,
                                     beta_hat = beta_hat,
                                     vcov_hat = vcov_hat,
                                     dof = dof,
                                     p_adj_method = p_adj_method,
                                     alpha = alpha)
    } else {
      res_global = .ancombc_global_LRT(full_model = para2$fits,
                                       fix_formula = fix_formula,
                                       rand_formula = rand_formula,
                                       control = lme_control,
                                       x = x, group = group,
                                       y = y_bias_crt,
                                       meta_data = meta_data,
                                       p_adj_method = p_adj_method,
                                       alpha = alpha)
    }
    if (pseudo_sens) {
      ss_flag_global = ss_flag["global"]
      ss_flag_global$global = (ss_flag_global$global == 0 & res_global$diff_abn == TRUE) |
        (ss_flag_global$global == 1 & res_global$diff_abn == FALSE)
      colnames(ss_flag_global) = "passed_ss"
      res_global = cbind(res_global, ss_flag_global)
    }
    rownames(res_global) = NULL
  } else { res_global = NULL }
  
  # 7. Results of multiple pairwise comparisons
  if (pairwise) {
    if (verbose) {
      message("ANCOM-BC2 multiple pairwise comparisons ...")
    }
    res_pair = .ancombc_pair(x = x, group = group,
                             beta_hat = beta_hat,
                             var_hat = var_hat,
                             vcov_hat = vcov_hat,
                             dof = dof,
                             fwer_ctrl_method = mdfdr_control$fwer_ctrl_method,
                             alpha = alpha,
                             full_model = para2$fits,
                             fix_formula = fix_formula,
                             rand_formula = rand_formula,
                             control = lme_control,
                             y = y_bias_crt,
                             meta_data = meta_data)
    beta_pair = data.frame(res_pair$beta, check.names = FALSE)
    se_pair = data.frame(res_pair$se, check.names = FALSE)
    W_pair = data.frame(res_pair$W, check.names = FALSE)
    p_pair = data.frame(res_pair$p_val, check.names = FALSE)
    q_pair = data.frame(res_pair$q_val, check.names = FALSE)
    diff_pair = data.frame(res_pair$diff_abn, check.names = FALSE)
    
    # Directional test summary
    colnames(beta_pair) = paste0("lfc_", colnames(beta_pair))
    colnames(se_pair) = paste0("se_", colnames(se_pair))
    colnames(W_pair) = paste0("W_", colnames(W_pair))
    colnames(p_pair) = paste0("p_", colnames(p_pair))
    colnames(q_pair) = paste0("q_", colnames(q_pair))
    colnames(diff_pair) = paste0("diff_", colnames(diff_pair))
    res_pair = do.call("cbind", list(data.frame(taxon = tax_name),
                                     beta_pair, se_pair, W_pair,
                                     p_pair, q_pair, diff_pair))
    pair_col_name = gsub("lfc_", "", colnames(beta_pair))
    if (pseudo_sens) {
      ss_flag_pair = ss_flag[grepl(group, colnames(ss_flag))]
      colnames(ss_flag_pair) = pair_col_name
      for (col in pair_col_name) {
        ss_flag_pair[[col]] = with(ss_flag_pair,
                                   (ss_flag_pair[[col]] == 0 & res_pair[[paste0("diff_", col)]] == TRUE) |
                                     (ss_flag_pair[[col]] == 1 & res_pair[[paste0("diff_", col)]] == FALSE))
      }
      
      colnames(ss_flag_pair) = paste0("passed_ss_", colnames(ss_flag_pair))
      res_pair = cbind(res_pair, ss_flag_pair)
    }
    rownames(res_pair) = NULL
  } else {
    res_pair = NULL
  }
  
  # 8. Results of Dunnet's type of test
  if (dunnet) {
    if (verbose) {
      message("ANCOM-BC2 multiple pairwise comparisons against the reference group ...")
    }
    res_dunn = .ancombc_dunn(x = x, group = group, beta_hat = beta_hat,
                             var_hat = var_hat, dof = dof,
                             fwer_ctrl_method = mdfdr_control$fwer_ctrl_method,
                             B = mdfdr_control$B, alpha = alpha)
    beta_dunn = data.frame(res_dunn$beta, check.names = FALSE)
    se_dunn = data.frame(res_dunn$se, check.names = FALSE)
    W_dunn = data.frame(res_dunn$W, check.names = FALSE)
    p_dunn = data.frame(res_dunn$p_val, check.names = FALSE)
    q_dunn = data.frame(res_dunn$q_val, check.names = FALSE)
    diff_dunn = data.frame(res_dunn$diff_abn, check.names = FALSE)
    
    # Directional test summary
    colnames(beta_dunn) = paste0("lfc_", colnames(beta_dunn))
    colnames(se_dunn) = paste0("se_", colnames(se_dunn))
    colnames(W_dunn) = paste0("W_", colnames(W_dunn))
    colnames(p_dunn) = paste0("p_", colnames(p_dunn))
    colnames(q_dunn) = paste0("q_", colnames(q_dunn))
    colnames(diff_dunn) = paste0("diff_", colnames(diff_dunn))
    res_dunn = do.call("cbind", list(data.frame(taxon = tax_name),
                                     beta_dunn, se_dunn, W_dunn,
                                     p_dunn, q_dunn, diff_dunn))
    dunn_col_name = gsub("lfc_", "", colnames(beta_dunn))
    if (pseudo_sens) {
      ss_flag_dunn = ss_flag[dunn_col_name]
      for (col in dunn_col_name) {
        ss_flag_dunn[[col]] = with(ss_flag_dunn,
                                   (ss_flag_dunn[[col]] == 0 & res_dunn[[paste0("diff_", col)]] == TRUE) |
                                     (ss_flag_dunn[[col]] == 1 & res_dunn[[paste0("diff_", col)]] == FALSE))
      }
      colnames(ss_flag_dunn) = paste0("passed_ss_", colnames(ss_flag_dunn))
      res_dunn = cbind(res_dunn, ss_flag_dunn)
    }
    rownames(res_dunn) = NULL
  } else {
    res_dunn = NULL
  }
  
  # 9. Results of pattern analysis
  if (trend) {
    if (verbose) {
      message("ANCOM-BC2 pattern analysis ...")
    }
    res_trend = .ancombc_trend(
      x = x, group = group, beta_hat = beta_hat,
      var_hat = var_hat, vcov_hat = vcov_hat,
      p_adj_method = p_adj_method, alpha = alpha,
      trend_control = trend_control)
    beta_trend = res_trend$beta
    se_trend = res_trend$se
    W_trend = res_trend$W
    p_trend = res_trend$p_val
    q_trend = res_trend$q_val
    diff_trend = res_trend$diff_abn
    
    # Directional test summary
    colnames(beta_trend) = paste0("lfc_", colnames(beta_trend))
    colnames(se_trend) = paste0("se_", colnames(se_trend))
    res_trend = cbind(data.frame(taxon = tax_name),
                      beta_trend, se_trend,
                      data.frame(W = W_trend, p_val = p_trend,
                                 q_val = q_trend, diff_abn = diff_trend))
    if (pseudo_sens) {
      ss_flag_trend = ss_flag["trend"]
      ss_flag_trend$trend = (ss_flag_trend$trend == 0 & res_trend$diff_abn == TRUE) |
        (ss_flag_trend$trend == 1 & res_trend$diff_abn == FALSE)
      colnames(ss_flag_trend) = "passed_ss"
      res_trend = cbind(res_trend, ss_flag_trend)
    }
    rownames(res_trend) = NULL
  } else {
    res_trend = NULL
  }
  
  # 10. Outputs
  out = list(feature_table = O2,
             bias_correct_log_table = y_bias_crt,
             ss_tab = ss_tab,
             zero_ind = zero_ind,
             samp_frac = theta_hat,
             delta_em = delta_em,
             delta_wls = delta_wls,
             res = res,
             res_global = res_global,
             res_pair = res_pair,
             res_dunn = res_dunn,
             res_trend = res_trend)
  
  if (n_cl > 1) {
    parallel::stopCluster(cl)
  }
  
  return(out)
}



# Iterative MLE
.iter_mle = function(x, y, meta_data, formula, theta = NULL,
                     tol, max_iter, verbose = FALSE) {
  tax_id = rownames(y)
  n_tax = nrow(y)
  samp_id = colnames(y)
  n_samp = ncol(y)
  fix_eff = colnames(x)
  n_fix_eff = length(fix_eff)
  tformula = formula(paste0("y_crt ~ ", formula))
  
  # Test for over-parameterization
  lm_smoke = stats::lm(formula = tformula,
                       data = data.frame(y_crt = rnorm(n = n_samp), meta_data))
  
  if (any(is.na(lm_smoke$coefficients))) {
    stop_txt = sprintf(paste("Estimation failed for the following covariates:",
                             paste(names(which(is.na(lm_smoke$coefficients))), collapse = ", "),
                             "Please ensure that these covariates do not have missing values and check for multicollinearity before re-estimating the model",
                             sep = "\n"))
    stop(stop_txt, call. = FALSE)
  }
  
  if (lm_smoke$df.residual == 0) {
    stop_txt = sprintf(paste("No residual degrees of freedom! The model is over-parameterized",
                             "Please consider a more parsimonious model",
                             sep = "\n"))
    stop(stop_txt, call. = FALSE)
  }
  
  if (is.null(theta)) {
    # Sampling fractions
    theta = rep(0, n_samp)
    
    # ML fits
    fits = lapply(seq_len(n_tax), function(i) {
      df = data.frame(y_crt = unlist(y[i, ]) - theta, meta_data)
      suppressWarnings(fit <- try(stats::lm(tformula, data = df),
                                  silent = TRUE))
      if (inherits(fit, "try-error")) {fit = NA}
      return(fit)
    })
    
    # Degree of freedom
    dof = NULL
    
    # Degree of freedom
    dof = NULL
    
    # Degree of freedom
    dof = NULL
    
    # Coefficients
    empty_coef = rep(NA, n_fix_eff)
    names(empty_coef) = fix_eff
    beta = lapply(fits, function(i) {
      beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
      coef_i = if (inherits(i, "lm")) {
        stats::coef(i)
      } else {
        empty_coef
      }
      beta_i[match(names(coef_i), fix_eff)] = coef_i
      return(beta_i)
    })
    beta = do.call("rbind", beta)
    
    # Iterative least square
    iterNum = 0
    epsilon = 100
    empty_fitted = rep(NA, n_samp)
    names(empty_fitted) = samp_id
    while (epsilon > tol & iterNum < max_iter) {
      # Updating beta
      fits = lapply(seq_len(n_tax), function(i) {
        df = data.frame(y_crt = unlist(y[i, ]) - theta, meta_data)
        suppressWarnings(fit <- try(stats::lm(tformula, data = df),
                                    silent = TRUE))
        if (inherits(fit, "try-error")) {fit = NA}
        return(fit)
      })
      
      beta_new = lapply(fits, function(i) {
        beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
        coef_i = if (inherits(i, "lm")) {
          stats::coef(i)
        } else {
          empty_coef
        }
        beta_i[match(names(coef_i), fix_eff)] = coef_i
        return(beta_i)
      })
      beta_new = do.call("rbind", beta_new)
      
      # Updating theta
      y_crt_hat = lapply(fits, function(i) {
        y_crt_hat_i = rep(0, n_samp)
        fitted_i = if (inherits(i, "lm")) {
          stats::fitted(i)
        } else {
          empty_fitted
        }
        y_crt_hat_i[match(names(fitted_i), samp_id)] = fitted_i
        return(y_crt_hat_i)
      })
      y_crt_hat = do.call("rbind", y_crt_hat)
      theta_new = colMeans(y - y_crt_hat, na.rm = TRUE)
      
      # Iteration
      epsilon = sqrt(sum((beta_new - beta)^2, na.rm = TRUE) +
                       sum((theta_new - theta)^2, na.rm = TRUE))
      iterNum = iterNum + 1
      beta = beta_new
      theta = theta_new
      
      if (verbose) {
        txt = sprintf(paste0("ML iteration = ", iterNum,
                             ", epsilon = ", signif(epsilon, 2)))
        message(txt)
      }
    }
    
    # Variance-covariance matrices
    y_crt_hat = lapply(fits, function(i) {
      y_crt_hat_i = rep(0, n_samp)
      fitted_i = if (inherits(i, "lm")) {
        stats::fitted(i)
      } else {
        empty_fitted
      }
      y_crt_hat_i[match(names(fitted_i), samp_id)] = fitted_i
      return(y_crt_hat_i)
    })
    y_crt_hat = do.call("rbind", y_crt_hat)
    eps = t(t(y - y_crt_hat) - theta)
    
    XTX_inv = MASS::ginv(t(x[complete.cases(x), ]) %*% x[complete.cases(x), ])
    vcov_hat = vector(mode = "list", length = n_tax)
    var_hat = matrix(NA, nrow = n_tax, ncol = n_fix_eff)
    for (i in seq_len(n_tax)) {
      sigma2_xxT = matrix(0, ncol = n_fix_eff, nrow = n_fix_eff)
      for (j in seq_len(n_samp)) {
        sigma2_xxT_j = eps[i, j]^2 * x[j, ] %*% t(x[j, ])
        sigma2_xxT_j[is.na(sigma2_xxT_j)] = 0.1
        sigma2_xxT = sigma2_xxT + sigma2_xxT_j
      }
      vcov_hat[[i]] = XTX_inv %*% sigma2_xxT %*% XTX_inv
      rownames(vcov_hat[[i]]) = fix_eff
      colnames(vcov_hat[[i]]) = fix_eff
      var_hat[i, ] = diag(vcov_hat[[i]])
    }
  } else {
    # ML fits
    fits = lapply(seq_len(n_tax), function(i) {
      df = data.frame(y_crt = unlist(y[i, ]) - theta, meta_data)
      suppressWarnings(fit <- try(stats::lm(tformula, data = df),
                                  silent = TRUE))
      if (inherits(fit, "try-error")) {fit = NA}
      return(fit)
    })
    
    # Degree of freedom
    dof = vapply(fits, function(i) {
      if (inherits(i, "lm")) {
        summary(i)$df[2]
      } else {
        999L
      }
    }, FUN.VALUE = integer(1))
    dof = matrix(rep(dof, n_fix_eff), ncol = n_fix_eff, byrow = FALSE)
    
    # Degree of freedom
    dof = vapply(fits, function(i) {
      if (inherits(i, "lm")) {
        summary(i)$df[2]
      } else {
        999L
      }
    }, FUN.VALUE = integer(1))
    dof = matrix(rep(dof, n_fix_eff), ncol = n_fix_eff, byrow = FALSE)
    
    # Degree of freedom
    dof = vapply(fits, function(i) {
      if (inherits(i, "lm")) {
        summary(i)$df[2]
      } else {
        999L
      }
    }, FUN.VALUE = integer(1))
    dof = matrix(rep(dof, n_fix_eff), ncol = n_fix_eff, byrow = FALSE)
    
    # Coefficients
    empty_coef = rep(NA, n_fix_eff)
    names(empty_coef) = fix_eff
    
    beta = lapply(fits, function(i) {
      beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
      coef_i = if (inherits(i, "lm")) {
        stats::coef(i)
      } else {
        empty_coef
      }
      beta_i[match(names(coef_i), fix_eff)] = coef_i
      return(beta_i)
    })
    beta = do.call("rbind", beta)
    
    # Variance-covariance matrices
    empty_fitted = rep(NA, n_samp)
    names(empty_fitted) = samp_id
    
    y_crt_hat = lapply(fits, function(i) {
      y_crt_hat_i = rep(0, n_samp)
      fitted_i = if (inherits(i, "lm")) {
        stats::fitted(i)
      } else {
        empty_fitted
      }
      y_crt_hat_i[match(names(fitted_i), samp_id)] = fitted_i
      return(y_crt_hat_i)
    })
    y_crt_hat = do.call("rbind", y_crt_hat)
    eps = t(t(y - y_crt_hat) - theta)
    
    XTX_inv = MASS::ginv(t(x[complete.cases(x), ]) %*% x[complete.cases(x), ])
    vcov_hat = vector(mode = "list", length = n_tax)
    var_hat = matrix(NA, nrow = n_tax, ncol = n_fix_eff)
    for (i in seq_len(n_tax)) {
      sigma2_xxT = matrix(0, ncol = n_fix_eff, nrow = n_fix_eff)
      for (j in seq_len(n_samp)) {
        sigma2_xxT_j = eps[i, j]^2 * x[j, ] %*% t(x[j, ])
        sigma2_xxT_j[is.na(sigma2_xxT_j)] = 0.1
        sigma2_xxT = sigma2_xxT + sigma2_xxT_j
      }
      vcov_hat[[i]] = XTX_inv %*% sigma2_xxT %*% XTX_inv
      rownames(vcov_hat[[i]]) = fix_eff
      colnames(vcov_hat[[i]]) = fix_eff
      var_hat[i, ] = diag(vcov_hat[[i]])
    }
  }
  
  if (!is.null(dof)) {
    colnames(dof) = fix_eff
    rownames(dof) = tax_id
  }
  colnames(beta) = fix_eff
  rownames(beta) = tax_id
  names(theta) = samp_id
  names(vcov_hat) = tax_id
  colnames(var_hat) = fix_eff
  rownames(var_hat) = tax_id
  
  output = list(dof = dof, beta = beta, theta = theta,
                vcov_hat = vcov_hat, var_hat = var_hat)
  return(output)
}




###############   ancombc_bias_correct    ###############
# Iterative REML
.iter_remle = function(x, y, meta_data, fix_formula, rand_formula,
                       lme_control = lme_control, theta = NULL,
                       tol, max_iter, verbose = FALSE) {
  tax_id = rownames(y)
  n_tax = nrow(y)
  samp_id = colnames(y)
  n_samp = ncol(y)
  fix_eff = colnames(x)
  n_fix_eff = length(fix_eff)
  tformula = formula(paste0("y_crt ~ ", fix_formula, "+ ", rand_formula))
  
  # Test for over-parameterization
  lm_smoke = stats::lm(formula = formula(paste0("y ~ ", fix_formula)),
                       data = data.frame(y = rnorm(n = n_samp), meta_data))
  
  if (any(is.na(lm_smoke$coefficients))) {
    stop_txt = sprintf(paste("Estimation failed for the following covariates:",
                             paste(names(which(is.na(lm_smoke$coefficients))), collapse = ", "),
                             "Please ensure that these covariates do not have missing values and check for multicollinearity before re-estimating the model",
                             sep = "\n"))
    stop(stop_txt, call. = FALSE)
  }
  
  if (lm_smoke$df.residual == 0) {
    stop_txt = sprintf(paste("No residual degrees of freedom! The model is over-parameterized",
                             "Please consider a more parsimonious model",
                             sep = "\n"))
    stop(stop_txt, call. = FALSE)
  }
  
  # Test for the fitting of linear mixed-effects model
  tryCatch({
    # Try to run the lmerTest model
    result <- lmerTest::lmer(formula = tformula,
                             data = data.frame(y_crt = y[1, ], meta_data),
                             control = lme_control)
  },
  error = function(e) {
    # This block will be executed if there's an error in the above code
    message <- sprintf(paste("Encountering the error for `lmerTest` package.",
                             "Please try to select one of your taxa and use its raw counts to fix the same linear mixed-effects model using `lmerTest` without the `ANCOMBC` package.",
                             "Load all necessary packages EXCEPT `ANCOMBC`, and see if the error arises due to package incompatibility or other issues.",
                             "The error message from `lmerTest` is as follows:",
                             e$message, sep = "\n"))
    stop(message, call. = FALSE)
  })
  
  # Estimate sample-specific biases
  if (is.null(theta)) {
    # Initial values
    theta = rep(0, n_samp)
    
    # REML fits
    fits = lapply(seq_len(n_tax), function(i) {
      df = data.frame(y_crt = unlist(y[i, ]) - theta, meta_data)
      fit = tryCatch(
        {
          suppressWarnings(suppressMessages(
            lmerTest::lmer(tformula, data = df, control = lme_control)
          ))
        },
        error = function(e) {
          NA
        }
      )
      return(fit)
    })
    
    # Degree of freedom
    dof = NULL
    
    # Coefficients
    empty_coef = rep(NA, n_fix_eff)
    names(empty_coef) = fix_eff
    beta = lapply(fits, function(i) {
      beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
      if (inherits(i, "lmerModLmerTest")) {
        summ_i = summary(i)
        coef_i = summ_i$coefficients[, "Estimate"]
      } else {
        coef_i = empty_coef
      }
      beta_i[match(names(coef_i), fix_eff)] = coef_i
      return(beta_i)
    })
    beta = do.call("rbind", beta)
    
    # Iterative REML
    iterNum = 0
    epsilon = 100
    empty_fitted = rep(NA, n_samp)
    names(empty_fitted) = samp_id
    while (epsilon > tol & iterNum < max_iter) {
      # Updating beta
      fits = lapply(seq_len(n_tax), function(i) {
        df = data.frame(y_crt = unlist(y[i, ]) - theta, meta_data)
        fit = tryCatch(
          {
            suppressWarnings(suppressMessages(
              lmerTest::lmer(tformula,
                             data = df,
                             control = lme_control)
            ))
          },
          error = function(e) {
            NA
          }
        )
        return(fit)
      })
      
      beta_new = lapply(fits, function(i) {
        beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
        if (inherits(i, "lmerModLmerTest")) {
          summ_i = summary(i)
          coef_i = summ_i$coefficients[, "Estimate"]
        } else {
          coef_i = empty_coef
        }
        beta_i[match(names(coef_i), fix_eff)] = coef_i
        return(beta_i)
      })
      beta_new = do.call("rbind", beta_new)
      
      # Updating theta
      y_crt_hat = lapply(fits, function(i) {
        y_crt_hat_i = rep(0, n_samp)
        fitted_i = if (inherits(i, "lmerModLmerTest")) {
          stats::fitted(i)
        } else {
          empty_fitted
        }
        y_crt_hat_i[match(names(fitted_i), samp_id)] = fitted_i
        return(y_crt_hat_i)
      })
      y_crt_hat = do.call("rbind", y_crt_hat)
      theta_new = colMeans(y - y_crt_hat, na.rm = TRUE)
      
      # Iteration
      epsilon = sqrt(sum((beta_new - beta)^2, na.rm = TRUE) +
                       sum((theta_new - theta)^2, na.rm = TRUE))
      iterNum = iterNum + 1
      beta = beta_new
      theta = theta_new
      
      if (verbose) {
        txt = sprintf(paste0("REML iteration = ", iterNum,
                             ", epsilon = ", signif(epsilon, 2)))
        message(txt)
      }
    }
    
    # Residuals
    empty_resid = rep(NA, n_samp)
    names(empty_resid) = samp_id
    eps = lapply(fits, function(i) {
      eps_i = rep(0, n_samp)
      if (inherits(i, "lmerModLmerTest")) {
        summ_i = summary(i)
        resid_i = summ_i$residuals
      } else {
        resid_i = empty_resid
      }
      eps_i[match(names(resid_i), samp_id)] = resid_i
      return(eps_i)
    })
    eps = do.call("rbind", eps)
    
    # Variance-covariance matrices
    empty_vcov = matrix(NA, nrow = n_fix_eff, ncol = n_fix_eff)
    colnames(empty_vcov) = fix_eff
    rownames(empty_vcov) = fix_eff
    vcov_hat = lapply(fits, function(i) {
      Sigma_hat_i = diag(0.1, nrow = n_fix_eff)
      colnames(Sigma_hat_i) = fix_eff
      rownames(Sigma_hat_i) = fix_eff
      if (inherits(i, "lmerModLmerTest")) {
        summ_i = summary(i)
        vcov_hat_i = as.matrix(summ_i$vcov)
      } else {
        vcov_hat_i = empty_vcov
      }
      Sigma_hat_i[match(rownames(vcov_hat_i), fix_eff),
                  match(colnames(vcov_hat_i), fix_eff)] = vcov_hat_i
      return(Sigma_hat_i)
    })
    var_hat = lapply(vcov_hat, function(i) {
      return(diag(i))
    })
    var_hat = do.call("rbind", var_hat)
  } else {
    # REML fits
    fits = lapply(seq_len(n_tax), function(i) {
      df = data.frame(y_crt = unlist(y[i, ]) - theta, meta_data)
      fit = tryCatch(
        {
          suppressWarnings(suppressMessages(
            lmerTest::lmer(tformula, data = df, control = lme_control)
          ))
        },
        error = function(e) {
          NA
        }
      )
      return(fit)
    })
    
    # Degree of freedom
    dof = lapply(fits, function(i) {
      if (inherits(i, "lmerModLmerTest")) {
        summary(i)$coefficients[, "df"]
      } else {
        rep(999, n_fix_eff)
      }
    })
    dof = do.call("rbind", dof)
    
    # Coefficients
    empty_coef = rep(NA, n_fix_eff)
    names(empty_coef) = fix_eff
    beta = lapply(fits, function(i) {
      beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
      if (inherits(i, "lmerModLmerTest")) {
        summ_i = summary(i)
        coef_i = summ_i$coefficients[, "Estimate"]
      } else {
        coef_i = empty_coef
      }
      beta_i[match(names(coef_i), fix_eff)] = coef_i
      return(beta_i)
    })
    beta = do.call("rbind", beta)
    
    # Residuals
    empty_resid = rep(NA, n_samp)
    names(empty_resid) = samp_id
    eps = lapply(fits, function(i) {
      eps_i = rep(0, n_samp)
      if (inherits(i, "lmerModLmerTest")) {
        summ_i = summary(i)
        resid_i = summ_i$residuals
      } else {
        resid_i = empty_resid
      }
      eps_i[match(names(resid_i), samp_id)] = resid_i
      return(eps_i)
    })
    eps = do.call("rbind", eps)
    
    # Variance-covariance matrices
    empty_vcov = matrix(NA, nrow = n_fix_eff, ncol = n_fix_eff)
    colnames(empty_vcov) = fix_eff
    rownames(empty_vcov) = fix_eff
    vcov_hat = lapply(fits, function(i) {
      Sigma_hat_i = diag(0.1, nrow = n_fix_eff)
      colnames(Sigma_hat_i) = fix_eff
      rownames(Sigma_hat_i) = fix_eff
      if (inherits(i, "lmerModLmerTest")) {
        summ_i = summary(i)
        vcov_hat_i = as.matrix(summ_i$vcov)
      } else {
        vcov_hat_i = empty_vcov
      }
      Sigma_hat_i[match(rownames(vcov_hat_i), fix_eff),
                  match(colnames(vcov_hat_i), fix_eff)] = vcov_hat_i
      return(Sigma_hat_i)
    })
    var_hat = lapply(vcov_hat, function(i) {
      return(diag(i))
    })
    var_hat = do.call("rbind", var_hat)
  }
  
  if (!is.null(dof)) {
    colnames(dof) = fix_eff
    rownames(dof) = tax_id
  }
  colnames(beta) = fix_eff
  rownames(beta) = tax_id
  names(theta) = samp_id
  names(vcov_hat) = tax_id
  colnames(var_hat) = fix_eff
  rownames(var_hat) = tax_id
  rownames(eps) = tax_id
  
  output = list(fits = fits, beta = beta, theta = theta, eps = eps,
                dof = dof, vcov_hat = vcov_hat, var_hat = var_hat)
  return(output)
}

# E-M algorithm
.bias_em = function(beta, var_hat, tol, max_iter) {
  beta = beta[!is.na(beta)]
  nu0 = var_hat
  nu0 = nu0[!is.na(nu0)]
  
  if (any(nu0 == 0)) {
    stop_txt = sprintf(paste("Zero variances have been detected for the following taxa:",
                             paste(names(which(nu0 == 0)), collapse = ", "),
                             "Please remove these taxa or select a more parsimonious model",
                             sep = "\n"))
    stop(stop_txt, call. = FALSE)
  }
  
  # Initials
  pi0_0 = 0.75
  pi1_0 = 0.125
  pi2_0 = 0.125
  delta_0 = mean(beta[beta >= quantile(beta, 0.25, na.rm = TRUE)&
                        beta <= quantile(beta, 0.75, na.rm = TRUE)],
                 na.rm = TRUE)
  if(is.na(delta_0)) delta_0 = mean(beta, na.rm = TRUE)
  l1_0 = mean(beta[beta < quantile(beta, 0.125, na.rm = TRUE)],
              na.rm = TRUE)
  if(is.na(l1_0)) l1_0 = min(beta, na.rm = TRUE)
  l2_0 = mean(beta[beta > quantile(beta, 0.875, na.rm = TRUE)],
              na.rm = TRUE)
  if(is.na(l2_0)) l2_0 = max(beta, na.rm = TRUE)
  kappa1_0 = var(beta[beta < quantile(beta, 0.125, na.rm = TRUE)],
                 na.rm = TRUE)
  if(is.na(kappa1_0)|kappa1_0 == 0) kappa1_0 = 1
  kappa2_0 = var(beta[beta > quantile(beta, 0.875, na.rm = TRUE)],
                 na.rm = TRUE)
  if(is.na(kappa2_0)|kappa2_0 == 0) kappa2_0 = 1
  
  # Apply E-M algorithm
  # Store all paras in vectors/matrices
  pi0_vec = pi0_0
  pi1_vec = pi1_0
  pi2_vec = pi2_0
  delta_vec = delta_0
  l1_vec = l1_0
  l2_vec = l2_0
  kappa1_vec = kappa1_0
  kappa2_vec = kappa2_0
  n_tax = length(beta)
  
  # E-M iteration
  iterNum = 0
  epsilon = 100
  while (epsilon > tol & iterNum < max_iter) {
    # Current value of paras
    pi0 = pi0_vec[length(pi0_vec)]
    pi1 = pi1_vec[length(pi1_vec)]
    pi2 = pi2_vec[length(pi2_vec)]
    delta = delta_vec[length(delta_vec)]
    l1 = l1_vec[length(l1_vec)]
    l2 = l2_vec[length(l2_vec)]
    kappa1 = kappa1_vec[length(kappa1_vec)]
    kappa2 = kappa2_vec[length(kappa2_vec)]
    
    # E-step
    pdf0 = vapply(seq(n_tax), function(i)
      dnorm(beta[i], delta, sqrt(nu0[i])),
      FUN.VALUE = double(1))
    pdf1 = vapply(seq(n_tax), function(i)
      dnorm(beta[i], delta + l1, sqrt(nu0[i] + kappa1)),
      FUN.VALUE = double(1))
    pdf2 = vapply(seq(n_tax), function(i)
      dnorm(beta[i], delta + l2, sqrt(nu0[i] + kappa2)),
      FUN.VALUE = double(1))
    r0i = pi0*pdf0/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2)
    r0i[is.na(r0i)] = 0
    r1i = pi1*pdf1/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2)
    r1i[is.na(r1i)] = 0
    r2i = pi2*pdf2/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2)
    r2i[is.na(r2i)] = 0
    
    # M-step
    pi0_new = mean(r0i, na.rm = TRUE)
    pi1_new = mean(r1i, na.rm = TRUE)
    pi2_new = mean(r2i, na.rm = TRUE)
    delta_new = sum(r0i*beta/nu0 + r1i*(beta-l1)/(nu0+kappa1) +
                      r2i*(beta-l2)/(nu0+kappa2), na.rm = TRUE)/
      sum(r0i/nu0 + r1i/(nu0+kappa1) + r2i/(nu0+kappa2), na.rm = TRUE)
    l1_new = min(sum(r1i*(beta-delta)/(nu0+kappa1), na.rm = TRUE)/
                   sum(r1i/(nu0+kappa1), na.rm = TRUE), 0)
    if (is.na(l1_new)) l1_new = 0
    l2_new = max(sum(r2i*(beta-delta)/(nu0+kappa2), na.rm = TRUE)/
                   sum(r2i/(nu0+kappa2), na.rm = TRUE), 0)
    if (is.na(l2_new)) l2_new = 0
    
    # Nelder-Mead simplex algorithm for kappa1 and kappa2
    obj_kappa1 = function(x){
      log_pdf = log(vapply(seq(n_tax), function(i)
        dnorm(beta[i], delta+l1, sqrt(nu0[i]+x)),
        FUN.VALUE = double(1)))
      log_pdf[is.infinite(log_pdf)] = 0
      -sum(r1i*log_pdf, na.rm = TRUE)
    }
    kappa1_new = nloptr::neldermead(x0 = kappa1,
                                    fn = obj_kappa1, lower = 0)$par
    
    obj_kappa2 = function(x){
      log_pdf = log(vapply(seq(n_tax), function(i)
        dnorm(beta[i], delta+l2, sqrt(nu0[i]+x)),
        FUN.VALUE = double(1)))
      log_pdf[is.infinite(log_pdf)] = 0
      -sum(r2i*log_pdf, na.rm = TRUE)
    }
    kappa2_new = nloptr::neldermead(x0 = kappa2,
                                    fn = obj_kappa2, lower = 0)$par
    
    # Merge to the paras vectors/matrices
    pi0_vec = c(pi0_vec, pi0_new)
    pi1_vec = c(pi1_vec, pi1_new)
    pi2_vec = c(pi2_vec, pi2_new)
    delta_vec = c(delta_vec, delta_new)
    l1_vec = c(l1_vec, l1_new)
    l2_vec = c(l2_vec, l2_new)
    kappa1_vec = c(kappa1_vec, kappa1_new)
    kappa2_vec = c(kappa2_vec, kappa2_new)
    
    # Calculate the new epsilon
    epsilon = sqrt((pi0_new-pi0)^2 + (pi1_new-pi1)^2 + (pi2_new-pi2)^2 +
                     (delta_new-delta)^2 + (l1_new-l1)^2 + (l2_new-l2)^2 +
                     (kappa1_new-kappa1)^2 + (kappa2_new-kappa2)^2)
    iterNum = iterNum + 1
  }
  
  # The EM estimator of bias
  delta_em = delta_new
  
  # The WLS estimator of bias
  pi1 = pi1_new
  pi2 = pi2_new
  l1 = l1_new
  l2 = l2_new
  kappa1 = kappa1_new
  kappa2 = kappa2_new
  # Cluster 0
  C0 = which(beta >= quantile(beta, pi1, na.rm = TRUE) &
               beta < quantile(beta, 1 - pi2, na.rm = TRUE))
  # Cluster 1
  C1 = which(beta < quantile(beta, pi1, na.rm = TRUE))
  # Cluster 2
  C2 = which(beta >= quantile(beta, 1 - pi2, na.rm = TRUE))
  # Numerator of the WLS estimator
  nu = nu0
  nu[C1] = nu[C1] + kappa1
  nu[C2] = nu[C2] + kappa2
  wls_deno = sum(1 / nu)
  # Denominator of the WLS estimator
  wls_nume = 1 / nu
  wls_nume[C0] = (wls_nume * beta)[C0]
  wls_nume[C1] = (wls_nume * (beta - l1))[C1]
  wls_nume[C2] = (wls_nume * (beta - l2))[C2]
  wls_nume = sum(wls_nume)
  
  delta_wls = wls_nume / wls_deno
  
  # Estimate the variance of bias
  var_delta = 1 / wls_deno
  if (is.na(var_delta)) var_delta = 0
  
  output = c(delta_em = delta_em,
             delta_wls = delta_wls,
             var_delta = var_delta)
}



################# ANCOMBC_mult ######################
# ANCOM-BC global test
.ancombc_global_F = function(x, group, beta_hat, vcov_hat,
                             dof = NULL, p_adj_method, alpha){
  tax_id = rownames(beta_hat)
  n_tax = nrow(beta_hat)
  covariates = colnames(x)
  output = data.frame(matrix(NA, nrow = n_tax, ncol = 5))
  colnames(output) = c("taxon", "W", "p_val", "q_val", "diff_abn")
  output$taxon = tax_id
  
  # Loop over the parameters of interest
  group_ind = grepl(group, covariates)
  beta_hat_sub = beta_hat[, group_ind, drop = FALSE]
  vcov_hat_sub = lapply(vcov_hat, function(x) {
    x = x[group_ind, group_ind, drop = FALSE]
  })
  
  if (is.null(dof)) {
    for (i in seq_len(n_tax)) {
      # Loop over taxa
      beta_hat_sub_i = beta_hat_sub[i, ]
      vcov_hat_sub_i = vcov_hat_sub[[i]]
      A = diag(x = 1, nrow = length(beta_hat_sub_i))
      
      suppressWarnings(W_global <- try(t(A %*% beta_hat_sub_i) %*%
                                         MASS::ginv(A %*% vcov_hat_sub_i %*% t(A)) %*%
                                         (A %*% beta_hat_sub_i),
                                       silent = TRUE))
      
      if (inherits(W_global, "try-error")) {
        output[i, "W"] = NA
        output[i, "p_val"] = 1
      } else {
        p_global = 2 * min(pchisq(W_global, df = length(beta_hat_sub_i),
                                  lower.tail = TRUE),
                           pchisq(W_global, df = length(beta_hat_sub_i),
                                  lower.tail = FALSE))
        output[i, "W"] = W_global
        output[i, "p_val"] = p_global
      }
    }
  } else {
    for (i in seq_len(n_tax)) {
      # Loop over taxa
      beta_hat_sub_i = beta_hat_sub[i, ]
      vcov_hat_sub_i = vcov_hat_sub[[i]]
      dof_i = unique(dof[i, ])
      A = diag(x = 1, nrow = length(beta_hat_sub_i))
      
      suppressWarnings(W_global <- try(t(A %*% beta_hat_sub_i) %*%
                                         MASS::ginv(A %*% vcov_hat_sub_i %*% t(A)) %*%
                                         (A %*% beta_hat_sub_i),
                                       silent = TRUE))
      
      if (inherits(W_global, "try-error")) {
        output[i, "W"] = NA
        output[i, "p_val"] = 1
      } else {
        p_global = 2 * min(pf(W_global,
                              df1 = length(beta_hat_sub_i),
                              df2 = dof_i,
                              lower.tail = TRUE),
                           pf(W_global,
                              df1 = length(beta_hat_sub_i),
                              df2 = dof_i,
                              lower.tail = FALSE))
        output[i, "W"] = W_global
        output[i, "p_val"] = p_global
      }
    }
  }
  
  # Model summary
  q_global = p.adjust(output[, "p_val"], method = p_adj_method)
  q_global[is.na(q_global)] = 1
  diff_global = q_global <= alpha & !is.na(q_global)
  
  output$q_val = q_global
  output$diff_abn = diff_global
  return(output)
}

.ancombc_global_LRT = function(full_model, fix_formula, rand_formula,
                               control, x, group,
                               y, meta_data, p_adj_method, alpha){
  tax_id = rownames(y)
  n_tax = nrow(y)
  covariates = colnames(x)
  output = data.frame(matrix(NA, nrow = n_tax, ncol = 5))
  colnames(output) = c("taxon", "W", "p_val", "q_val", "diff_abn")
  output$taxon = tax_id
  
  # Perform LRT
  reduce_fix_formula = gsub(pattern = paste0(" \\+ ", group),
                            replacement = "", x = fix_formula)
  reduce_formula = formula(paste0("y ~ ",
                                  reduce_fix_formula,
                                  "+ ", rand_formula))
  
  reduced_model = lapply(seq_len(n_tax), function(i) {
    df = data.frame(y = unlist(y[i, ]), meta_data)
    fit = try(suppressMessages(lmerTest::lmer(reduce_formula,
                                              data = df,
                                              control = control)),
              silent = TRUE)
    if (inherits(fit, "try-error")) {fit = NA}
    return(fit)
  })
  
  W_p_global = lapply(seq_len(n_tax), function(i) {
    model_comparison = try(suppressMessages(anova(full_model[[i]], reduced_model[[i]])),
                           silent = TRUE)
    if (inherits(model_comparison, "try-error")) {
      output = c(W = NA, p = 1)
    } else {
      output = c(W = model_comparison$Chisq[2],
                 p = model_comparison$`Pr(>Chisq)`[2])
    }
    return(output)
  })
  W_p_global = do.call("rbind", W_p_global)
  W_global = W_p_global[, "W"]
  p_global = W_p_global[, "p"]
  
  # Model summary
  q_global = p.adjust(p_global, method = p_adj_method)
  q_global[is.na(q_global)] = 1
  diff_global = q_global <= alpha & !is.na(q_global)
  
  output$W = W_global
  output$p_val = p_global
  output$q_val = q_global
  output$diff_abn = diff_global
  return(output)
}

# ANCOM-BC multiple pairwise comparisons
.ancombc_pair = function(x, group, beta_hat, var_hat, vcov_hat, dof,
                         fwer_ctrl_method, alpha, full_model,
                         fix_formula, rand_formula, control, y, meta_data) {
  covariates = colnames(x)
  
  # Subset the parameters of interest
  group_ind = grepl(group, covariates)
  beta_hat_sub = beta_hat[, group_ind, drop = FALSE]
  vcov_hat_sub = lapply(vcov_hat, function(x) {
    x[group_ind, group_ind, drop = FALSE]
  })
  dof_group = dof[, group_ind]
  
  # Run the combination function to obtain pairwise comparison results
  beta_hat_pair = t(apply(beta_hat_sub, 1, function(x)
    .combn_fun(x, fun = base::diff, sep = "_")))
  var_hat_pair = t(vapply(vcov_hat_sub, function(x)
    .combn_fun2(x, fun = .var_diff, sep = "_"),
    FUN.VALUE = double(ncol(beta_hat_pair))))
  rownames(var_hat_pair) = rownames(beta_hat_pair)
  se_hat_pair = sqrt(var_hat_pair)
  W_pair = beta_hat_pair/se_hat_pair
  
  # Obtain p-values and mdFDR adjusted p-values
  p_q_pair = .mdfdr(global_test = "pairwise",
                    W = W_pair,
                    dof = dof_group,
                    fwer_ctrl_method = fwer_ctrl_method,
                    x = x, group = group,
                    beta_hat = beta_hat,
                    vcov_hat = vcov_hat,
                    alpha = alpha,
                    full_model = full_model,
                    fix_formula = fix_formula,
                    rand_formula = rand_formula,
                    control = control,
                    y = y,
                    meta_data = meta_data)
  p_hat_pair = p_q_pair$p_val
  q_hat_pair = p_q_pair$q_val
  diff_pair = ifelse(q_hat_pair <= alpha, TRUE, FALSE)
  
  output = list(beta = beta_hat_pair, se = se_hat_pair,
                W = W_pair, p_val = p_hat_pair,
                q_val = q_hat_pair, diff_abn = diff_pair)
  return(output)
}

# ANCOM-BC Dunnet's type of test
.dunn_global = function(x, group, W, B, dof, p_adj_method, alpha) {
  covariates = colnames(x)
  group_ind = grepl(group, covariates)
  n_group = sum(group_ind)
  n_tax = nrow(W)
  tax_id = rownames(W)
  output = data.frame(matrix(NA, nrow = n_tax, ncol = 5))
  colnames(output) = c("taxon", "W", "p_val", "q_val", "diff_abn")
  output$taxon = tax_id
  
  suppressWarnings(W_global <- apply(W, 1, function(x) max(abs(x), na.rm = TRUE)))
  
  W_global_null = matrix(NA, nrow = n_tax, ncol = B)
  for (b in seq_len(B)) {
    W_null_b = matrix(unlist(apply(dof, seq_len(2), function(df) rt(1, df = df))),
                      nrow = nrow(dof), ncol = ncol(dof))
    W_global_null_b = apply(W_null_b, 1, function(x)
      max(abs(x), na.rm = TRUE))
    W_global_null[, b] = W_global_null_b
  }
  p_global = 1/B * apply(W_global_null > W_global, 1, function(x)
    sum(x, na.rm = TRUE))
  
  q_global = p.adjust(p_global, method = p_adj_method)
  q_global[is.na(q_global)] = 1
  diff_global = q_global <= alpha & !is.na(q_global)
  
  output$W = W_global
  output$p_val = p_global
  output$q_val = q_global
  output$diff_abn = diff_global
  return(output)
}

.ancombc_dunn = function(x, group, beta_hat, var_hat, dof,
                         B, fwer_ctrl_method, alpha) {
  covariates = colnames(x)
  
  # Subset the parameters of interest
  group_ind = grepl(group, covariates)
  beta_hat_dunn = beta_hat[, group_ind]
  var_hat_dunn = var_hat[, group_ind]
  se_hat_dunn = sqrt(var_hat_dunn)
  W_dunn = beta_hat_dunn/se_hat_dunn
  dof_dunn = dof[, group_ind]
  
  # Obtain p-values and mdFDR adjusted p-values
  p_q_dunn = .mdfdr(global_test = "dunnet", W = W_dunn, dof = dof_dunn,
                    fwer_ctrl_method = fwer_ctrl_method,
                    x = x, group = group, B = B, alpha = alpha)
  p_hat_dunn = p_q_dunn$p_val
  q_hat_dunn = p_q_dunn$q_val
  diff_dunn = ifelse(q_hat_dunn <= alpha, TRUE, FALSE)
  
  output = list(beta = beta_hat_dunn, se = se_hat_dunn,
                W = W_dunn, p_val = p_hat_dunn,
                q_val = q_hat_dunn, diff_abn = diff_dunn)
  return(output)
}

# ANCOM-BC pattern analysis
.ancombc_trend = function(x, group, beta_hat, var_hat, vcov_hat,
                          p_adj_method, alpha,
                          trend_control = list(contrast = NULL,
                                               node = NULL,
                                               solver = "ECOS",
                                               B = 100)){
  tax_id = rownames(beta_hat)
  n_tax = nrow(beta_hat)
  covariates = colnames(x)
  
  group_ind = grepl(group, covariates)
  n_group = sum(group_ind)
  beta_hat_sub = beta_hat[, group_ind, drop = FALSE]
  var_hat_sub = var_hat[, group_ind, drop = FALSE]
  vcov_hat_sub = lapply(vcov_hat, function(x) {
    x = x[group_ind, group_ind, drop = FALSE]
  })
  
  contrast = trend_control$contrast
  node = trend_control$node
  solver = trend_control$solver
  B = trend_control$B
  
  n_trend = length(contrast)
  trend_name = names(contrast)
  
  fun_list = list(.constrain_est, .l_infty)
  
  beta_hat_opt_all = foreach(i = seq_len(n_tax), .combine = rbind) %dorng% {
    beta_hat_opt = unlist(lapply(X = contrast,
                                 FUN = fun_list[[1]],
                                 beta_hat = beta_hat_sub[i, ],
                                 vcov_hat = vcov_hat_sub[[i]],
                                 solver = solver))
  }
  
  l = matrix(NA, nrow = n_tax, ncol = n_trend)
  for (i in seq_len(n_trend)) {
    beta_hat_opt_i = beta_hat_opt_all[, grepl(trend_name[i], colnames(beta_hat_opt_all))]
    node_i = node[[i]]
    l[, i] = apply(beta_hat_opt_i, 1, function(x) fun_list[[2]](x, node_i))
  }
  W_trend = apply(l, 1, function(x) max(x, na.rm = TRUE))
  names(W_trend) = tax_id
  opt_trend = apply(l, 1, function(x) trend_name[which.max(x)])
  beta_hat_trend = matrix(NA, nrow = n_tax, ncol = n_group)
  for (i in seq_len(n_tax)) {
    beta_hat_trend[i, ] = beta_hat_opt_all[i, grepl(opt_trend[i], colnames(beta_hat_opt_all))]
  }
  rownames(beta_hat_trend) = tax_id
  colnames(beta_hat_trend) = colnames(beta_hat_sub)
  
  ident_mat = diag(1, nrow = n_group)
  var_hat_sub_dup = var_hat_sub[, rep(seq_len(ncol(var_hat_sub)), n_trend)]
  var_hat_sub_dup[is.na(var_hat_sub_dup)] = 1
  b = NULL
  W_trend_null = foreach(b = seq_len(B), .combine = 'cbind') %dorng%
    {
      set.seed(b)
      beta_null = matrix(rnorm(n_group * n_tax), nrow = n_tax)
      beta_null_opt = t(apply(beta_null, 1, function(x) {
        beta_null_opt_x = unlist(lapply(X = contrast,
                                        FUN = fun_list[[1]],
                                        beta_hat = x,
                                        vcov_hat = ident_mat,
                                        solver = solver))
        return(beta_null_opt_x)
      }))
      beta_null_opt = beta_null_opt * sqrt(var_hat_sub_dup)
      
      beta_null_opt_list = split(beta_null_opt, row(beta_null_opt))
      beta_null_opt_names = colnames(beta_null_opt)
      l_null = lapply(beta_null_opt_list, function(x) {
        l_null_x = unlist(lapply(X = seq_len(n_trend),
                                 FUN = function(j) {
                                   l = fun_list[[2]](beta_opt = x[grepl(trend_name[j], beta_null_opt_names)],
                                                     node = node[[j]])
                                   return(l)
                                 }))
        return(l_null_x)
      })
      l_null = do.call(rbind, l_null)
      W_null = apply(l_null, 1, function(x) max(x, na.rm = TRUE))
    }
  W_trend_null = as.matrix(W_trend_null)
  
  p_trend = 1/B * apply(W_trend_null > W_trend, 1, function(x)
    sum(x, na.rm = TRUE))
  
  q_trend = p.adjust(p_trend, method = p_adj_method)
  q_trend[is.na(q_trend)] = 1
  diff_trend = q_trend <= alpha & !is.na(q_trend)
  
  output = list(beta = beta_hat_trend,
                se = sqrt(var_hat_sub),
                W = W_trend, p_val = p_trend,
                q_val = q_trend, diff_abn = diff_trend)
  return(output)
}



################# ANCOMBC_prep ####################
# Construct TSE object
.tse_construct = function(data, assay_name, tax_level, phyloseq) {
  if (!is.null(data)) {
    # Check data types
    if (!inherits(data, c("phyloseq", "SummarizedExperiment",
                          "TreeSummarizedExperiment"))) {
      stop_txt = paste0("The input data should be in one of the ",
                        "following format: ",
                        "`phyloseq`, `SummarizedExperiment`, ",
                        "`TreeSummarizedExperiment`")
      stop(stop_txt, call. = FALSE)
    }
    
    if (inherits(data, "phyloseq")) {
      tse = mia::makeTreeSummarizedExperimentFromPhyloseq(data)
      assay_name = "counts"
    } else {
      tse = data
      assay_name = assay_name
    }
    
    # Check feature metadata
    if (ncol(SummarizedExperiment::rowData(tse)) == 0) {
      tax_tab = matrix(rownames(tse), ncol = 1)
      rownames(tax_tab) = rownames(tse)
      colnames(tax_tab) = c("Species")
      tax_tab = S4Vectors::DataFrame(tax_tab)
      SummarizedExperiment::rowData(tse) = tax_tab
    }
    
    # Check if agglomeration should be performed
    if (is.null(tax_level)) {
      tax_level = "ASV"
      tse_alt = tse
    } else {
      tse_alt = .merge_features(tse, tax_level)
    }
    SingleCellExperiment::altExp(tse, tax_level) = tse_alt
  } else if (!is.null(phyloseq)) {
    tse = mia::makeTreeSummarizedExperimentFromPhyloseq(phyloseq)
    assay_name = "counts"
    
    if (ncol(SummarizedExperiment::rowData(tse)) == 0) {
      tax_tab = matrix(rownames(tse), ncol = 1)
      rownames(tax_tab) = rownames(tse)
      colnames(tax_tab) = c("Species")
      tax_tab = S4Vectors::DataFrame(tax_tab)
      SummarizedExperiment::rowData(tse) = tax_tab
    }
    
    if (is.null(tax_level)) {
      tax_level = "ASV"
      tse_alt = tse
    } else {
      tse_alt = .merge_features(tse, tax_level)
    }
    SingleCellExperiment::altExp(tse, tax_level) = tse_alt
  } else {
    stop_txt = paste0("The input data are missing. ",
                      "Please specify either `data` or `phyloseq`")
    stop(stop_txt, call. = FALSE)
  }
  
  tse_obj = list(tse = tse, assay_name = assay_name,
                 tax_level = tax_level, tse_alt = tse_alt)
  
  return(tse_obj)
}

# Filter data by prevalence and library size
.data_core = function(tse = tse, tax_level, assay_name = assay_name,
                      alt = FALSE, prv_cut, lib_cut,
                      tax_keep = NULL, samp_keep = NULL) {
  if (alt) {
    tse_alt = SingleCellExperiment::altExp(tse, tax_level)
    feature_table = SummarizedExperiment::assay(tse_alt, assay_name)
    meta_data = SummarizedExperiment::colData(tse_alt)
  } else {
    feature_table = SummarizedExperiment::assay(tse, assay_name)
    meta_data = SummarizedExperiment::colData(tse)
  }
  
  # Discard taxa with prevalences < prv_cut
  if (is.null(tax_keep)) {
    prevalence = apply(feature_table, 1, function(x)
      sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
    tax_keep = which(prevalence >= prv_cut)
  }else if (length(tax_keep) == 0) {
    stop("All taxa contain structural zeros", call. = FALSE)
  } else {
    # Discard taxa with structural zeros
    feature_table = feature_table[tax_keep, , drop = FALSE]
    prevalence = apply(feature_table, 1, function(x)
      sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
    tax_keep = which(prevalence >= prv_cut)
  }
  
  if (length(tax_keep) > 0) {
    feature_table = feature_table[tax_keep, , drop = FALSE]
  } else {
    stop("No taxa remain under the current cutoff", call. = FALSE)
  }
  
  # Discard samples with library sizes < lib_cut
  if (is.null(samp_keep)) {
    lib_size = colSums(feature_table, na.rm = TRUE)
    samp_keep = which(lib_size >= lib_cut)
  }
  if (length(samp_keep) > 0){
    feature_table = feature_table[, samp_keep, drop = FALSE]
    meta_data = meta_data[samp_keep, , drop = FALSE]
  } else {
    stop("No samples remain under the current cutoff", call. = FALSE)
  }
  
  output = list(feature_table = feature_table,
                meta_data = meta_data,
                tax_keep = tax_keep,
                samp_keep = samp_keep)
  return(output)
}

# Metadata and arguments check
.data_qc = function(meta_data, formula, group, struc_zero,
                    global, pairwise, dunnet,
                    mdfdr_control, trend, trend_control) {
  # Drop unused levels
  meta_data[] = lapply(meta_data, function(x)
    if(is.factor(x)) factor(x) else x)
  
  # Check if all covariates specified in the formula are columns in meta_data
  vars = unlist(strsplit(formula, split = "\\s*\\+\\s*"))
  missing_vars = vars[!vars %in% colnames(meta_data)]
  if(length(missing_vars) > 0) {
    stop("The following variables specified are not in the meta data: ",
         paste(missing_vars, collapse = ", "))
  }
  
  # Check the group variable
  if (is.null(group)) {
    if (any(c(global, pairwise, dunnet, trend))) {
      stop_txt = paste0("Group variable is required for the multi-group comparison \n",
                        "`group` is `NULL` while some of the arguments ",
                        "(`global`, `pairwise`, `dunnet`, `trend`) are `TRUE`")
      stop(stop_txt, call. = FALSE)
    }
    if (struc_zero) {
      stop_txt = paste0("Please specify the group variable for detecting structural zeros \n",
                        "Otherwise, set `struc_zero = FALSE` to proceed")
      stop(stop_txt, call. = FALSE)
    }
  } else {
    meta_data[, group] = as.factor(meta_data[, group])
    # Check the number of groups
    n_level = nlevels(meta_data[, group])
    if (n_level < 2) {
      stop("The group variable should have >= 2 categories",
           call. = FALSE)
    } else if (n_level < 3) {
      global = FALSE
      pairwise = FALSE
      dunnet = FALSE
      trend = FALSE
      warn_txt = paste0("The group variable has < 3 categories \n",
                        "The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated")
      warning(warn_txt, call. = FALSE)
    }
    
    # Check the mdfdr setting for pairwise and dunnet's tests
    if (pairwise | dunnet) {
      if (is.null(mdfdr_control)) {
        stop("Please specify `mdfdr_control` for pairwise or dunnet's test",
             call. = FALSE)
      }
    }
    
    # Check contrast matrices and nodes for trend test
    if (trend) {
      if (is.null(trend_control)) {
        stop("Please specify the `trend_control` parameter for the trend test.",
             call. = FALSE)
      }
      if (is.null(trend_control$contrast)) {
        stop("Please specify the contrast matrices for the trend test.",
             call. = FALSE)
      }
      if (is.null(trend_control$node)) {
        stop("Please specify the nodes for the trend test",
             call. = FALSE)
      }
      if (length(trend_control$contrast) != length(trend_control$node)) {
        stop("The number of nodes should match the number of contrast matrices",
             call. = FALSE)
      }
      sq_mat_check = vapply(trend_control$contrast, function(x)
        nrow(x) == ncol(x), FUN.VALUE = logical(1))
      if (any(sq_mat_check == FALSE)) {
        stop("The contrast matrices for the trend test should be square matrices",
             call. = FALSE)
      }
      dim_mat_check = vapply(trend_control$contrast, function(x)
        nrow(x), FUN.VALUE = integer(1))
      if (any(dim_mat_check != (n_level - 1))) {
        stop_txt = paste0("The contrast matrices for the trend test should be square matrices ",
                          "with dimension #group - 1 \n",
                          "The number of groups in current data is: ",
                          n_level)
        
        stop(stop_txt, call. = FALSE)
      }
      
      n_trend = length(trend_control$contrast)
      if (is.null(names(trend_control$contrast))) {
        names(trend_control$contrast) = paste0("trend", seq_len(n_trend))
        names(trend_control$node) = paste0("trend", seq_len(n_trend))
      }
    }
    
    # Check the sample size per group
    size_per_group = tapply(meta_data[, group], meta_data[, group], length)
    if (any(size_per_group < 2)) {
      stop_txt = sprintf(paste("Sample size per group should be >= 2",
                               "Small sample size detected for the following group(s): ",
                               paste(names(size_per_group)[which(size_per_group < 2)], collapse = ", "),
                               sep = "\n"))
      stop(stop_txt, call. = FALSE)
    } else if (any(size_per_group < 5)) {
      warn_txt = sprintf(paste("Small sample size detected for the following group(s): ",
                               paste(names(size_per_group)[which(size_per_group < 5)], collapse = ", "),
                               "Variance estimation would be unstable when the sample size is < 5 per group",
                               sep = "\n"))
      warning(warn_txt, call. = FALSE)
    }
  }
  
  output = list(meta_data = meta_data,
                global = global,
                pairwise = pairwise,
                dunnet = dunnet,
                trend = trend,
                trend_control = trend_control)
  return(output)
}

# Identify structural zeros
.get_struc_zero = function(tse, tax_level, assay_name,
                           alt = FALSE, group, neg_lb) {
  if (alt) {
    tse_alt = SingleCellExperiment::altExp(tse, tax_level)
    feature_table = SummarizedExperiment::assay(tse_alt, assay_name)
    meta_data = SummarizedExperiment::colData(tse_alt)
    tax_name = rownames(tse_alt)
  } else {
    feature_table = SummarizedExperiment::assay(tse, assay_name)
    meta_data = SummarizedExperiment::colData(tse)
    tax_name = rownames(tse)
  }
  group_data = factor(meta_data[, group])
  present_table = as.matrix(feature_table)
  present_table[is.na(present_table)] = 0
  present_table[present_table != 0] = 1
  n_tax = nrow(feature_table)
  n_group = nlevels(group_data)
  
  p_hat = matrix(NA, nrow = n_tax, ncol = n_group)
  rownames(p_hat) = rownames(feature_table)
  colnames(p_hat) = levels(group_data)
  for (i in seq_len(n_tax)) {
    p_hat[i, ] = tapply(present_table[i, ], group_data,
                        function(x) mean(x, na.rm = TRUE))
  }
  
  samp_size = matrix(NA, nrow = n_tax, ncol = n_group)
  rownames(samp_size) = rownames(feature_table)
  colnames(samp_size) = levels(group_data)
  for (i in seq_len(n_tax)) {
    samp_size[i, ] = tapply(as.matrix(feature_table)[i, ], group_data,
                            function(x) length(x[!is.na(x)]))
  }
  
  p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)
  
  output = (p_hat == 0)
  # Shall we classify a taxon as a structural zero by its negative lower bound?
  if (neg_lb) output[p_hat_lo <= 0] = TRUE
  
  output = cbind(tax_name, output)
  colnames(output) = c("taxon",
                       paste0("structural_zero (", group,
                              " = ", colnames(output)[-1], ")"))
  output = data.frame(output, check.names = FALSE, row.names = NULL)
  output[, -1] = apply(output[, -1], 2, as.logical)
  return(output)
}


############## ANCOMBC1 ##############
ancombc = function(data = NULL, assay.type = NULL, assay_name = "counts",
                   rank = NULL, tax_level = NULL, phyloseq = NULL,
                   formula, p_adj_method = "holm", prv_cut = 0.10,
                   lib_cut = 0, group = NULL, struc_zero = FALSE,
                   neg_lb = FALSE, tol = 1e-05, max_iter = 100,
                   conserve = FALSE, alpha = 0.05, global = FALSE,
                   n_cl = 1, verbose = FALSE){
  message("'ancombc' has been fully evolved to 'ancombc2'. \n",
          "Explore the enhanced capabilities of our refined method!")
  
  if (n_cl > 1) {
    cl = parallel::makeCluster(n_cl)
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }
  
  # 1. Data pre-processing
  # Check for aliases
  if (!is.null(assay.type)) {
    assay_name = assay.type
  }
  
  if (!is.null(rank)) {
    tax_level = rank
  }
  
  # TSE data construction
  tse_obj = .tse_construct(data = data, assay_name = assay_name,
                           tax_level = tax_level, phyloseq = phyloseq)
  tse = tse_obj$tse
  assay_name = tse_obj$assay_name
  tax_level = tse_obj$tax_level
  
  # Filter data by prevalence and library size
  core = .data_core(tse = tse, tax_level = tax_level, assay_name = assay_name,
                    alt = TRUE, prv_cut = prv_cut, lib_cut = lib_cut,
                    tax_keep = NULL, samp_keep = NULL)
  feature_table = core$feature_table
  meta_data = core$meta_data
  tax_keep = core$tax_keep
  n_tax = nrow(feature_table)
  tax_name = rownames(feature_table)
  if (n_tax < 10) {
    warn_txt = sprintf(paste("ANCOM-BC results would be unreliable when the number of taxa is too small (e.g. < 10)",
                             "The number of taxa in the current dataset is: ",
                             n_tax, sep = "\n"))
    warning(warn_txt, call. = FALSE)
  }
  
  # Metadata and arguments check
  qc = .data_qc(meta_data = meta_data,
                formula = formula, group = group,
                struc_zero = struc_zero, global = global,
                pairwise = FALSE, dunnet = FALSE,
                mdfdr_control = NULL, trend = FALSE, trend_control = NULL)
  meta_data = qc$meta_data
  global = qc$global
  
  # Add pseudocount (1) and take logarithm.
  y = log(feature_table + 1)
  options(na.action = "na.pass") # Keep NA's in rows of x
  x = model.matrix(formula(paste0("~", formula)), data = meta_data)
  options(na.action = "na.omit") # Switch it back
  covariates = colnames(x)
  n_covariates = length(covariates)
  
  # 2. Identify taxa with structural zeros
  if (struc_zero) {
    if (is.null(group)) {
      stop_txt = paste("Please specify the group variable for",
                       "detecting structural zeros.",
                       "Otherwise, set struc_zero = FALSE to proceed")
      stop(stop_txt, call. = FALSE)
    }
    zero_ind = .get_struc_zero(tse = tse, tax_level = tax_level,
                               assay_name = assay_name,
                               alt = TRUE, group = group, neg_lb = neg_lb)
    zero_ind = zero_ind[tax_keep, ]
    rownames(zero_ind) = NULL
  }else{ zero_ind = NULL }
  
  # 3. Estimation of parameters
  if (verbose) {
    message("Obtaining initial estimates ...")
  }
  para = .iter_mle(x = x, y = y, meta_data = meta_data,
                   formula = formula, theta = NULL, tol = tol,
                   max_iter = max_iter, verbose = FALSE)
  beta = para$beta
  vcov_hat = para$vcov_hat
  var_hat = para$var_hat
  
  # 4. Estimation of the sample-specific bias
  if (verbose) {
    message("Estimating sample-specific biases ...")
  }
  fun_list = list(.bias_em)
  bias = foreach(i = seq_len(ncol(beta)), .combine = rbind) %dorng% {
    output = fun_list[[1]](beta = beta[, i],
                           var_hat = var_hat[, i],
                           tol = tol,
                           max_iter = max_iter)
  }
  bias = data.frame(bias, row.names = covariates, check.names = FALSE)
  delta_em = bias$delta_em
  delta_wls = bias$delta_wls
  var_delta = bias$var_delta
  
  # 5. Obtain coefficients, standard errors, and sampling fractions
  beta_hat = beta
  beta_hat = t(t(beta_hat) - delta_em)
  
  # Account for the variance of delta_hat
  if (conserve) {
    var_hat = sweep(var_hat, 2, var_delta, "+") +
      2 * sqrt(sweep(var_hat, 2, var_delta, "*"))
    
    vcov_hat = lapply(seq_len(n_tax), function(i) {
      diag(vcov_hat[[i]]) = var_hat[i, ]
      return(vcov_hat[[i]])
    })
    se_hat = sqrt(var_hat)
  }else{ se_hat = sqrt(var_hat) }
  
  theta_hat = matrix(NA, nrow = n_tax, ncol = ncol(y))
  for (i in seq_len(n_tax)) {
    theta_hat[i, ] = y[i, ] - x %*% beta_hat[i, ]
  }
  theta_hat = colMeans(theta_hat, na.rm = TRUE)
  
  # 6. Primary results
  if (verbose) {
    message("ANCOM-BC primary results ...")
  }
  W = beta_hat/se_hat
  p = 2 * pnorm(abs(W), mean = 0, sd = 1, lower.tail = FALSE)
  q = apply(p, 2, function(x) p.adjust(x, method = p_adj_method))
  diff_abn = q <= alpha & !is.na(q)
  
  beta_prim = cbind(taxon = data.frame(taxon = tax_name),
                    data.frame(beta_hat, check.names = FALSE,
                               row.names = NULL))
  se_prim = cbind(taxon = data.frame(taxon = tax_name),
                  data.frame(se_hat, check.names = FALSE,
                             row.names = NULL))
  W_prim = cbind(taxon = data.frame(taxon = tax_name),
                 data.frame(W, check.names = FALSE,
                            row.names = NULL))
  p_prim = cbind(taxon = data.frame(taxon = tax_name),
                 data.frame(p, check.names = FALSE,
                            row.names = NULL))
  q_prim = cbind(taxon = data.frame(taxon = tax_name),
                 data.frame(q, check.names = FALSE,
                            row.names = NULL))
  diff_prim = cbind(taxon = data.frame(taxon = tax_name),
                    data.frame(diff_abn, check.names = FALSE,
                               row.names = NULL))
  
  res = list(lfc = beta_prim,
             se = se_prim,
             W = W_prim,
             p_val = p_prim,
             q_val = q_prim,
             diff_abn = diff_prim)
  
  # 7. Global test results
  if (global) {
    if (verbose) {
      message("ANCOM-BC global test ...")
    }
    res_global = .ancombc_global_F(x = x, group = group,
                                   beta_hat = beta_hat,
                                   vcov_hat = vcov_hat,
                                   p_adj_method = p_adj_method,
                                   alpha = alpha)
  } else { res_global = NULL }
  
  # 8. Combine the information of structural zeros
  # Set p/q-values and SEs of structural zeros to be 0s.
  if (struc_zero) {
    if (verbose) {
      txt = paste0("Merge the information of structural zeros ... \n",
                   "Note that taxa with structural zeros will have ",
                   "0 p/q-values and SEs")
      message(txt)
    }
    zero_idx = as.matrix(zero_ind[, -1])
    group_ind = grepl(group, c("taxon", covariates))
    zero_mask = 1 - abs((zero_idx - zero_idx[, 1]))
    zero_mask = zero_mask[, -1, drop = FALSE]
    res$se[, group_ind] = res$se[, group_ind] * zero_mask
    res$p_val[, group_ind] = res$p_val[, group_ind] * zero_mask
    res$q_val[, group_ind] = res$q_val[, group_ind] * zero_mask
    res$diff_abn[, group_ind] = data.frame(res$q_val[, group_ind] <= alpha &
                                             !is.na(res$q_val[, group_ind]),
                                           check.names = FALSE)
    
    # Global test
    if (global) {
      zero_mask = 1 - apply(zero_idx, 1, function(x)
        sum(x) > 0 & sum(x) < ncol(zero_idx))
      res_global[, "p_val"] = res_global[, "p_val"] * zero_mask
      res_global[, "q_val"] = res_global[, "q_val"] * zero_mask
      res_global[, "diff_abn"] = res_global[, "q_val"] <= alpha &
        !is.na(res_global[, "q_val"])
    }
  }
  
  # 9. Outputs
  out = list(feature_table = feature_table, zero_ind = zero_ind,
             samp_frac = theta_hat, delta_em = delta_em,
             delta_wls = delta_wls, res = res, res_global = res_global)
  
  if (n_cl > 1) {
    parallel::stopCluster(cl)
  }
  
  return(out)
}



######################## utils #########################
# The function to extract off-diagonol elements
.odiag = function(x) x[col(x) != row(x)]

# The function to calculate the variance of difference
.var_diff = function(x) {
  sum(diag(x)) - sum(.odiag(x))
}

# The function of combination
.combn_fun = function(x, fun, sep) {
  y = c(x, utils::combn(x, 2, FUN = fun))
  combn_mat = utils::combn(names(x), 2)
  combn_name = paste(combn_mat[2, ], combn_mat[1, ], sep = sep)
  names(y) = c(names(x), combn_name)
  return(y)
}

.combn_fun2 = function(x, fun, sep) {
  combn_mat = utils::combn(colnames(x), 2)
  y = vector(mode = "numeric")
  for (i in seq(ncol(combn_mat))) {
    idx = c(combn_mat[2, i], combn_mat[1, i])
    y = c(y, fun(x[idx, idx]))
  }
  y = c(diag(x), y)
  combn_name = paste(combn_mat[2, ], combn_mat[1, ], sep = sep)
  names(y) = c(colnames(x), combn_name)
  return(y)
}

# The mdFDR correction
.mdfdr = function(global_test = c("pairwise", "dunnet"),
                  W, dof, fwer_ctrl_method, ...) {
  
  input_list = list(...)
  
  # The total number of null hypotheses rejected in the global test
  if (global_test == "pairwise") {
    if (is.null(input_list$rand_formula)) {
      res_screen = .ancombc_global_F(x = input_list$x,
                                     group = input_list$group,
                                     beta_hat = input_list$beta_hat,
                                     vcov_hat = input_list$vcov_hat,
                                     dof = dof,
                                     p_adj_method = "BH",
                                     alpha = input_list$alpha)
    } else {
      res_screen = .ancombc_global_LRT(full_model = input_list$full_model,
                                       fix_formula = input_list$fix_formula,
                                       rand_formula = input_list$rand_formula,
                                       control = input_list$control,
                                       x = input_list$x,
                                       group = input_list$group,
                                       y = input_list$y,
                                       meta_data = input_list$meta_data,
                                       p_adj_method = "BH",
                                       alpha = input_list$alpha)
    }
  } else {
    res_screen = .dunn_global(x = input_list$x, group = input_list$group,
                              W = W,
                              B = input_list$B,
                              dof = dof,
                              p_adj_method = "BH",
                              alpha = input_list$alpha)
  }
  R = sum(res_screen$diff_abn)
  
  # P-values for pairwise tests
  p_val = 2 * (pt(abs(W), df = dof, lower.tail = FALSE))
  
  # Only consider R significant taxa with regards to the global test
  screen_ind = res_screen$diff_abn
  p_val = p_val * screen_ind
  p_val[p_val == 0] = 1
  p_val[is.na(p_val)] = 1
  
  # Adjust pairwise p-values at level of R * alpha / d
  n_tax = nrow(W)
  q_val = t(apply(p_val, 1, function(x)
    p.adjust(x, method = fwer_ctrl_method, n = length(x) * n_tax / R)))
  
  output = list(p_val = p_val, q_val = q_val)
  return(output)
}

# Estimate coefficients under constraints
.constrain_est = function(beta_hat, vcov_hat, contrast, solver) {
  beta_opt = CVXR::Variable(rows = length(beta_hat), cols = 1, name = "beta")
  obj = CVXR::Minimize(CVXR::matrix_frac(beta_opt - beta_hat, vcov_hat))
  cons = suppressMessages(contrast %*% beta_opt >= 0)
  problem = CVXR::Problem(objective = obj, constraints = list(cons))
  
  suppressMessages(result <- try(CVXR::solve(problem, solver = solver),
                                 silent = TRUE))
  
  if (inherits(result, "try-error")) {
    beta_opt = rep(0, length(beta_hat))
  } else {
    beta_opt = as.numeric(result$getValue(beta_opt))
  }
  return(beta_opt)
}

# Compute the l_infty norm for a pattern
.l_infty = function(beta_opt, node) {
  l = max(abs(beta_opt[node]),
          abs(beta_opt[node] - beta_opt[length(beta_opt)]),
          na.rm = TRUE)
  return(l)
}

# Generate random variables from the poisson log-normal distribution
.rplnm = function(mu, sigma, n, N) {
  d = length(mu)
  y = MASS::mvrnorm(n = n, mu = mu, Sigma = sigma)
  x = N * exp(y)
  otu_table = matrix(rpois(n = n * d, lambda = x), nrow = n)
  return(otu_table)
}

# Get the p-values for the sensitivity analysis
.get_p = function(y, data, formula, group, n_levels, pairwise, global, trend) {
  tformula = paste0("y ~ ", formula)
  df = data.frame(y = y, data)
  lm_fit = stats::lm(formula(tformula), data = df)
  summ = summary(lm_fit)
  p_val = summ$coefficients[, "Pr(>|t|)"]
  p_val[p_val == 0] = 2e-16
  names(p_val) = rownames(summ$coefficients)
  
  if (pairwise) {
    mcp_arg = paste0(group, ' = "Tukey"')
    comparison = multcomp::glht(lm_fit, linfct = eval(parse(text = paste0("multcomp::mcp(", mcp_arg, ")"))))
    summ = summary(comparison, test = multcomp::adjusted("none"))
    pair_p_val = summ$test$pvalues
    pair_p_val[pair_p_val == 0] = 2e-16
    names(pair_p_val) = paste0(summ$focus, names(pair_p_val))
    pair_p_val = pair_p_val[-(seq_len(n_levels - 1))]
    p_val = c(p_val, pair_p_val)
  }
  
  if (global | trend) {
    anova_fit = anova(lm_fit)
    group_p_val = anova_fit$`Pr(>F)`[grepl(group, rownames(anova_fit))]
    if (group_p_val == 0) group_p_val = 2e-16
    if (global) p_val = c(p_val, global = group_p_val)
    if (trend) p_val = c(p_val, trend = group_p_val)
  }
  
  return(p_val)
}

# Internal wrappers for mia::agglomerateByRank/mergeRows
.merge_features = function(x, merge.by, ...) {
  # Check if merge.by parameter belongs to taxonomyRanks
  if (is.character(merge.by) && length(merge.by) == 1 && merge.by %in% mia::taxonomyRanks(x)) {
    # Merge using agglomerateByRank
    x = mia::agglomerateByRank(x, rank = merge.by, ...)
  } else {
    # Merge using mia::mergeRows
    f = factor(SummarizedExperiment::rowData(x)[, merge.by])
    x = mia::mergeRows(x, f = f, ...)
  }
  return(x)
}


