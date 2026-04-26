# ShiVa v1.1.0 — Algorithm core from ShiVa_me.R, clean package interface.

# ---- Internal utilities ----

updatevcv = function(ape.tre, new.rates) {
  vv = vcv(updatetree(ape.tre, new.rates))
  return(vv)
}

updatetree = function(ape.tre, new.rates) {
  ape.tre$edge.length = ape.tre$edge.length * new.rates
  return(ape.tre)
}

#' @title OU.vcv
#' @description Generate covariance matrix for OU process.
#' @param tree phylogenetic tree
#' @param alpha selection strength
#' @return Covariance matrix V.
OU.vcv = function(tree, alpha) {
  times = updatevcv(tree, 1)
  V = matrix(nrow = nrow(times), ncol = ncol(times))
  for (i in 1:nrow(V)) {
    for (j in 1:ncol(V)) {
      V[i, j] = 1/(2*alpha) * exp(-2*alpha*(times[i,i] - times[i,j])) *
                (1 - exp(-2*alpha*times[i,j]))
    }
  }
  return(V)
}

#' @title Generate Design Matrix
#' @description Constructs a design matrix for a given phylogenetic tree.
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param type A character string: \code{"simpX"} or \code{"orgX"}.
#' @param alpha Selection strength (only for \code{type = "orgX"}).
#' @return A design matrix X.
#' @export
#' @importFrom igraph graph.edgelist get.shortest.paths
generate_design_matrix = function(tree, type = "simpX", alpha = 0) {
  stopifnot(is.ultrametric(tree, tol = 0.001, option = 2))
  stopifnot(sum(1:length(tree$tip.label) %in% tree$edge[, 1]) == 0)
  nTips = length(tree$tip.label)
  rNode = nTips + 1
  nEdges = Nedge(tree)
  g = graph.edgelist(tree$edge, directed = TRUE)
  X = matrix(0, nTips, nEdges)
  root2tip = get.shortest.paths(g, rNode, to = 1:nTips, mode = "out",
                                 output = "epath")$epath
  stopifnot(all(lapply(root2tip, length) > 0))
  Tval = sum(tree$edge.length[root2tip[[1]]])
  if (type == "orgX") {
    for (i in 1:nTips) {
      lvec = c(0, tree$edge.length[root2tip[[i]]])
      timeVec = Tval - cumsum(lvec)
      timeVec = timeVec[1:length(timeVec) - 1]
      X[i, root2tip[[i]]] = 1 - exp(-alpha * timeVec)
    }
  } else if (type == "simpX") {
    for (i in 1:nTips) {
      X[i, root2tip[[i]]] = 1
    }
  } else stop("Undefined design matrix type")
  return(X)
}

#' @title Soft Thresholding
#' @param z Numeric value.
#' @param lambda Threshold level.
#' @return Thresholded value.
#' @export
soft_thresholding = function(z, lambda) {
  if (z > lambda) return(z - lambda)
  else if (z < -lambda) return(z + lambda)
  else return(0)
}

# ---- Core algorithm (from ShiVa_me.R) ----

#' @title Update Step for Gamma
#' @param gamma_k Current gamma value.
#' @param X_k Column of design matrix.
#' @param Sigma Current covariance matrix.
#' @param r Residual vector.
#' @param lambda2 Penalty parameter.
#' @param t Step size.
#' @param penalty \code{"L1"} or \code{"None"}.
#' @param V Baseline covariance matrix.
#' @param q_k Element of design vector q.
#' @return List with updated gamma_k and Sigma.
#' @export
#' @importFrom psych tr
update_step_gamma = function(gamma_k, X_k, Sigma, r, lambda2, t, penalty, V, q_k) {
  XXT = X_k %*% t(X_k)
  M_k = XXT * V

  chol_result = tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(chol_result)) {
    warning("Cholesky failed: skipping gamma update.")
    return(list(gamma_k = gamma_k, Sigma = Sigma))
  }
  InvSigma = chol2inv(chol_result)
  InvSigmar = InvSigma %*% r

  delta = -0.5 * (t(InvSigmar) %*% M_k %*% InvSigmar) -
          0.5 * q_k * (t(X_k) %*% InvSigmar)^2 +
          0.5 * tr(M_k %*% InvSigma) +
          0.5 * q_k * (t(X_k) %*% InvSigma %*% X_k)
  delta = as.numeric(delta)

  if (abs(delta) < 1e-6 || !is.finite(delta)) {
    return(list(gamma_k = gamma_k, Sigma = Sigma))
  }

  z = gamma_k - t * delta
  new_value = if (penalty == "L1") soft_thresholding(z, lambda2) else z
  new_value = max(0, min(new_value, 10))

  for (i in 1:5) {
    diag_shift = (new_value - gamma_k) * (q_k + V[1, 1])
    if (all(diag(Sigma)[X_k != 0] + diag_shift > 0)) break
    t = t * 0.5
    z = gamma_k - t * delta
    new_value = if (penalty == "L1") soft_thresholding(z, lambda2) else z
    new_value = max(0, min(new_value, 10))
  }

  Sigma = Sigma + (new_value - gamma_k) * (M_k + q_k * XXT)
  return(list(gamma_k = new_value, Sigma = Sigma))
}

#' @title Fit OU Model with Shifts in Mean and Variance
#' @description Refit with known shift locations using \code{optim} (L-BFGS-B).
#' @param tree Phylogenetic tree.
#' @param Y Trait values.
#' @param alpha Selection strength.
#' @param sv_mean Branch indices with mean shifts.
#' @param sv_var Branch indices with variance shifts.
#' @param max.steps Max iterations. Default 1000.
#' @param t Step size (unused; kept for backward compatibility).
#' @param thres Convergence threshold (unused; kept for backward compatibility).
#' @param measurement_error Logical.
#' @param max.num.shifts Maximum number of shifts. Default \code{Inf}.
#' @param verbose Print convergence info. Default FALSE.
#' @return Fitted model list.
#' @export
#' @import ape
#' @importFrom stats optim
fit_OU_mean_var = function(tree, Y, alpha, sv_mean, sv_var,
                            max.steps = 1000, t = 0.01, thres = 0.01,
                            measurement_error = FALSE, max.num.shifts = Inf,
                            verbose = FALSE) {
  X = generate_design_matrix(tree, "simpX")
  n = nrow(X); p = ncol(X)

  max.num.shifts = min(max.num.shifts, p)
  if (length(sv_mean) > max.num.shifts || length(sv_var) > max.num.shifts) {
    return(list(loglik = -Inf, BIC = Inf, mBIC = Inf, pBIC = Inf))
  }

  V = if (alpha == 0) vcv(tree) else OU.vcv(tree, alpha)
  tb = node.depth.edgelength(tree)[tree$edge[, 1]]
  q = exp(-2*alpha*max(node.depth.edgelength(tree))) / (2*alpha) * (1 - exp(2*alpha*tb))

  par0 = c(rep(0, length(sv_mean)), rep(0, length(sv_var)), 0, log(1))
  if (measurement_error) par0 = c(par0, log(1))

  nll_fn = function(par) {
    m = length(sv_mean); v = length(sv_var)
    beta = rep(0, p); gamma = rep(0, p)
    if (m > 0) beta[sv_mean] = par[1:m]
    if (v > 0) gamma[sv_var] = par[(m+1):(m+v)]
    b0 = par[m + v + 1]
    sigma2 = exp(par[m + v + 2])
    sigma2_error = if (measurement_error) exp(par[m + v + 3]) else 0
    mu = X %*% beta + b0
    r = Y - mu
    Sigma = V*sigma2 + (X %*% diag(gamma) %*% t(X))*V + X %*% diag(gamma*q) %*% t(X) +
            sigma2_error * diag(n)
    cholSigma = tryCatch(chol(Sigma), error = function(e) NULL)
    if (is.null(cholSigma)) return(1e10)
    logdet = 2 * sum(log(diag(cholSigma)))
    InvSigma_r = backsolve(cholSigma, forwardsolve(t(cholSigma), r))
    0.5 * (n*log(2*pi) + logdet + sum(r * InvSigma_r))
  }

  opt_result = optim(par = par0, fn = nll_fn, method = "L-BFGS-B",
                     control = list(maxit = max.steps))
  opt_par = opt_result$par
  m = length(sv_mean); v = length(sv_var)

  beta = rep(0, p); gamma = rep(0, p)
  if (m > 0) beta[sv_mean] = opt_par[1:m]
  if (v > 0) gamma[sv_var] = opt_par[(m+1):(m+v)]
  b0 = opt_par[m + v + 1]
  sigma2 = exp(opt_par[m + v + 2])
  sigma2_error = if (measurement_error) exp(opt_par[m + v + 3]) else 0

  mu = X %*% beta + b0
  r = Y - mu
  Sigma = V*sigma2 + (X %*% diag(gamma) %*% t(X))*V + X %*% diag(gamma*q) %*% t(X) +
          sigma2_error * diag(n)
  InvSigma = solve(Sigma, tol = exp(-100))
  loglik = -nll_fn(opt_par)

  k_beta = sum(beta != 0); k_gamma = sum(gamma != 0)
  BIC = -2*loglik + log(n)*(2*k_beta + 2*k_gamma + 3)

  logn_mean = if (length(sv_mean) > 1) sum(log(colSums(X[, sv_mean, drop = FALSE])))
              else if (length(sv_mean) == 1) log(sum(X[, sv_mean])) else 0
  logn_var = if (length(sv_var) > 1) sum(log(colSums(X[, sv_var, drop = FALSE])))
             else if (length(sv_var) == 1) log(sum(X[, sv_var])) else 0
  active = c(sv_mean, sv_var)
  logn0 = if (length(active) == 0) log(n)
          else if (length(active) == 1) log(n - sum(X[, active] > 0))
          else {
            z = rowSums(X[, active, drop = FALSE]) > 0
            if (sum(z) < n) log(n - sum(z)) else 0
          }
  mBIC = -2*loglik + (2*k_beta + 2*k_gamma - 1)*log(n) + logn_mean + logn_var + logn0

  pBIC = tryCatch({
    design_matrix = cbind(1, X[, sv_mean, drop = FALSE])
    Xt_Si_X = t(design_matrix) %*% InvSigma %*% design_matrix
    logdet_proj = as.numeric(determinant(Xt_Si_X, logarithm = TRUE)$modulus)
    -2*loglik + 2*(k_beta + k_gamma)*log(2*n - 3) + 2*log(n) + logdet_proj
  }, error = function(e) Inf)

  if (!is.finite(loglik)) loglik = -Inf
  if (!is.finite(BIC))    BIC    = Inf
  if (!is.finite(mBIC))   mBIC   = Inf
  if (!is.finite(pBIC))   pBIC   = Inf

  result = list(tree = tree, Y = Y, alpha = alpha,
                shifts_mean = which(beta != 0), shifts_var = which(gamma != 0),
                beta = beta, gamma = gamma, sigma2 = sigma2, b0 = b0,
                loglik = loglik, BIC = BIC, mBIC = mBIC, pBIC = pBIC,
                fitted.values = mu, Sigma = Sigma)
  if (measurement_error) result$sigma2_error = sigma2_error
  class(result) <- "ShiftModel"
  return(result)
}

#' @title Estimate Shifts in Mean and Variance (LASSO)
#' @description L1-penalised estimation of shifts.
#' @param Y Trait values.
#' @param tree Phylogenetic tree.
#' @param alpha Selection strength.
#' @param lambda1 Penalty for mean shifts.
#' @param lambda2 Penalty for variance shifts.
#' @param max.steps Max iterations.
#' @param t Step size.
#' @param penalty \code{"L1"} or \code{"None"}.
#' @param thres Convergence threshold.
#' @param sigma2 Optional initial sigma2.
#' @param measurement_error Logical.
#' @param verbose Print convergence info. Default FALSE.
#' @return List with detected shifts and parameter estimates.
#' @export
#' @importFrom glmnet glmnet
get_mean_var_shifts = function(Y, tree, alpha, lambda1, lambda2,
                                max.steps = 1000, t = 0.01, penalty = 'L1',
                                thres = 0.01, sigma2 = NULL,
                                measurement_error = FALSE, verbose = FALSE) {
  X = generate_design_matrix(tree, type = 'simpX')
  internal_list = (1:ncol(X))[colSums(X) > 1 & colSums(X) != (nrow(X) - 1)]
  V = if (alpha == 0) vcv(tree) else OU.vcv(tree, alpha)
  n = nrow(X); p = ncol(X)

  sigma2_0 = if (is.null(sigma2)) 1 else sigma2
  sigma2_error = if (measurement_error) 1 else 0
  tau_0 = log(sigma2_0); tau_error = log(sigma2_error)
  tb = node.depth.edgelength(tree)[tree$edge[, 1]]
  q = exp(-2*alpha*max(node.depth.edgelength(tree))) / (2*alpha) * (1 - exp(2*alpha*tb))

  b0 = 0; beta = rep(0, p); gamma = rep(0, p); r = Y
  loss = Inf; loss_list = c()
  Sigma = V*sigma2_0 + (X %*% diag(gamma) %*% t(X))*V + X %*% diag(gamma*q) %*% t(X)
  InvSigma = solve(Sigma)

  for (s in 1:max.steps) {
    last_tau = tau_0; last_gamma = gamma; last_beta = beta
    last_b0 = b0; last_loss = loss; last_tau_error = tau_error

    # update beta (glmnet)
    svd_result = svd(InvSigma)
    SqrtInvSigma = svd_result$u %*% diag(sqrt(svd_result$d)) %*% t(svd_result$v)
    YY = SqrtInvSigma %*% Y
    XX = SqrtInvSigma %*% X
    beta = as.vector(glmnet(XX, YY, 'gaussian', lambda = lambda1, intercept = FALSE)$beta)
    r = Y - X %*% beta

    # update b0
    tmp = t(rep(1, n)) %*% InvSigma
    b0 = as.vector((tmp %*% (r + b0)) / (tmp %*% rep(1, n)))
    r = r - (b0 - last_b0)

    # update gamma
    for (k in internal_list) {
      result = update_step_gamma(gamma[k], X[, k], Sigma, r, lambda2, t, penalty, V, q[k])
      gamma[k] = result[['gamma_k']]
      Sigma = result[['Sigma']]
    }

    # update sigma0
    if (is.null(sigma2)) {
      InvSigma = solve(Sigma, tol = exp(-100))
      InvSigmar = InvSigma %*% r
      tau_delta = (-1/2*t(InvSigmar) %*% V %*% InvSigmar + 1/2*tr(V %*% InvSigma)) * exp(tau_0)
      tau_0 = (last_tau - t*tau_delta)[1]
      Sigma = Sigma + (exp(tau_0) - exp(last_tau)) * V
      InvSigma = solve(Sigma, tol = exp(-100))
      InvSigmar = InvSigma %*% r
    }

    if (measurement_error) {
      tau_error_delta = (-1/2*t(InvSigmar) %*% InvSigmar + 1/2*tr(InvSigma)) * exp(tau_error)
      tau_error = (last_tau_error - t*tau_error_delta)[1]
      Sigma = Sigma + (exp(tau_error) - exp(last_tau_error)) * diag(1, n)
      InvSigma = solve(Sigma)
    }

    loglik = -1/2*(n*log(2*pi) + t(r) %*% (InvSigma %*% r) + determinant(Sigma)$modulu)
    loss = -loglik
    if (lambda1 != Inf) loss = loss + sum(lambda1*abs(beta))
    if (lambda2 != Inf) loss = loss + sum(lambda2*abs(gamma))
    loss_list = c(loss_list, loss)

    if (is.na(loglik) || is.nan(loglik)) {
      if (verbose) warning("NaN loglik at step ", s)
      break
    }
    if (((abs(loglik) != Inf) & (abs(loss - last_loss) < thres)) | loglik == -Inf) {
      if (verbose) cat(s, "steps to converge\n")
      break
    }
  }

  sigma2_0 = exp(tau_0)
  sigma2_error = exp(tau_error)
  result = list(sv_mean = (1:ncol(X))[beta != 0], sv_var = (1:ncol(X))[gamma != 0],
                gamma = gamma, beta = beta, sigma2 = sigma2_0, b0 = b0)
  if (measurement_error) result$sigma2_error = sigma2_error
  return(result)
}

#' @title Backward Selection for Shift Correction
#' @param tree Phylogenetic tree.
#' @param Y Trait values.
#' @param alpha Selection strength.
#' @param sv_mean Mean shift branches.
#' @param sv_var Variance shift branches.
#' @param criterion \code{"BIC"}, \code{"mBIC"}, or \code{"pBIC"}.
#' @param original_model Optional pre-fitted model.
#' @param measurement_error Logical.
#' @param verbose Logical.
#' @return Refined model.
#' @export
backward_correction = function(tree, Y, alpha, sv_mean, sv_var,
                                criterion = 'BIC', original_model = NULL,
                                measurement_error = FALSE, max.num.shifts = Inf,
                                verbose = FALSE) {
  if (is.null(original_model)) {
    OModel = fit_OU_mean_var(tree, Y, alpha, sv_mean, sv_var,
                              measurement_error = measurement_error,
                              max.num.shifts = max.num.shifts, verbose = verbose)
  } else {
    OModel = original_model
  }
  original_score = OModel[[criterion]][[1]]
  best_score = original_score
  n1 = length(sv_mean); n2 = length(sv_var)
  if (n1 + n2 == 0) return(OModel)

  after_score_list = rep(0, n1 + n2)
  for (i in 1:(n1 + n2)) {
    if (i <= length(sv_mean)) {
      new_Model = fit_OU_mean_var(tree, Y, alpha, setdiff(sv_mean, sv_mean[i]), sv_var,
                                   measurement_error = measurement_error,
                                   max.num.shifts = max.num.shifts, verbose = verbose)
    } else {
      new_Model = fit_OU_mean_var(tree, Y, alpha, sv_mean, setdiff(sv_var, sv_var[i - n1]),
                                   measurement_error = measurement_error,
                                   max.num.shifts = max.num.shifts, verbose = verbose)
    }
    score_i = new_Model[[criterion]][[1]]
    if (!is.finite(score_i)) score_i = Inf
    after_score_list[i] = score_i
    if (isTRUE(score_i < best_score)) {
      OModel = new_Model
      best_score = score_i
    }
  }
  index_list = order(after_score_list)[sort(after_score_list) < original_score]
  if (length(index_list) == 0) return(OModel)
  remove_list = c(index_list[1])
  if (length(index_list) >= 2) {
    for (i in index_list[2:length(index_list)]) {
      new_Model = fit_OU_mean_var(tree, Y, alpha,
                                   setdiff(sv_mean, sv_mean[c(i, remove_list)[c(i, remove_list) <= n1]]),
                                   setdiff(sv_var, sv_var[c(i, remove_list)[c(i, remove_list) > n1] - n1]),
                                   measurement_error = measurement_error,
                                   max.num.shifts = max.num.shifts, verbose = verbose)
      score_new = new_Model[[criterion]][[1]]
      if (!is.finite(score_new)) score_new = Inf
      if (isTRUE(score_new < best_score)) {
        OModel = new_Model
        best_score = score_new
        remove_list = c(remove_list, i)
      }
    }
  }
  return(OModel)
}

#' @title Model Selection for Shifts in Mean and Variance
#' @param Y Trait values.
#' @param tree Phylogenetic tree.
#' @param alpha Selection strength.
#' @param lambda1_list Candidate lambda1 values.
#' @param lambda2_list Candidate lambda2 values.
#' @param criterion \code{"BIC"}, \code{"mBIC"}, or \code{"pBIC"}.
#' @param max.steps Max iterations.
#' @param nfolds CV folds.
#' @param top_k Top candidates for backward correction.
#' @param lambda.type \code{"lambda.1se"} (default, conservative) or \code{"lambda.min"} (less conservative).
#' @param measurement_error Logical.
#' @param verbose Logical. Default FALSE.
#' @return List with best_model and score_summary.
#' @export
#' @importFrom glmnet cv.glmnet
get_mean_var_shifts_model_selection = function(Y, tree, alpha,
                                               lambda1_list = NULL, lambda2_list = exp(1:10 - 6),
                                               criterion = "BIC", max.steps = 300,
                                               nfolds = 5, top_k = 3,
                                               lambda.type = "lambda.1se",
                                               max.num.shifts = Inf,
                                               measurement_error = FALSE, verbose = FALSE) {
  X = generate_design_matrix(tree, 'simpX')
  Y = as.vector(Y)
  V = if (alpha == 0) vcv(tree) else OU.vcv(tree, alpha)
  n = nrow(X); p = ncol(X)
  tb = node.depth.edgelength(tree)[tree$edge[, 1]]
  q = exp(-2*alpha*max(node.depth.edgelength(tree))) / (2*alpha) * (1 - exp(2*alpha*tb))

  score_summary = data.frame(matrix(nrow = 0, ncol = 8))
  names(score_summary) = c('lambda1','lambda2','sv_mean','sv_var','loglik','BIC','mBIC','pBIC')
  Y1 = Y
  best_score = Inf

  for (lambda2 in lambda2_list) {
    if (verbose) cat("Trying lambda2 =", lambda2, "\n")
    ret_pre = get_mean_var_shifts(Y, tree, alpha, Inf, lambda2, max.steps = max.steps, verbose = verbose)
    Sigma = V*ret_pre$sigma2 + (X %*% diag(ret_pre$gamma) %*% t(X))*V +
            X %*% diag(ret_pre$gamma*q) %*% t(X)
    InvSigma = solve(Sigma, tol = exp(-100))
    svd_result = svd(InvSigma)
    SqrtInvSigma = svd_result$u %*% diag(sqrt(svd_result$d)) %*% t(svd_result$v)
    YY = SqrtInvSigma %*% Y
    XX = SqrtInvSigma %*% X
    lambda1 = cv.glmnet(XX, YY, lambda = lambda1_list, intercept = FALSE, nfolds = nfolds)[[lambda.type]]
    ret = get_mean_var_shifts(Y, tree, alpha, lambda1, lambda2, max.steps = max.steps,
                              measurement_error = measurement_error, verbose = verbose)
    OModel = fit_OU_mean_var(tree, Y1, alpha, ret$sv_mean, ret$sv_var,
                              measurement_error = measurement_error,
                              max.num.shifts = max.num.shifts, verbose = verbose)
    score_summary[nrow(score_summary) + 1, ] = c(lambda1, lambda2,
                                                   paste(ret$sv_mean, collapse = ';'),
                                                   paste(ret$sv_var, collapse = ';'),
                                                   OModel$loglik, OModel$BIC, OModel$mBIC, OModel$pBIC)
    score_o = OModel[[criterion]][[1]]
    if (!is.finite(score_o)) score_o = Inf
    if (isTRUE(score_o < best_score)) {
      best_score = score_o
      best_Model = OModel
      best_Model$lambda1 = lambda1
      best_Model$lambda2 = lambda2
    }
  }

  if (!exists("best_Model")) {
    best_Model = OModel
    best_Model$lambda1 = lambda1
    best_Model$lambda2 = lambda2
  }

  score_summary = score_summary[!duplicated(score_summary[, c('sv_mean','sv_var')]), ]
  score_summary = cbind(score_summary, NA, NA, NA, NA, NA, NA)
  names(score_summary)[9:14] = c('sv_mean_corrected','sv_var_corrected',
                                  'loglik_corrected','BIC_corrected','mBIC_corrected','pBIC_corrected')

  if (verbose) cat("Backward correction (top", top_k, ")\n")
  for (i in 1:min(top_k, nrow(score_summary))) {
    ind = order(as.numeric(score_summary[[criterion]]))[i]
    sv_mean = as.numeric(strsplit(score_summary$sv_mean[ind], ';')[[1]])
    sv_var = as.numeric(strsplit(score_summary$sv_var[ind], ';')[[1]])
    OModel = backward_correction(tree, Y1, alpha, sv_mean, sv_var,
                                  measurement_error = measurement_error,
                                  max.num.shifts = max.num.shifts, verbose = verbose)
    score_summary[ind, 9:14] = c(paste(which(OModel$beta != 0), collapse = ';'),
                                  paste(which(OModel$gamma != 0), collapse = ';'),
                                  OModel$loglik, OModel$BIC, OModel$mBIC, OModel$pBIC)
    score_o = OModel[[criterion]][[1]]
    if (!is.finite(score_o)) score_o = Inf
    if (isTRUE(score_o < best_score)) {
      best_score = score_o
      best_Model = OModel
      best_Model$lambda1 = score_summary$lambda1[ind]
      best_Model$lambda2 = score_summary$lambda2[ind]
    }
  }

  if (verbose) cat("Best", criterion, "=", best_score, "\n")
  return(list(best_model = best_Model, score_summary = score_summary))
}

# ---- Main wrapper ----

#' @title ShiVa: Automatic Shift Detection in Mean and Variance
#' @description Detects evolutionary shifts in both optimal trait values and diffusion variance
#' under an Ornstein-Uhlenbeck process. Optionally refits alpha after shift detection.
#' @param Y Trait values at tips.
#' @param tree Phylogenetic tree of class \code{phylo}.
#' @param alpha Selection strength. If \code{NULL}, estimated via \code{phylolm}.
#' @param t Step size. Default 0.01.
#' @param lambda1_list Candidate lambda1 values.
#' @param lambda2_list Candidate lambda2 values. Default \code{exp(1:10 - 6)}.
#' @param criterion \code{"BIC"}, \code{"mBIC"}, or \code{"pBIC"}. Default \code{"BIC"}.
#' @param max.steps Max iterations. Default 300.
#' @param nfolds CV folds. Default 8.
#' @param top_k Top candidates for backward correction. Default 10.
#' @param lambda.type \code{"lambda.1se"} (default) or \code{"lambda.min"}.
#' @param measurement_error Logical. Default FALSE.
#' @param refit_alpha Logical. If TRUE, refit alpha after detection. Default TRUE.
#' @param verbose Logical. Default FALSE.
#' @return List with best_model and score_summary.
#' @export
#' @importFrom phylolm phylolm
ShiVa = function(Y, tree, alpha = NULL, t = 0.01,
                  lambda1_list = NULL, lambda2_list = exp(1:10 - 6),
                  criterion = "BIC", max.steps = 300,
                  nfolds = 5, top_k = 3,
                  lambda.type = "lambda.1se",
                  max.num.shifts = Inf,
                  measurement_error = FALSE,
                  refit_alpha = TRUE, verbose = FALSE) {

  X = generate_design_matrix(tree, 'simpX')

  # Step 1: estimate alpha if not provided
  if (is.null(alpha)) {
    alpha = phylolm(Y ~ 1, phy = tree, model = "OUfixedRoot")$optpar
    if (verbose) cat("Estimated alpha =", alpha, "\n")
  }

  # Step 2: run model selection
  result = get_mean_var_shifts_model_selection(Y, tree, alpha,
              lambda1_list, lambda2_list, criterion, max.steps,
              nfolds, top_k, lambda.type, max.num.shifts,
              measurement_error, verbose)
  result$alpha = alpha

  # Step 3: optionally refit alpha
  if (refit_alpha) {
    sv_mean = which(result$best_model$beta != 0)
    sv_var  = which(result$best_model$gamma != 0)

    refit_fit = if (length(sv_mean) > 0) {
      phylolm(Y ~ X[, sv_mean, drop = FALSE], phy = tree, model = "OUfixedRoot")
    } else {
      phylolm(Y ~ 1, phy = tree, model = "OUfixedRoot")
    }
    alpha_refit = refit_fit$optpar

    refit_model = fit_OU_mean_var(tree, Y, alpha_refit, sv_mean, sv_var,
                                   measurement_error = measurement_error,
                                   max.num.shifts = max.num.shifts, verbose = verbose)

    score_refit = refit_model[[criterion]][[1]]
    score_orig  = result$best_model[[criterion]][[1]]
    if (!is.finite(score_refit)) score_refit = Inf
    if (!is.finite(score_orig))  score_orig  = Inf
    if (isTRUE(score_refit < score_orig)) {
      result$best_model = refit_model
      result$alpha_refit = alpha_refit
      if (verbose) cat("Refit improved: alpha_refit =", alpha_refit, "\n")
    } else {
      result$alpha_refit = alpha
      if (verbose) cat("Refit did not improve, kept original alpha\n")
    }
  }

  return(result)
}

# ---- Utility functions ----

#' @importFrom MASS mvrnorm
get_test_data = function(tree, Sigma, n_test, alpha, beta) {
  X = generate_design_matrix(tree, type = 'simpX')
  eps = mvrnorm(n = 1, mu = rep(0, nrow(X)), Sigma = Sigma)
  ret = as.data.frame(t(X %*% beta + eps))
  for (i in 2:n_test) {
    eps = mvrnorm(n = 1, mu = rep(0, nrow(X)), Sigma = Sigma)
    Y = X %*% beta + eps
    ret[nrow(ret) + 1, ] = Y
  }
  return(as.matrix(ret))
}

get_prediction_likelihood = function(Y_pred, Sigma, test_data) {
  InvSigma = solve(Sigma)
  n = length(Y_pred)
  cmp_likelihood = function(target, prediction, Sigma, InvSigma) {
    return(-1/2*n*log(2*pi) - 1/2*t(target - prediction) %*% (InvSigma %*% (target - prediction)) -
           1/2 * determinant(Sigma)$modulu)
  }
  loglik = apply(test_data, 1, cmp_likelihood, prediction = Y_pred, Sigma = Sigma, InvSigma = InvSigma)
  return(mean(loglik))
}

# ---- S3 methods ----

# Internal helper: extract shift info from a ShiftModel object
extract_shifts <- function(model) {
  tree = model$tree
  edge_nodes <- tree$edge[, 2]
  mean_shifts <- which(model$beta != 0)
  var_shifts <- which(model$gamma != 0)
  list(
    mean = data.frame(edge = mean_shifts, node = edge_nodes[mean_shifts],
                      size = model$beta[mean_shifts]),
    var  = data.frame(edge = var_shifts, node = edge_nodes[var_shifts],
                      size = model$gamma[var_shifts]),
    sigma2 = signif(model$sigma2, 2),
    alpha  = signif(model$alpha, 2),
    loglik = round(model$loglik, 2),
    BIC    = round(model$BIC, 2)
  )
}

#' @title Plot Method for ShiftModel Objects
#' @description Plots a phylogenetic tree with detected shifts highlighted.
#' @param x A ShiftModel object.
#' @param title Plot title.
#' @param ... Additional arguments passed to plot.phylo.
#' @return No return value (side effect: plot).
#' @method plot ShiftModel
#' @importFrom grDevices heat.colors
#' @importFrom graphics legend
#' @export
plot.ShiftModel <- function(x, title = "", ...) {
  tree <- x$tree
  trait <- x$Y
  if (is.null(names(trait)) || length(names(trait)) == 0)
    names(trait) <- tree$tip.label
  shift_info <- extract_shifts(x)
  plot(tree, show.tip.label = TRUE, main = title, cex = 0.6, ...)
  tiplabels(pch = 15, adj = 0.48,
            col = rev(heat.colors(100))[as.numeric(cut(trait, 100))],
            cex = 1.5)
  if (nrow(shift_info$mean) > 0) {
    edgelabels("*", shift_info$mean$edge, col = "green", cex = 2,
               adj = c(0.5, 0.8), frame = "none")
    edgelabels(round(shift_info$mean$size, 3), shift_info$mean$edge,
               adj = c(0.6, -0.5), frame = "none", cex = 0.8, col = "darkgreen")
  }
  if (nrow(shift_info$var) > 0) {
    edgelabels("*", shift_info$var$edge, col = "red", cex = 2,
               adj = c(0.5, 0.8), frame = "none")
    edgelabels(round(shift_info$var$size, 3), shift_info$var$edge,
               adj = c(0.6, -0.5), frame = "none", cex = 0.8, col = "darkred")
  }
  legend("topright", legend = c(
    paste0("sigma2 = ", shift_info$sigma2),
    paste0("alpha = ", shift_info$alpha),
    paste0("logLik = ", shift_info$loglik),
    paste0("BIC = ", shift_info$BIC)
  ), bty = "n", cex = 0.8)
}

#' @title Summary of a ShiVa Shift Model
#' @description Summary of detected shifts and fitted parameters.
#' @param object A ShiftModel object.
#' @param ... Unused.
#' @return A summary.ShiftModel object.
#' @method summary ShiftModel
#' @export
summary.ShiftModel <- function(object, ...) {
  shift_info <- extract_shifts(object)
  structure(list(
    alpha = shift_info$alpha, sigma2 = shift_info$sigma2,
    loglik = shift_info$loglik, BIC = shift_info$BIC,
    mean_shifts = shift_info$mean, var_shifts = shift_info$var
  ), class = "summary.ShiftModel")
}

#' @title Print Method for Summary of ShiftModel
#' @param x A summary.ShiftModel object.
#' @param ... Unused.
#' @return No return value (side effect: console output).
#' @method print summary.ShiftModel
#' @export
print.summary.ShiftModel <- function(x, ...) {
  cat("ShiVa Model Summary: Shift Detection in Mean and Variance\n")
  cat("----------------------------------------------------------\n")
  cat("alpha:", x$alpha, "\n")
  cat("sigma2:", x$sigma2, "\n")
  cat("Log-likelihood:", x$loglik, "\n")
  cat("BIC:", x$BIC, "\n\n")
  if (nrow(x$mean_shifts) > 0) {
    cat("Shifts in Optimal Value (theta):\n")
    print(data.frame(Edge = x$mean_shifts$edge, Node = x$mean_shifts$node,
                     Magnitude = round(x$mean_shifts$size, 4)), row.names = FALSE)
    cat("\n")
  } else cat("No shifts in optimal value detected.\n\n")
  if (nrow(x$var_shifts) > 0) {
    cat("Shifts in Variance (sigma^2):\n")
    print(data.frame(Edge = x$var_shifts$edge, Node = x$var_shifts$node,
                     Magnitude = round(x$var_shifts$size, 4)), row.names = FALSE)
    cat("\n")
  } else cat("No shifts in variance detected.\n")
}
