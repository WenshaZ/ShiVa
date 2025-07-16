updatevcv = function(ape.tre, new.rates) {
  vv = vcv(updatetree(ape.tre, new.rates))
  return(vv)
}

updatetree = function(ape.tre, new.rates) {
  ape.tre$edge.length = ape.tre$edge.length * new.rates
  return(ape.tre)
}

#' @title OU.vcv
#' @description generate covariance matrix for OU process
#' @param tree phylogenetic tree
#' @param alpha selective force
#' @return
#' \item{V}{covariance matrix}
OU.vcv = function(tree,alpha) {
  #theta is STRENGTH OF SELECTION
  times = updatevcv(tree,1)
  V = matrix(nrow=nrow(times),ncol=ncol(times))
  for (i in 1:nrow(V)) {
    for (j in 1:ncol(V)) {
      V[i,j] = 1/(2*alpha)*exp(-2*alpha*(times[i,i]-times[i,j]))*(1-exp(-2*alpha*times[i,j]))
    }
  }
  return(V)
}




#' @title Generate Design Matrix
#' @description Constructs a design matrix for a given phylogenetic tree, used in Ornstein-Uhlenbeck modeling.
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param type A character string specifying the type of design matrix to generate. Options are \code{"simpX"} or \code{"orgX"}.
#' @param alpha The selection strength parameter (only used when \code{type = "orgX"}).
#' @return A design matrix \code{X}, where each row corresponds to a tip and each column to an edge in the tree.
#' @export
#' @importFrom igraph graph.edgelist get.shortest.paths

generate_design_matrix = function(tree, type = "simpX", alpha = 0) 
{
    stopifnot(is.ultrametric(tree))
    stopifnot(sum(1:length(tree$tip.label) %in% tree$edge[, 1]) == 
        0)
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
    }
    else if (type == "simpX") {
        for (i in 1:nTips) {
            X[i, root2tip[[i]]] = 1
        }
    }
    else stop("Undefined design matrix type")
    return(X)
}


#' @title Soft Thresholding
#' @description Applies the soft thresholding operation to a numeric input, commonly used in Lasso and sparse modeling to induce shrinkage and sparsity.
#' @param z A numeric value to be thresholded.
#' @param lambda A non-negative numeric value indicating the threshold level (degree of shrinkage).
#' @return A numeric value after applying soft thresholding:
#' \eqn{\text{sign}(z) \cdot \max(|z| - \lambda, 0)}.
#' @export

soft_thresholding = function(z, lambda){
  if(z > lambda){
    new_value = z-lambda
  }else if(z < -lambda){
    new_value = z+lambda
  }else{
    new_value = 0
  }
  return(new_value)
}

#' @title Update Step for Gamma
#' @description Performs one update step for the \code{gamma_k} parameter in an iterative optimization routine, potentially applying L1 shrinkage.
#' @param gamma_k Current value of the \code{gamma_k} parameter to be updated.
#' @param X_k The \eqn{k}th column of the design matrix \code{X}.
#' @param Sigma Current covariance matrix.
#' @param r Residual vector, typically \code{Y - X \%*\% beta}.
#' @param lambda2 Non-negative tuning parameter controlling the degree of L1 shrinkage applied to \code{gamma_k}.
#' @param t Step size for the gradient update.
#' @param penalty Penalty type; either \code{"L1"} for soft thresholding or \code{"None"} for unpenalized updates.
#' @param V Baseline covariance matrix when \eqn{\sigma^2 = 1}.
#' @param q_k The \eqn{k}th element of the design vector \code{q}, used in the update.
#' @return A list containing:
#' \item{gamma_k}{The updated value of \code{gamma_k}.}
#' \item{Sigma}{The updated covariance matrix \code{Sigma}.}
#' @export
#' @importFrom psych tr

update_step_gamma <- function(gamma_k, X_k, Sigma, r, lambda2, t, penalty, V, q_k) {
  XXT <- X_k %*% t(X_k)
  M_k <- XXT * V

  # Use Cholesky-based inversion (faster and more stable)
  chol_result <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(chol_result)) {
    warning("Cholesky failed: skipping gamma update.")
    return(list(gamma_k = gamma_k, Sigma = Sigma))
  }
  InvSigma <- chol2inv(chol_result)

  InvSigmar <- InvSigma %*% r

  # Compute derivative
  delta <- -0.5 * (t(InvSigmar) %*% M_k %*% InvSigmar) -
           0.5 * q_k * (t(X_k) %*% InvSigmar)^2 +
           0.5 * tr(M_k %*% InvSigma) +
           0.5 * q_k * (t(X_k) %*% InvSigma %*% X_k)
  delta <- as.numeric(delta)

  # Skip update if gradient is very small (save time)
  if (abs(delta) < 1e-6 || !is.finite(delta)) {
    return(list(gamma_k = gamma_k, Sigma = Sigma))
  }

  # Initialize update
  z <- gamma_k - t * delta
  if (penalty == "L1") {
    new_value <- soft_thresholding(z, lambda2)
  } else {
    new_value <- z
  }

  # Clamp gamma for stability
  new_value <- max(0, min(new_value, 10))

  # Limit to 5 backtracking steps max
  for (i in 1:5) {
    diag_shift <- (new_value - gamma_k) * (q_k + V[1, 1])
    if (all(diag(Sigma)[X_k != 0] + diag_shift > 0)) {
      break
    }
    t <- t * 0.5
    z <- gamma_k - t * delta
    if (penalty == "L1") {
      new_value <- soft_thresholding(z, lambda2)
    } else {
      new_value <- z
    }
    new_value <- max(0, min(new_value, 10))
  }

  # Final Sigma update
  Sigma <- Sigma + (new_value - gamma_k) * (M_k + q_k * XXT)

  return(list(gamma_k = new_value, Sigma = Sigma))
}

#' @title Estimate Shifts in Optimal Trait Values and Variance
#' @description Estimates shifts in both the optimal trait values (mean) and evolutionary variance along a phylogeny under an Ornstein-Uhlenbeck (OU) process, using an \eqn{\ell_1}-penalized optimization procedure. Optionally accounts for measurement error in the observed trait data.
#'
#' @param Y A numeric vector of continuous trait values for the species at the tips of the tree.
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param alpha The selection strength parameter in the OU process.
#' @param lambda1 Non-negative penalty for \eqn{\beta} (shifts in optimal trait values).
#' @param lambda2 Non-negative penalty for \eqn{\gamma} (shifts in evolutionary variance).
#' @param max.steps Maximum number of optimization steps. Default is 1000.
#' @param t Step size for the gradient-based updates. Default is 0.01.
#' @param penalty Type of penalty to apply. Options are \code{"L1"} (default) or \code{"None"}.
#' @param thres Convergence threshold for the change in loss between steps. Default is 0.01.
#' @param sigma2 Optional initial value for the base evolutionary variance. If \code{NULL}, it is initialized to 1.
#' @param measurement_error Logical. If \code{TRUE}, the method estimates additional measurement error variance.
#'
#' @return A list containing:
#' \item{shifts_mean}{Indices of branches with detected shifts in optimal trait values (\eqn{\beta \neq 0}).}
#' \item{shifts_var}{Indices of branches with detected shifts in evolutionary variance (\eqn{\gamma \neq 0}).}
#' \item{beta}{Estimated shift coefficients for optimal trait values.}
#' \item{gamma}{Estimated shift coefficients for evolutionary variance.}
#' \item{sigma2}{Estimated base variance (\eqn{\sigma^2}) of the OU process.}
#' \item{b0}{Estimated intercept (root state).}
#' \item{sigma2_error}{Estimated measurement error variance (only returned if \code{measurement_error = TRUE}).}
#'
#' @export
#' @import ape
#' @importFrom glmnet glmnet


 get_mean_var_shifts = function(Y, tree, alpha, lambda1, lambda2, max.steps=1000, t = 0.01, penalty = 'L1', thres = 0.01,sigma2= NULL,measurement_error=FALSE){
   X = generate_design_matrix(tree,type='simpX')
   internal_list = (1:ncol(X))[colSums(X)>1 & colSums(X)!=(nrow(X)-1)]

   if(alpha == 0){
    V = vcv(tree)
   }else{
    V = OU.vcv(tree,alpha)
  }
   n = nrow(X)
   p = ncol(X)
   if(is.null(sigma2)){
     sigma2_0 = 1
   }else{
     sigma2_0 = sigma2
   }
   sigma2_error = if(measurement_error) 1 else 0
   tau_0 = log(sigma2_0)
   tau_error = log(sigma2_error)
   tb = node.depth.edgelength(tree)[tree$edge[,1]]
   q = exp(-2*alpha*max(node.depth.edgelength(tree)))/(2*alpha)*(1-exp(2*alpha*tb))

   b0 = 0
   beta = rep(0,p)
   gamma = rep(0,p)
   r = Y
   #tlist = rep(t,p)
    loss = Inf
#    best_loss = Inf
#    best_sigma2 = sigma2_0
#    best_sigma2_error = sigma2_error
#    best_gamma = gamma
#    best_beta = beta
#    best_b0 = b0
    Sigma = V*sigma2_0 +  (X%*%diag(gamma)%*%t(X))*V+ X%*%diag(gamma*q)%*%t(X)
    InvSigma = solve(Sigma)
    loss_list = c()

   for(s in 1:max.steps){
     
     #print(paste('step:',s,sep=''))
     #print(paste('sigma2:',sigma2_0,sep=''))
     
     last_tau = tau_0
     last_gamma = gamma
     last_beta = beta
     last_b0 = b0
     last_loss = loss
     last_tau_error = tau_error
     
     # update beta
     #for(k in 1:p){
     #  result = update_step_beta(beta[k], X[,k], InvSigma, r, lambda1, t,penalty)
     #  r = r - X[,k]*(result[['beta_k']] - beta[k]) 
     #  beta[k] = result[['beta_k']]
     #}
     svd_result = svd(InvSigma)
     SqrtInvSigma = svd_result$u %*% diag(sqrt(svd_result$d))  %*% t(svd_result$v)
     YY = SqrtInvSigma %*% Y
     XX = SqrtInvSigma %*% X
     beta = as.vector(glmnet(XX,YY,'gaussian',lambda=lambda1,intercept = FALSE)$beta)
     r = Y - X %*% beta
     #b0 = 0
     # update b0
     tmp = t(rep(1,n)) %*% InvSigma
     b0 = as.vector((tmp %*% (r+b0))/(tmp %*% rep(1,n)))
     r = r - (b0-last_b0)
     
     # update gamma
     for(k in internal_list){
       #print(k)
       result = update_step_gamma(gamma[k],X[,k],Sigma,r,lambda2,t,penalty,V,q[k])
       #print(min(diag(Sigma)))
       gamma[k] = result[['gamma_k']]
       Sigma = result[['Sigma']]
       #tlist[k] = result[['t_k']]

     } 

     # update sigma0
     if(is.null(sigma2)){
        InvSigma = solve(Sigma,tol=exp(-100))
        InvSigmar = InvSigma%*%r
        tau_delta =  (-1/2*t(InvSigmar)%*%V%*%InvSigmar + 1/2*tr(V%*%InvSigma))*exp(tau_0)
        #st = abs(2/(t(InvSigmar) %*% V %*% InvSigma %*% V %*% r-tr(V%*%InvSigma%*%V%*%InvSigma)))
        #print(paste("step size:",st,sep=""))
        #while((last_sigma2 - st*sigma2_delta)[1]*V[1,1]+min(gamma)<0){

        #print(paste("value:",min(diag(Sigma))  - st*sigma2_delta[1],sep=""))
        tau_0 = (last_tau - t*tau_delta)[1]
        Sigma = Sigma + (exp(tau_0)-exp(last_tau)) * V
        InvSigma = solve(Sigma,tol=exp(-100))
        InvSigmar = InvSigma%*%r
     }

    # update sigma_error
    if(measurement_error){  
      tau_error_delta = (-1/2*t(InvSigmar)%*%InvSigmar + 1/2*tr(InvSigma))*exp(tau_error)
      tau_error = (last_tau_error - t*tau_error_delta)[1]
      Sigma = Sigma + (exp(tau_error) - exp(last_tau_error)) * diag(1,n)
      InvSigma = solve(Sigma)
    }

      if (any(!is.finite(Sigma))) {
      stop("Sigma has NaN or Inf - likely from unstable gamma or bad matrix ops")
      }
     #print(paste("beta:",paste(beta,collapse=',')))
     #print(paste("gamma:",paste(gamma,collapse=',')))
     #print(paste('b0:',b0))
      max_update = max(abs(c(last_beta-beta,last_gamma-gamma)))
     #print(paste('max_update:',max_update,sep=' '))
      loglik = - 1/2*(n*log(2*pi) + t(r) %*% (InvSigma %*% r) +  determinant(Sigma)$modulu)
      loss = -loglik
      if(lambda1 != Inf) loss = loss + sum(lambda1*abs(beta))
      if(lambda2 != Inf) loss = loss + sum(lambda2*abs(gamma))
      #print(paste('sigma2:',exp(tau_0),sep=''))
      #print(paste('sigma2_error:',exp(tau_error),sep=''))
      #print(paste("det(Sigma):",determinant(Sigma)$modulu,sep=""))
      #print(paste("loss:",loss,sep=""))
#      if(loss < best_loss){
#          best_loss = loss
#          best_sigma2 = sigma2_0
#          best_sigma2_error = sigma2_error
#          best_gamma = gamma
#          best_beta = beta
#          best_b0 = b0
#      }
      loss_list = c(loss_list,loss)
      if(((abs(loglik)!=Inf)&(abs(loss-last_loss)<thres))|loglik==-Inf) {
        print(paste(s,'steps to converge',sep=' '))
        break
      }
    }

    sigma2_0 = exp(tau_0)
    sigma2_error = exp(tau_error)
#    sigma2_0 = best_sigma2
#    sigma2_error = best_sigma2_error
#    gamma = best_gamma 
#    beta = best_beta
#    b0 = best_b0
   
   result = list('shifts_mean' = which(beta!=0),'shifts_var' = which(gamma!=0),
  'gamma'=gamma,'beta'=beta,'sigma2'=sigma2_0,'b0'=b0)
  if(measurement_error){
    result[['sigma2_error']] = sigma2_error
  }
  #plot(1:length(loss_list),loss_list,type="l")
  return(result)
}

#' @title Fit OU Model with Shifts in Mean and Variance
#' @description Fits an Ornstein-Uhlenbeck (OU) model with user-specified shifts in both optimal trait values (mean) and evolutionary variance along a phylogeny. The method uses numerical optimization to estimate shift magnitudes, base variance, and intercept, and can optionally incorporate measurement error in trait values.
#'
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param Y A numeric vector of continuous trait values at the tips of the tree.
#' @param alpha A non-negative numeric value specifying the strength of selection in the OU process.
#' @param shifts_mean An integer vector indicating the indices of branches where shifts in the optimal trait value occur.
#' @param shifts_var An integer vector indicating the indices of branches where shifts in evolutionary variance occur.
#' @param max.steps Maximum number of optimization steps. Default is 1000.
#' @param t Step size for the optimizer. Default is 0.01.
#' @param thres Convergence threshold for change in log-likelihood. Default is 0.01.
#' @param measurement_error Logical. If \code{TRUE}, a separate measurement error variance is estimated and added to the diagonal of the covariance matrix.
#' @param max.num.shifts Maximum allowed number of shifts (combined for mean and variance). Default is \code{Inf}.
#'
#' @return A list containing:
#' \item{shifts_mean}{Indices of branches with non-zero shifts in optimal trait value.}
#' \item{shifts_var}{Indices of branches with non-zero shifts in evolutionary variance.}
#' \item{beta}{Estimated shift magnitudes for optima values.}
#' \item{gamma}{Estimated shift magnitudes for variance.}
#' \item{sigma2}{Estimated base evolutionary variance.}
#' \item{b0}{Estimated intercept (ancestral trait value).}
#' \item{sigma2_error}{Estimated measurement error variance (only returned if \code{measurement_error = TRUE}).}
#' \item{loglik}{Log-likelihood of the fitted model.}
#' \item{BIC}{BBIC for model selection.}
#' \item{mBIC}{mBIC for accounting shift sparsity and shared support.}
#' \item{pBIC}{pBIC incorporating determinant of projected design matrix.}
#' \item{fitted.values}{Fitted trait values based on the estimated model.}
#' \item{Sigma}{Estimated trait covariance matrix under the fitted model.}
#'
#' @export
#' @import ape
#' @importFrom stats optim

fit_OU_mean_var = function(tree, Y, alpha, shifts_mean, shifts_var,
                            max.steps = 1000, t = 0.01, thres = 0.01,
                            measurement_error = FALSE, max.num.shifts = Inf) {
    X <- generate_design_matrix(tree, "simpX")
    n <- nrow(X)
    p <- ncol(X)
    
    max.num.shifts = min(max.num.shifts, p)
    if (length(shifts_mean) > max.num.shifts || length(shifts_var) > max.num.shifts) {
      return(list(loglik = -Inf, BIC = Inf, mBIC = Inf, pBIC = Inf))
    }

    if (alpha == 0) {
        V <- vcv(tree)
    } else {
        V <- OU.vcv(tree, alpha)
    }
    
    tb <- node.depth.edgelength(tree)[tree$edge[, 1]]
    q <- exp(-2 * alpha * max(node.depth.edgelength(tree))) / (2 * alpha) * (1 - exp(2 * alpha * tb))
    
    par0 <- c(rep(0, length(shifts_mean)), rep(0, length(shifts_var)), 0, log(1))
    if (measurement_error) {
        par0 <- c(par0, log(1))
    }
    
    nll_fn <- function(par) {
        m <- length(shifts_mean)
        v <- length(shifts_var)
        beta <- rep(0, p)
        gamma <- rep(0, p)
        
        if (m > 0) beta[shifts_mean] <- par[1:m]
        if (v > 0) gamma[shifts_var] <- par[(m+1):(m+v)]
        
        b0 <- par[m + v + 1]
        log_sigma2 <- par[m + v + 2]
        sigma2 <- exp(log_sigma2)
        sigma2_error <- if (measurement_error) exp(par[m + v + 3]) else 0
        
        mu <- X %*% beta + b0
        r <- Y - mu
        
        Sigma <- V * sigma2 +
            (X %*% diag(gamma) %*% t(X)) * V +
            X %*% diag(gamma * q) %*% t(X) +
            sigma2_error * diag(n)
        
        cholSigma <- tryCatch(chol(Sigma), error = function(e) return(NULL))
        if (is.null(cholSigma)) return(1e10)
        
        logdet <- 2 * sum(log(diag(cholSigma)))
        InvSigma_r <- backsolve(cholSigma, forwardsolve(t(cholSigma), r))
        0.5 * (n * log(2 * pi) + logdet + sum(r * InvSigma_r))
    }
    
    opt_result <- optim(par = par0, fn = nll_fn, method = "L-BFGS-B", control = list(maxit = max.steps))
    
    opt_par <- opt_result$par
    m <- length(shifts_mean)
    v <- length(shifts_var)
    
    beta <- rep(0, p)
    gamma <- rep(0, p)
    if (m > 0) beta[shifts_mean] <- opt_par[1:m]
    if (v > 0) gamma[shifts_var] <- opt_par[(m+1):(m+v)]
    
    b0 <- opt_par[m + v + 1]
    sigma2 <- exp(opt_par[m + v + 2])
    sigma2_error <- if (measurement_error) exp(opt_par[m + v + 3]) else 0
    
    mu <- X %*% beta + b0
    r <- Y - mu
    Sigma <- V * sigma2 +
        (X %*% diag(gamma) %*% t(X)) * V +
        X %*% diag(gamma * q) %*% t(X) +
        sigma2_error * diag(n)
    InvSigma <- solve(Sigma, tol = exp(-100))
    loglik <- -nll_fn(opt_par)
    
    # Criteria
    k_beta <- sum(beta != 0)
    k_gamma <- sum(gamma != 0)
    BIC <- -2 * loglik + log(n) * (2 * k_beta + 2 * k_gamma + 3)
    
    logn_mean <- if (length(shifts_mean) > 1) sum(log(colSums(X[, shifts_mean, drop = FALSE])))
    else if (length(shifts_mean) == 1) log(sum(X[, shifts_mean]))
    else 0
    
    logn_var <- if (length(shifts_var) > 1) sum(log(colSums(X[, shifts_var, drop = FALSE])))
    else if (length(shifts_var) == 1) log(sum(X[, shifts_var]))
    else 0
    
    active <- c(shifts_mean, shifts_var)
    logn0 <- if (length(active) == 0) log(n)
    else if (length(active) == 1) log(n - sum(X[, active] > 0))
    else {
        z <- rowSums(X[, active, drop = FALSE]) > 0
        if (sum(z) < n) log(n - sum(z)) else 0
    }
    
    mBIC <- -2 * loglik + (2 * k_beta + 2 * k_gamma - 1) * log(n) + logn_mean + logn_var + logn0
    
    design_matrix <- cbind(1, X[, shifts_mean, drop = FALSE])
    Xt_Si_X <- t(design_matrix) %*% InvSigma %*% design_matrix
    logdet_proj <- as.numeric(determinant(Xt_Si_X, logarithm = TRUE)$modulus)
    
    pBIC <- -2 * loglik + 2 * (k_beta + k_gamma) * log(2 * n - 3) + 2 * log(n) + logdet_proj
    
    result <- list(
        tree = tree,
        shifts_mean = which(beta!=0),
        shifts_var = which(gamma!=0),
        alpha = alpha,
        beta = beta,
        gamma = gamma,
        sigma2 = sigma2,
        b0 = b0,
        loglik = loglik,
        BIC = BIC,
        mBIC = mBIC,
        pBIC = pBIC,
        fitted.values = mu,
        Sigma = Sigma
    )
    
    if (measurement_error) {
        result$sigma2_error <- sigma2_error
    }
    
    return(result)
}

#' @title Backward Selection for OU Model Shift Correction
#' @description Performs backward stepwise selection on a given set of candidate shifts in optimal trait values (mean) and evolutionary variance under an Ornstein-Uhlenbeck (OU) model. This function iteratively removes individual shifts to improve model fit based on a specified selection criterion.
#'
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param Y A numeric vector of trait values corresponding to the tips of the tree.
#' @param alpha A non-negative numeric value specifying the strength of selection in the OU process.
#' @param shifts_mean A vector of branch indices with candidate shifts in optimal trait values.
#' @param shifts_var A vector of branch indices with candidate shifts in evolutionary variance.
#' @param criterion A model selection criterion to guide backward elimination. Options include \code{"BIC"}, \code{"mBIC"}, or \code{"pBIC"}. Default is \code{"BIC"}.
#' @param original_model (Optional) A previously fitted OU model returned by \code{fit_OU_mean_var}. If \code{NULL}, the model is refit using the provided shifts.
#' @param measurement_error Logical. If \code{TRUE}, the model accounts for measurement error by estimating an additional variance term. Default is \code{FALSE}.
#' @param max.num.shifts An integer specifying the maximum number of total shifts (mean and variance combined) allowed in the model. Default is \code{Inf}.
#'
#' @return A fitted OU model object (a list), as returned by \code{fit_OU_mean_var}, with a potentially reduced set of shifts that minimizes the specified criterion.
#'
#' @export


backward_correction = function(tree,Y,alpha,shifts_mean,shifts_var,criterion='BIC',original_model=NULL,measurement_error = FALSE,
                               max.num.shifts = Inf){

  if(is.null(original_model)){
    OModel = fit_OU_mean_var(tree,Y,alpha,shifts_mean,shifts_var,measurement_error = measurement_error,
                              max.num.shifts = max.num.shifts)
  }else{
    original_model = OModel
  }
  original_score = OModel[[criterion]][[1]]
  best_score = original_score
  n1 = length(shifts_mean)
  n2 = length(shifts_var)
  if(n1+n2==0){
    return(OModel)
  }
  after_score_list = rep(0,n1+n2)
  for(i in 1:(n1+n2)){
    if(i <= length(shifts_mean)){
      new_Model = fit_OU_mean_var(tree,Y,alpha,setdiff(shifts_mean,c(shifts_mean[i])),shifts_var,measurement_error = measurement_error,
                              max.num.shifts = max.num.shifts)
      after_score_list[i] = new_Model[[criterion]][[1]]
    } else{
      new_Model = fit_OU_mean_var(tree,Y,alpha,shifts_mean,setdiff(shifts_var,c(shifts_var[i-n1])),measurement_error = measurement_error,
                              max.num.shifts = max.num.shifts)
      after_score_list[i] = new_Model[[criterion]][[1]]      
    }
    #print(paste("index:",i,"score:",new_Model[[criterion]][[1]],sep=' '))
    if(after_score_list[i]<best_score){
      OModel = new_Model
      best_score = after_score_list[i]
    }
  }
  index_list = order(after_score_list)[sort(after_score_list)<original_score]
  #print(paste('index_list:',paste(index_list,collapse = ';'),sep=''))
  remove_list = c(index_list[1])
  for(i in index_list[2:length(index_list)]){
    new_Model = fit_OU_mean_var(tree,Y,alpha,
                                setdiff(shifts_mean,shifts_mean[c(i,remove_list)[c(i,remove_list)<=n1]]),
                                setdiff(shifts_var,shifts_var[c(i,remove_list)[c(i,remove_list)>n1]-n1]),
                                measurement_error = measurement_error,
                                max.num.shifts = max.num.shifts)
   if(new_Model[[criterion]][[1]]<best_score){
      OModel = new_Model
      best_score = new_Model[[criterion]][[1]]
      remove_list = c(remove_list,i)
    }
  }
  return(OModel)
  }

#' @title Model Selection for OU Shifts in Optimal value and Variance
#' @description Performs model selection to estimate the locations and magnitudes of evolutionary shifts in optimal trait values (mean) and diffusion variance under an Ornstein-Uhlenbeck (OU) process. This function searches across user-defined grids of shrinkage parameters for both types of shifts, uses cross-validation for selecting \code{lambda1}, and applies backward correction to refine top candidate models.
#'
#' @param Y A numeric vector of trait values for the species at the tips of the phylogenetic tree.
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param alpha A non-negative numeric value representing the selection strength in the OU process.
#' @param lambda1_list A numeric vector of candidate \eqn{\lambda_1} values controlling shrinkage for shifts in optimal values.
#' @param lambda2_list A numeric vector of candidate \eqn{\lambda_2} values controlling shrinkage for shifts in variance. Default is exp(1:10*0.4-6)
#' @param criterion Model selection criterion to optimize. Options include \code{"BIC"}, \code{"mBIC"}, or \code{"pBIC"}. Default is \code{"BIC"}.
#' @param max.steps Maximum number of optimization steps. Default is 300.
#' @param nfolds Number of cross-validation folds for tuning \code{lambda1}. Default is 8.
#' @param top_k Number of top candidate models (ranked by criterion) to further refine using backward correction. Default is 10.
#' @param t Step size for iterative optimization. Default is 0.01.
#' @param measurement_error Logical. If \code{TRUE}, estimates a separate measurement error variance component. Default is \code{FALSE}.
#' @param lambda.type A character string specifying the cross-validation rule used to select \code{lambda1} from \code{lambda1_list}. Options are \code{"lambda.min"} (minimum CV error) and \code{"lambda.1se"} (1-SE rule, higher penalty). Default is \code{"lambda.1se"}.
#' @param max.num.shifts An integer specifying the maximum number of allowed shifts (combined across mean and variance). Default is \code{Inf}.
#'
#' @return A list containing:
#' \item{best_model}{The final selected OU model object, with estimated shifts and parameters.}
#' \item{score_summary}{A data frame summarizing the model selection results, including pre- and post-correction scores and shift locations.}
#'
#' @export
#' @import ape
#' @importFrom glmnet glmnet cv.glmnet

get_mean_var_shifts_model_selection <- function(Y, tree, alpha, t = 0.01,
                                                lambda1_list = NULL, lambda2_list = exp(1:10*0.4-6),
                                                criterion = "BIC", max.steps = 300,
                                                nfolds = 8, top_k = 10,
                                                measurement_error = FALSE,
                                                lambda.type = "lambda.1se",
                                                max.num.shifts = Inf
                                                ) {

  top_k = top_k
  X <- generate_design_matrix(tree, 'simpX')
  Y <- as.vector(Y)
  V <- if (alpha == 0) vcv(tree) else OU.vcv(tree, alpha)

  n <- nrow(X)
  p <- ncol(X)
  tb <- node.depth.edgelength(tree)[tree$edge[,1]]
  q <- exp(-2 * alpha * max(node.depth.edgelength(tree))) / (2 * alpha) * (1 - exp(2 * alpha * tb))

  max.num.shifts = min(max.num.shifts, p)
  score_summary <- data.frame(matrix(nrow = 0, ncol = 8))
  names(score_summary) <- c('lambda1', 'lambda2', 'shifts_mean', 'shifts_var',
                            'loglik', 'BIC', 'mBIC', 'pBIC')

  Y1 <- Y
  best_score <- Inf
  model_counter <- 1

  for (lambda2 in lambda2_list) {
    cat("====== Model Selection Round", model_counter, "======\n")
    cat("Trying lambda2 =", lambda2, "...\n")

    ret_pre <- get_mean_var_shifts(Y, tree, alpha, Inf, lambda2, max.steps = max.steps)
    Sigma <- V * ret_pre$sigma2 +
             (X %*% diag(ret_pre$gamma) %*% t(X)) * V +
             X %*% diag(ret_pre$gamma * q) %*% t(X)
    InvSigma <- solve(Sigma, tol = exp(-100))
    svd_result <- svd(InvSigma)
    SqrtInvSigma <- svd_result$u %*% diag(sqrt(svd_result$d)) %*% t(svd_result$v)

    YY <- SqrtInvSigma %*% Y
    XX <- SqrtInvSigma %*% X
    lambda1 <- cv.glmnet(XX, YY, lambda = lambda1_list, intercept = FALSE, nfolds = nfolds)[[lambda.type]]
    cat("Selected lambda1 from CV:", lambda1, "\n")

    ret <- get_mean_var_shifts(Y, tree, alpha, lambda1, lambda2,
                               t = t, max.steps = max.steps, measurement_error = measurement_error)

    cat("  shifts_mean =", if (length(ret$shifts_mean)) paste(ret$shifts_mean, collapse = ",") else "none", "\n")
    cat("  shifts_var  =", if (length(ret$shifts_var)) paste(ret$shifts_var, collapse = ",") else "none", "\n")

    OModel <- fit_OU_mean_var(tree, Y1, alpha, ret$shifts_mean, ret$shifts_var,
                              t = t, measurement_error = measurement_error,
                                                max.num.shifts = max.num.shifts)
    cat("  log-likelihood =", OModel$loglik, "\n\n")

    score_summary[nrow(score_summary) + 1, ] <- c(lambda1, lambda2,
                                                  paste(ret$shifts_mean, collapse = ';'),
                                                  paste(ret$shifts_var, collapse = ';'),
                                                  OModel$loglik,
                                                  OModel$BIC,
                                                  OModel$mBIC,
                                                  OModel$pBIC)
    if (OModel[[criterion]][[1]] < best_score) {
      best_score <- OModel[[criterion]][[1]]
      best_Model <- OModel
      best_Model$lambda1 <- lambda1
      best_Model$lambda2 <- lambda2
    }

    model_counter <- model_counter + 1
  }

  score_summary <- score_summary[!duplicated(score_summary[, c('shifts_mean', 'shifts_var')]), ]
  score_summary <- cbind(score_summary, NA, NA, NA, NA, NA, NA)
  names(score_summary)[9:14] <- c('shifts_mean_corrected', 'shifts_var_corrected',
                                  'loglik_corrected', 'BIC_corrected',
                                  'mBIC_corrected', 'pBIC_corrected')

  cat("\n====== Backward Correction (Top", top_k, ") ======\n")
  for (i in 1:min(top_k, nrow(score_summary))) {
    ind <- order(as.numeric(score_summary[[criterion]]))[i]
    shifts_mean <- as.numeric(strsplit(score_summary$shifts_mean[ind], ';')[[1]])
    shifts_var <- as.numeric(strsplit(score_summary$shifts_var[ind], ';')[[1]])

    cat("Correcting model", i, "with shifts_mean =", paste(shifts_mean, collapse = ","),
        "shifts_var =", paste(shifts_var, collapse = ","), "...\n")

    OModel <- backward_correction(tree, Y1, alpha, shifts_mean, shifts_var, measurement_error = measurement_error,
                                                max.num.shifts = max.num.shifts)

    score_summary[ind, 9:14] <- c(paste(which(OModel$beta != 0), collapse = ';'),
                                  paste(which(OModel$gamma != 0), collapse = ';'),
                                  OModel$loglik,
                                  OModel$BIC,
                                  OModel$mBIC,
                                  OModel$pBIC)

    if (OModel[[criterion]][[1]] < best_score) {
      best_score <- OModel[[criterion]][[1]]
      best_Model <- OModel
      best_Model$lambda1 <- score_summary$lambda1[ind]
      best_Model$lambda2 <- score_summary$lambda2[ind]
    }
  }

  cat("====== Selection Finished. Best", criterion, "=", best_score, "======\n")
  return(list('best_model' = best_Model, 'score_summary' = score_summary))
}

#' @title ShiVa: Automatic Shift Detection in Mean and Variance
#' @description Performs automatic detection of evolutionary shifts in both optimal trait values (mean) and diffusion variance under an Ornstein-Uhlenbeck (OU) process. This function serves as a wrapper for \code{get_mean_var_shifts_model_selection}, with the added ability to automatically estimate the selection strength parameter \code{alpha} if not provided.
#'
#' @param Y A numeric vector of trait values at the tips of the phylogenetic tree.
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param alpha (Optional) A non-negative numeric value specifying the OU selection strength. If \code{NULL}, it is estimated via maximum likelihood using \code{phylolm()}.
#' @param t Step size for optimization. Default is 0.01.
#' @param lambda1_list A numeric vector of candidate \eqn{\lambda_1} values for penalizing mean shifts.
#' @param lambda2_list A numeric vector of candidate \eqn{\lambda_2} values for penalizing variance shifts. Default is \code{exp(1:10 * 0.4 - 6)}.
#' @param criterion Model selection criterion to use. Options are \code{"BIC"}, \code{"mBIC"}, or \code{"pBIC"}. Default is \code{"BIC"}.
#' @param max.steps Maximum number of optimization steps. Default is 300.
#' @param nfolds Number of folds for cross-validation in tuning \code{lambda1}. Default is 8.
#' @param top_k Number of top candidate models (based on criterion) to refine using backward correction. Default is 10.
#' @param measurement_error Logical. If \code{TRUE}, estimates a measurement error variance term. Default is \code{FALSE}.
#' @param lambda.type Cross-validation rule for selecting \code{lambda1}. Options are \code{"lambda.min"} or \code{"lambda.1se"}. Default is \code{"lambda.1se"}.
#' @param max.num.shifts Maximum number of allowed shifts (in both mean and variance). Default is \code{Inf}.
#'
#' @return A list with the same structure as \code{get_mean_var_shifts_model_selection}:
#' \item{best_model}{The final selected OU model object.}
#' \item{score_summary}{A data frame summarizing candidate models and backward-corrected scores.}
#'
#' @export
#' @importFrom phylolm phylolm

ShiVa = function(Y, tree, alpha=NULL, t = 0.01,
                  lambda1_list = NULL, lambda2_list = exp(1:10*0.4-6),
                  criterion = "BIC", max.steps = 300,
                  nfolds = 8, top_k = 10,
                  measurement_error = FALSE,
                  lambda.type = "lambda.1se",
                  max.num.shifts = Inf
                  ){
  if(is.null(alpha)){
   alpha = phylolm(Y~1, phy= tree,model="OUfixedRoot")$optpar
  }
  result = get_mean_var_shifts_model_selection(Y, tree, alpha, t,
                  lambda1_list, lambda2_list,
                  criterion, max.steps,
                  nfolds, top_k ,
                  measurement_error,
                  lambda.type,
                  max.num.shifts)
  return(result)
}

#' @importFrom MASS mvrnorm
get_test_data = function(tree, Sigma, n_test, alpha,beta){
  X = generate_design_matrix(tree,type='simpX')

  eps = mvrnorm(n = 1, mu = rep(0,nrow(X)), Sigma = Sigma)
  ret = as.data.frame(t(X%*%beta+eps))

  for(i in 2:n_test){
      eps = mvrnorm(n = 1, mu = rep(0,nrow(X)), Sigma = Sigma)
      Y = X %*% beta + eps
      ret[nrow(ret)+1,] = Y
  }
  return(as.matrix(ret))

}

get_prediction_likelihood = function(Y_pred,Sigma,test_data){

    #Y_pred = OModel$fitted.values
    #Sigma = OModel$Sigma
    InvSigma = solve(Sigma)
    n = length(Y_pred)

    cmp_likelihood = function(target,prediction,Sigma,InvSigma){
        return(-1/2*n*log(2*pi)-1/2*t(target-prediction)%*%(InvSigma %*% (target-prediction)) - 1/2* determinant(Sigma)$modulu)
            }

    loglik = apply(test_data,1,cmp_likelihood,prediction=Y_pred,Sigma = Sigma,InvSigma = InvSigma)

    return(mean(loglik))

}


# Extract shift info from ShiVa-style model
extract_shifts <- function(model, tree) {
  edge_nodes <- tree$edge[, 2]
  mean_shifts <- which(model$beta != 0)
  var_shifts <- which(model$gamma != 0)
  
  list(
    mean = data.frame(
      edge = mean_shifts,
      node = edge_nodes[mean_shifts],
      size = model$beta[mean_shifts]
    ),
    var = data.frame(
      edge = var_shifts,
      node = edge_nodes[var_shifts],
      size = model$gamma[var_shifts]
    ),
    sigma2 = signif(model$sigma2, 2),
    alpha = signif(model$alpha, 2),
    loglik = round(model$loglik, 2),
    BIC = round(model$BIC, 2)
  )
}



#' @title plot_model_shifts
#' @description Plots a phylogenetic tree with trait values and highlights branches with shifts in mean (green circles) and variance (pink squares).
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param trait A named numeric vector of trait values matching the tree tips.
#' @param model A fitted ShiVa model.
#' @param title Title for the plot.
#' @return A plot showing trait-colored tips and shift annotations.
#' @import ape
#' @importFrom grDevices heat.colors
#' @importFrom graphics legend
#' @export
plot_model_shifts <- function(tree, trait, model, title = "") {
  if (is.null(names(trait)) || length(names(trait)) == 0) {
    names(trait) <- tree$tip.label
  }

  shift_info <- extract_shifts(model, tree)


  plot(tree, show.tip.label = TRUE, main = title ,cex = 0.6)

  # Add colored boxes at tip positions
  tiplabels(pch = 15, adj = 0.48,
            col = rev(heat.colors(100))[as.numeric(cut(trait, 100))], 
            cex = 1.5)

  edgelabels("*",shift_info$mean$edge,col = "green", cex = 2,adj = c(0.5,0.8), bg=NULL,frame="none")
  edgelabels(round(shift_info$mean$size,3),shift_info$mean$edge,adj = c(0.6,-0.5),bg=NULL,frame="none",cex = 0.8, col = "darkgreen")
  edgelabels("*",shift_info$var$edge,col = "red",cex = 2,adj = c(0.5,0.8), bg=NULL,frame="none")
  edgelabels(round(shift_info$var$size,3),shift_info$var$edge,adj = c(0.6,-0.5),bg=NULL,frame="none", cex = 0.8, col = "darkred")



  # Add model summary
  legend("topright", legend = c(
    paste0("sigma2 = ", shift_info$sigma2),
    paste0("alpha = ", shift_info$alpha),
    paste0("logLik = ", shift_info$loglik),
    paste0("BIC = ", shift_info$BIC)
  ), bty = "n", cex = 0.8)
}


#' @title Summary of Shifts from a ShiVa Model Result
#' @description Print a concise summary of the best model selected by \code{ShiVa()}, including shift locations (edges) and magnitudes, as well as estimated parameters.
#'
#' @param result A list returned by the \code{ShiVa()} function, containing at least \code{best_model}.
#'
#' @export
summary_shifts <- function(result) {
  model <- result$best_model
  tree <- result$best_model$tree
  shift_info <- extract_shifts(model, tree)

  cat("Model Summary (ShiVa)\n")
  cat("---------------------\n")
  cat("alpha:", shift_info$alpha, "\n")
  cat("sigma2: ", shift_info$sigma2, "\n")
  cat("Log-likelihood:", shift_info$loglik, "\n")
  cat("BIC:", shift_info$BIC, "\n\n")

  if (nrow(shift_info$mean) > 0) {
    cat("Shifts in Optimal Value:\n")
    print(data.frame(
      Edge = shift_info$mean$edge,
      Node = shift_info$mean$node,
      Magnitude = round(shift_info$mean$size, 4)
    ), row.names = FALSE)
    cat("\n")
  } else {
    cat("No shifts in optimal value detected.\n\n")
  }

  if (nrow(shift_info$var) > 0) {
    cat("Shifts in Variance:\n")
    print(data.frame(
      Edge = shift_info$var$edge,
      Node = shift_info$var$node,
      Magnitude = round(shift_info$var$size, 4)
    ), row.names = FALSE)
    cat("\n")
  } else {
    cat("No shifts in variance detected.\n")
  }
}
