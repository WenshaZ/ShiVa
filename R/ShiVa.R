#' @title generate_design_matrix
#' @description generate design matrix
#' @param tree phylogenetic tree
#' @export
#' @importFrom igraph graph.edgelist
#' @importFrom igraph get.shortest.paths
#' @return
#' \item{X}{design matrix}

generate_design_matrix = function(tree){
  nTips = length(tree$tip.label)
  nEdges = Nedge(tree)
  X = matrix(0, nTips, nEdges)
  path = nodepath(tree)
  edges = graph.edgelist(tree$edge, directed = TRUE)
  X = matrix(0, nTips, nEdges)
  path = get.shortest.paths(edges, nTips+1, to = 1:nTips, mode = "out",
                                output = "epath")$epath
  for(i in 1:nTips){
    X[i,path[[i]]] = 1
  }
  return(X)
}

#' @title soft_thresholding
#' @description Soft thresholding function
#' @param z input value
#' @param lambda degree of shrinkage
#' @return
#' \item{new_value}{The result of soft thresholding function}
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

#' @title update_step_gamma
#' @description This function calculate the updated gamma for each step
#' @param gamma_k The gamma value to be updated
#' @param X_k The kth column of the design matrix
#' @param Sigma The covariance matrix
#' @param r The residual vector (Y-X*beta)
#' @param lambda2 The degree of shrinkage for gamma
#' @param t step size
#' @param penalty "L1" or "None"
#' @param V The original covariance matrix when sigma^2 = 1
#' @param q_k The kth element of designed vector q
#' @return
#' \item{gamma_k}{The updated gamma value}
#' \item{Sigma}{The updated covariance matrix}
#' @export
#' @importFrom psych tr

update_step_gamma  = function(gamma_k, X_k, Sigma, r, lambda2, t,penalty,V,q_k){
  XXT = (X_k%*%t(X_k))
  M_k = XXT*V
  InvSigma = solve(Sigma)
  InvSigmar = InvSigma %*% r
  delta = -1/2*(t(InvSigmar)%*%M_k%*%InvSigmar) - 1/2*q_k*(t(X_k)%*%InvSigmar)^2 +1/2*tr(M_k%*%InvSigma)+1/2*q_k*(t(X_k)%*%InvSigma%*%X_k)
  #gx = 1/2* t(r)%*% InvSigmar + 1/2*log(det(Sigma))
  z = (gamma_k - t*delta)[1]

  if(penalty == 'L1'){
    new_value = soft_thresholding(z,lambda2)
  }else if(penalty == 'None'){
    new_value = z
  }

  Sigma = Sigma + (new_value - gamma_k)*M_k + (new_value-gamma_k)*q_k*XXT
  #Gt = (gamma_k - new_value)/t_k
  #gx_new =  1/2* t(r)%*% (solve(Sigma) %*%r) + 1/2*log(det(Sigma))

  #if(gx_new > gx - t_k*delta*Gt){t_k=t_k*line_search}
  return(list('gamma_k'= new_value,
              'Sigma' = Sigma
  ))
}

#' @title get_mean_var_shifts
#' @description This function provide estimation of shifts in optimal value and shifts in variance given fixed lambda1 and lambda2
#' @param Y a vector of trait values
#' @param tree a phylogenetic tree
#' @param alpha the parameter for selection strength
#' @param lambda1 The degree of shrinkage for beta (shifts in optimal value)
#' @param lambda2 The degree of shrinkage for gamma (shifts in variance)
#' @param max.steps maximum number of steps, default 1000
#' @param t step size, default 0.01
#' @param penalty "L1" or "None"
#' @param thres threshold for stopping optimization
#' @param sigma2 the original sigma2, default is unknown.
#' @return
#' \item{beta}{Shifts in optimal value}
#' \item{gamma}{Shifs in variance}
#' \item{sigma2}{Original variance without shifts}
#' \item{b0}{The intercept term}
#' @export
#' @import ape
#' @importFrom glmnet glmnet
#' @importFrom PIGShift OU.vcv

get_mean_var_shifts = function(Y, tree, alpha, lambda1, lambda2, max.steps=1000, t = 0.01, penalty = 'L1', thres = 0.001,sigma2= NULL){
  X = generate_design_matrix(tree)
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
  tb = node.depth.edgelength(tree)[tree$edge[,1]]
  q = exp(-2*alpha*max(node.depth.edgelength(tree)))/(2*alpha)*(1-exp(2*alpha*tb))

  b0 = 0
  beta = rep(0,p)
  gamma = rep(0,p)
  r = Y
  #tlist = rep(t,p)

  for(s in 1:max.steps){

    #print(paste('step:',s,sep=''))
    #print(paste('sigma2:',sigma2_0,sep=''))

    last_sigma2 = sigma2_0
    last_gamma = gamma
    last_beta = beta
    last_b0 = b0
    Sigma = V*sigma2_0 +  (X%*%diag(gamma)%*%t(X))*V+ X%*%diag(gamma*q)%*%t(X)
    InvSigma = solve(Sigma)

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
    for(k in 1:p){
      #print(k)

      result = update_step_gamma(gamma[k],X[,k],Sigma,r,lambda2,t,penalty,V,q[k])
      gamma[k] = result[['gamma_k']]
      Sigma = result[['Sigma']]
      #tlist[k] = result[['t_k']]
    }

    # update sigma0
    if(is.null(sigma2)){
      InvSigma = solve(Sigma)
      InvSigmar = InvSigma%*%r
      sigma2_delta =  -1/2*t(InvSigmar)%*%V%*%InvSigmar + 1/2*tr(V%*%InvSigma)
      sigma2_0 = max(0.01,(last_sigma2 - t*sigma2_delta)[1])
    }

    #print(paste("beta:",paste(beta,collapse=',')))
    #print(paste("gamma:",paste(gamma,collapse=',')))
    #print(paste('b0:',b0))
    max_update = max(abs(c(last_beta-beta,last_gamma-gamma)))
    #print(paste('max_update:',max_update,sep=' '))
    if(max_update<thres) {
      #print(paste(s,'steps to converge',sep=' '))
      break
    }
  }
  return(list('beta'=beta,'gamma'=gamma,'sigma2'=sigma2_0,'b0'=b0))
}

#' @title fit_OU_mean_var
#' @description This function provide estimation of shift sizes in optimal value and shifts in variance given estimated shifts
#' @param tree a phylogenetic tree
#' @param Y a vector of trait values
#' @param alpha the parameter for selection strength
#' @param sv_mean shifts in optimal value
#' @param sv_var shifts in variance
#' @param max.steps maximum number of steps, default 1000
#' @param t step size, default 0.01
#' @param thres threshold for stopping optimization
#' @return
#' \item{beta}{Shift sizes for optimal value}
#' \item{gamma}{Shift sizes for variance}
#' \item{sigma2}{Original variance without shifts}
#' \item{b0}{The intercept term}
#' \item{loglik}{log-likelihood of the fitted model}
#' \item{BIC}{BIC of the fitted model}
#' \item{fitted.values}{fitted trait values}
#' \item{Sigma}{estimated covariance matrix}
#' @export
#' @import ape
#' @importFrom PIGShift OU.vcv
#' @importFrom stats coef lm
fit_OU_mean_var = function(tree, Y ,alpha, sv_mean, sv_var,max.steps=2000,t = 0.01, thres=0.01){
  X = generate_design_matrix(tree)
  if(alpha == 0){
    V = vcv(tree)
  }else{
    V = OU.vcv(tree,alpha)
  }
  n = nrow(X)
  p = ncol(X)
  tb = node.depth.edgelength(tree)[tree$edge[,1]]
  q = exp(-2*alpha*max(node.depth.edgelength(tree)))/(2*alpha)*(1-exp(2*alpha*tb))

  sigma2_0 = 1
  b0 = 0
  beta = rep(0,p)
  gamma = rep(0,p)

  r = Y
  for(s in 1:max.steps){

    #print(paste('step:',s,sep=''))
    #print(paste('sigma2:',sigma2_0,sep=''))

    last_sigma2 = sigma2_0
    last_gamma = gamma
    last_beta = beta
    last_b0 = b0
    Sigma = V*sigma2_0 +  (X%*%diag(gamma)%*%t(X))*V+ X%*%diag(gamma*q)%*%t(X)
    InvSigma = solve(Sigma)

    # update beta
    #for(k in 1:p){
    #  result = update_step_beta(beta[k], X[,k], InvSigma, r, lambda1, t,penalty)
    #  r = r - X[,k]*(result[['beta_k']] - beta[k])
    #  beta[k] = result[['beta_k']]
    #}
    svd_result = svd(InvSigma)
    SqrtInvSigma = svd_result$u %*% diag(sqrt(svd_result$d))  %*% t(svd_result$v)
    YY = SqrtInvSigma %*% Y
    if(length(sv_mean)>0){
      XX = SqrtInvSigma %*% X[,sv_mean]
      beta[sv_mean] = coef(lm(YY~XX))[-1]
      beta[is.na(beta)] = 0
    }
    r = Y - X %*% beta - b0
    #b0 = 0
    # update b0
    tmp = t(rep(1,n)) %*% InvSigma
    b0 = as.vector((tmp %*% (r+b0))/(tmp %*% rep(1,n)))
    r = r - (b0-last_b0)

    # update gamma
    for(k in sv_var){
      #print(k)
      result = update_step_gamma(gamma[k],X[,k],Sigma,r,0,t,penalty='None',V,q[k])
      gamma[k] = result[['gamma_k']]
      Sigma = result[['Sigma']]
    }

    # update sigma0
    InvSigma = solve(Sigma)
    InvSigmar = InvSigma%*%r
    sigma2_delta =  -1/2*t(InvSigmar)%*%V%*%InvSigmar + 1/2*tr(V%*%InvSigma)
    sigma2_0 = max(0,(last_sigma2 - t*sigma2_delta)[1])


    #print(paste("beta:",paste(beta,collapse=',')))
    #print(paste("gamma:",paste(gamma,collapse=',')))
    #print(paste('b0:',b0))
    max_update = max(abs(c(last_beta-beta,last_gamma-gamma,last_sigma2-sigma2_0)))
    #print(paste('max_update:',max_update,sep=' '))
    if(max_update<thres) {
      #print(paste(s,'steps to converge',sep=' '))
      break
    }
  }

  Sigma = V*sigma2_0 +  (X%*%diag(gamma)%*%t(X))*V+ X%*%diag(gamma*q)%*%t(X)
  InvSigma = solve(Sigma)
  loglik = - 1/2*(n*log(2*pi) + t(r) %*% (InvSigma %*% r) + log(det(Sigma)))
  BIC = -2*loglik + log(n)*(2*sum(beta!=0)+2*sum(gamma!=0)+3)
  #mBIC = -2*loglik + (2*sum(beta!=0)+2*sum(gamma!=0)-1)*log(n) + sum(log(colSums(X[,sv_mean])))+ sum(log(colSums(X[,sv_var])))
  #return(list('beta'=beta,'gamma'=gamma,'sigma2'=sigma2_0,'b0'=b0,'loglik' = loglik,'BIC'=BIC,'mBIC'=mBIC))
  return(list('beta'=beta,'gamma'=gamma,'sigma2'=sigma2_0,
              'b0'=b0,'loglik' = loglik,'BIC'=BIC,
              'fitted.values'=X %*% beta+b0,
              'Sigma'=Sigma))
}

#' @title get_mean_var_shifts_model_selection
#' @description This function provide estimation of shift sizes in optimal value and shifts in variance (including model selection)
#' @param Y a vector of trait values
#' @param tree a phylogenetic tree
#' @param alpha the parameter for selection strength
#' @param lambda1_list a list of degree of shrinkage parameters for shifts in optimal value
#' @param lambda2_list a list of degree of shrinkage parameters for shifts in variance
#' @param max.steps maximum number of steps, default 1000
#' @param t step size, default 0.01
#' @param penalty "L1" or "None"
#' @param thres threshold for stopping optimization, default 0.005
#' @param sigma2 the original sigma2, default is unknown.
#' @return
#' \item{beta}{Shift sizes for optimal value}
#' \item{gamma}{Shift sizes for variance}
#' \item{sigma2}{Original variance without shifts}
#' \item{b0}{The intercept term}
#' \item{loglik}{log-likelihood of the fitted model}
#' \item{BIC}{BIC of the fitted model}
#' \item{fitted.values}{fitted trait values}
#' \item{Sigma}{estimated covariance matrix}
#' \item{lambda1}{selected lambda for shifts in optimal value}
#' \item{lambda2}{selected lambda for shifts in variance}
#' @export
#' @import ape
#' @importFrom glmnet glmnet 
#' @importFrom glmnet cv.glmnet
#' @importFrom PIGShift OU.vcv
#' @importFrom graphics par
#' @examples
#' require(ape)
#' require(MASS)
#' require(PIGShift)
#' alpha = 1
#' sigma2_0 = alpha*2
#' true_shifts_var = c(10) #shift configuration in variance
#' size_var = 10 #shift size in variance
#' true_shifts_mean = 0 #shift configuration in mean
#' size_mean = 0 #shift size in mean
#' set.seed(123)
#' # generate a tree
#' tree = rcoal(20)
#' tree$edge.length = tree$edge.length/max(node.depth.edgelength(tree))
#' # generate design matrix
#' X = generate_design_matrix(tree)
#' V = OU.vcv(tree,alpha)
#' tb = node.depth.edgelength(tree)[tree$edge[,1]]
#' q = exp(-2*alpha*max(node.depth.edgelength(tree)))/(2*alpha)*(1-exp(2*alpha*tb))
#' # simulate shifts in variance
#' gamma = rep(0,ncol(X))
#' gamma[true_shifts_var]  = size_var
#' Sigma = V*sigma2_0 +  X%*%diag(gamma)%*%t(X)*V+ X%*%diag(gamma*q)%*%t(X)
#' # simulate shifts in mean
#' beta = rep(0,ncol(X))
#' beta[true_shifts_mean] = size_mean
#' # generate simulated data
#' eps = mvrnorm(n = 1, mu = rep(0,nrow(X)), Sigma = Sigma)
#' Y = X%*%beta+eps
#' result = get_mean_var_shifts_model_selection(Y,tree,alpha,NULL,exp(1:10-6))

get_mean_var_shifts_model_selection = function(Y, tree, alpha, lambda1_list=NULL, lambda2_list, max.steps=1000, t = 0.01, penalty = 'L1', thres = 0.005,sigma2= NULL){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  X = generate_design_matrix(tree)
  if(alpha == 0){
    V = vcv(tree)
  }else{
    V = OU.vcv(tree,alpha)
  }
  n = nrow(X)
  p = ncol(X)

  tb = node.depth.edgelength(tree)[tree$edge[,1]]
  q = exp(-2*alpha*max(node.depth.edgelength(tree)))/(2*alpha)*(1-exp(2*alpha*tb))
  minBIC = Inf
  ret = NULL

  BIC_list = c()
  loglik_list = c()

  for(lambda2 in lambda2_list){
    if(is.null(sigma2)){
      sigma2_0 = 1
    }else{
      sigma2_0 = sigma2
    }
    b0 = 0
    beta = rep(0,p)
    gamma = rep(0,p)
    r = Y
    #tlist = rep(t,p)
    Sigma = V*sigma2_0 +  (X%*%diag(gamma)%*%t(X))*V+ X%*%diag(gamma*q)%*%t(X)
    InvSigma = solve(Sigma)
    svd_result = svd(InvSigma)
    SqrtInvSigma = svd_result$u %*% diag(sqrt(svd_result$d))  %*% t(svd_result$v)
    YY = SqrtInvSigma %*% Y
    XX = SqrtInvSigma %*% X
    lambda1 = cv.glmnet(XX,YY,lambda=lambda1_list,intercept=FALSE)$lambda.min


    for(s in 1:max.steps){

      #print(paste('step:',s,sep=''))
      #print(paste('sigma2:',sigma2_0,sep=''))

      last_sigma2 = sigma2_0
      last_gamma = gamma
      last_beta = beta
      last_b0 = b0
      Sigma = V*sigma2_0 +  (X%*%diag(gamma)%*%t(X))*V+ X%*%diag(gamma*q)%*%t(X)
      InvSigma = solve(Sigma)

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
      for(k in 1:p){
        #print(k)

        result = update_step_gamma(gamma[k],X[,k],Sigma,r,lambda2,t,penalty,V,q[k])
        gamma[k] = result[['gamma_k']]
        Sigma = result[['Sigma']]
        #tlist[k] = result[['t_k']]
      }

      # update sigma0
      if(is.null(sigma2)){
        InvSigma = solve(Sigma)
        InvSigmar = InvSigma%*%r
        sigma2_delta =  -1/2*t(InvSigmar)%*%V%*%InvSigmar + 1/2*tr(V%*%InvSigma)
        sigma2_0 = max(0.01,(last_sigma2 - t*sigma2_delta)[1])
      }

      #print(paste("beta:",paste(beta,collapse=',')))
      #print(paste("gamma:",paste(gamma,collapse=',')))
      #print(paste('b0:',b0))
      max_update = max(abs(c(last_beta-beta,last_gamma-gamma)))
      #print(paste('max_update:',max_update,sep=' '))
      if(max_update<thres) {
        #print(paste(s,'steps to converge',sep=' '))
        break
      }
    }

    sv_mean = (1:ncol(X))[beta!=0]
    sv_var = (1:ncol(X))[gamma!=0]
    OModel = fit_OU_mean_var(tree,Y,alpha,sv_mean,sv_var)
    loglik_list =  c(loglik_list,OModel$loglik)
    BIC_list = c(BIC_list,OModel$BIC)
    if(OModel$BIC < minBIC){
      minBIC = OModel$BIC
      ret = OModel
      ret$lambda1 = lambda1
      ret$lambda2 = lambda2
    }
  }

  par(mfrow=c(1,2))
  plot(log(lambda2_list),loglik_list,type='l',main='loglik')
  plot(log(lambda2_list),BIC_list,type='l',main='BIC')

  return(ret)
}




