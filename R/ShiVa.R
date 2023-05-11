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
  edges = igraph::graph.edgelist(tree$edge, directed = TRUE)
  X = matrix(0, nTips, nEdges)
  path = igraph::get.shortest.paths(edges, nTips+1, to = 1:nTips, mode = "out",
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

update_step_gamma = function (gamma_k, X_k, Sigma, r, lambda2, t, penalty, V, q_k) {
    XXT = (X_k %*% t(X_k))
    M_k = XXT * V
    InvSigma = solve(Sigma,tol=exp(-100))
    InvSigmar = InvSigma %*% r
    delta = -1/2 * (t(InvSigmar) %*% M_k %*% InvSigmar) - 1/2 * 
        q_k * (t(X_k) %*% InvSigmar)^2 + 1/2 * tr(M_k %*% InvSigma) + 
        1/2 * q_k * (t(X_k) %*% InvSigma %*% X_k)
    new_value = -Inf
    #while(new_value+sigma2_0*V[1,1]<0){
    while(min(diag(Sigma)[X_k!=0])+(new_value-gamma_k)*(q_k+V[1,1])<0){
      #print(min(diag(Sigma)[X_k!=0]))
      #print(paste("value:",min(diag(Sigma)[X_k!=0])+(new_value-gamma_k)*(q_k+V[1,1]),sep=""))
      #print(t)
      #print(paste("before:",gamma_k," after:",new_value," t:",t,sep=""))
      z = (gamma_k - t * delta)[1]
      if (penalty == "L1") {
        new_value = soft_thresholding(z, lambda2)
      }
      else if (penalty == "None") {
        new_value = z
      }
      t = t*0.5
      if(t<exp(-10)) {
        new_value = gamma_k
        break
      }
    }
    #print(paste("value:",min(diag(Sigma)[X_k!=0])+(new_value-gamma_k)*(q_k+V[1,1]),sep=""))
    #print(paste("out:",new_value,paste=""))
    Sigma = Sigma + (new_value - gamma_k) * M_k + (new_value - 
        gamma_k) * q_k * XXT
    return(list(gamma_k = new_value, Sigma = Sigma))
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
#' \item{sv_mean}{Shifts in optimal value}
#' \item{sv_var}{Shifs in variance}
#' \item{beta}{Shift coefficients in optimal value}
#' \item{gamma}{Shift coefficients in variance}
#' \item{sigma2}{Original variance without shifts}
#' \item{b0}{The intercept term}
#' @export
#' @import ape
#' @importFrom glmnet glmnet
#' @importFrom PIGShift OU.vcv

 get_mean_var_shifts = function(Y, tree, alpha, lambda1, lambda2, max.steps=1000, t = 0.01, penalty = 'L1', thres = 0.01,sigma2= NULL){
   X = generate_design_matrix(tree)
   internal_list = (1:ncol(X))[colSums(X)>1]
   
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
    loss = Inf
    best_loss = Inf
    best_sigma2 = sigma2_0
    best_gamma = gamma
    best_beta = beta
    best_b0 = b0
    Sigma = V*sigma2_0 +  (X%*%diag(gamma)%*%t(X))*V+ X%*%diag(gamma*q)%*%t(X)
    InvSigma = solve(Sigma)

   for(s in 1:max.steps){
     
     #print(paste('step:',s,sep=''))
     #print(paste('sigma2:',sigma2_0,sep=''))
     
     last_sigma2 = sigma2_0
     last_gamma = gamma
     last_beta = beta
     last_b0 = b0
     last_loss = loss

     
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
        sigma2_delta =  -1/2*t(InvSigmar)%*%V%*%InvSigmar + 1/2*tr(V%*%InvSigma)
        st = t
        #st = abs(2/(t(InvSigmar) %*% V %*% InvSigma %*% V %*% r-tr(V%*%InvSigma%*%V%*%InvSigma)))
        #print(paste("step size:",st,sep=""))
        #while((last_sigma2 - st*sigma2_delta)[1]*V[1,1]+min(gamma)<0){
        while(min(diag(Sigma))  - st*sigma2_delta[1]<0){
          #print(st)
          #print(paste("value:",min(diag(Sigma))  - st*sigma2_delta[1],sep=""))
          st = st*0.5
        }
        #print(paste("value:",min(diag(Sigma))  - st*sigma2_delta[1],sep=""))
        sigma2_0 = (last_sigma2 - st*sigma2_delta)[1]
        Sigma = Sigma + (sigma2_0-last_sigma2) * V
        InvSigma = solve(Sigma,tol=exp(-100))
     }

    
     #print(paste("beta:",paste(beta,collapse=',')))
     #print(paste("gamma:",paste(gamma,collapse=',')))
     #print(paste('b0:',b0))
      max_update = max(abs(c(last_beta-beta,last_gamma-gamma)))
     #print(paste('max_update:',max_update,sep=' '))
      loglik = if(det(Sigma)>0) - 1/2*(n*log(2*pi) + t(r) %*% (InvSigma %*% r) + log(det(Sigma))) else -Inf
      loss = -loglik
      if(lambda1 != Inf) loss = loss + sum(lambda1*abs(beta))
      if(lambda2 != Inf) loss = loss + sum(lambda2*abs(gamma))
      #print(paste('sigma2:',sigma2_0,sep=''))
      #print(paste("loss:",loss,sep=""))
      if(loss < best_loss){
          best_loss = loss
          best_sigma2 = sigma2_0
          best_gamma = gamma
          best_beta = beta
          best_b0 = b0
      }

      if(((abs(loglik)!=Inf)&(abs(max_update)<thres))|sigma2_0<exp(-10)|loglik==-Inf) {
        #print(paste(s,'steps to converge',sep=' '))
        break
      }
    }

    sigma2_0 = best_sigma2
    gamma = best_gamma 
    beta = best_beta
    b0 = best_b0
   
 return(list('sv_mean' = (1:ncol(X))[beta!=0],'sv_var' = (1:ncol(X))[gamma!=0],
  'gamma'=gamma,'beta'=beta,'sigma2'=sigma2_0,'b0'=b0))

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
fit_OU_mean_var = function(tree, Y ,alpha, sv_mean, sv_var,max.steps=1000,t = 0.01, thres=0.01){
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
  loglik = Inf
  sigma2_0 = 1
  sigma2_error = 0
  b0 = 0
  beta = rep(0,p)
  gamma = rep(0,p)
  r = Y
  Sigma = V*sigma2_0 +  (X%*%diag(gamma)%*%t(X))*V+ X%*%diag(gamma*q)%*%t(X) + sigma2_error * diag(1,n)
  InvSigma = solve(Sigma)

  best_loglik = -Inf
  best_sigma2 = sigma2_0
  best_sigma2_error = sigma2_error
  best_gamma = gamma
  best_beta = beta
  best_b0 = b0
  best_Sigma = Sigma

  for(s in 1:max.steps){

    #print(paste('step:',s,sep=''))
    last_loglik = loglik
    last_sigma2 = sigma2_0
    last_sigma2_error = sigma2_error
    last_gamma = gamma
    last_beta = beta
    last_b0 = b0

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
      #print(min(diag(Sigma)))
    }

    # update sigma0
    InvSigma = solve(Sigma,tol=exp(-100))
    InvSigmar = InvSigma%*%r
    sigma2_delta =  -1/2*t(InvSigmar)%*%V%*%InvSigmar + 1/2*tr(V%*%InvSigma)
    st = abs(2/(t(InvSigmar) %*% V %*% InvSigma %*% V %*% r-tr(V%*%InvSigma%*%V%*%InvSigma)))
    #while((last_sigma2 - st*sigma2_delta)[1]*V[1,1]+min(gamma)<0){
    while(min(diag(Sigma))  - st*sigma2_delta[1]<0){
      st = st*0.5
    }
    sigma2_0 = (last_sigma2 - st*sigma2_delta)[1]
    Sigma = Sigma + (sigma2_0-last_sigma2) * V
    InvSigma = solve(Sigma,tol=exp(-100))
    
    #print(paste('sigma2:',sigma2_0,sep=''))
    #print(paste("beta:",paste(beta,collapse=',')))
    #print(paste("gamma:",paste(round(gamma[gamma!=0],3),collapse=',')))
    #print(paste('b0:',b0))
    #detSigma = if(det(Sigma)<0) 0 else det(Sigma)
    loglik = if(det(Sigma)>0) - 1/2*(n*log(2*pi) + t(r) %*% (InvSigma %*% r) + log(det(Sigma))) else -Inf
    #print(paste("loglik:", loglik,sep=""))
    #print(paste('max_update:',max_update,sep=' '))
    if(loglik > best_loglik){
        best_loglik = loglik
        best_sigma2 = sigma2_0
        best_sigma2_error = sigma2_error
        best_gamma = gamma
        best_beta = beta
        best_b0 = b0
    }
    
    if(((abs(loglik)!=Inf)&(abs(loglik-last_loglik)<thres))|sigma2_0<exp(-10)|loglik==-Inf) {
      #print(sigma2_0)
      #print(paste(s,'steps to converge',sep=' '))
      break
    }
  }

  loglik  =  best_loglik 
  sigma2_0 = best_sigma2
  sigma2_error = best_sigma2_error 
  gamma = best_gamma 
  beta = best_beta
  b0 = best_b0
  #Sigma = best_Sigma

  Sigma = V*sigma2_0 +  (X%*%diag(gamma)%*%t(X))*V+ X%*%diag(gamma*q)%*%t(X) + + sigma2_error * diag(1,n)
  InvSigma = solve(Sigma,tol=exp(-100))
  BIC = -2*loglik + log(n)*(2*sum(beta!=0)+2*sum(gamma!=0)+3)
  logn_mean = if(length(sv_mean)>1) sum(log(colSums(X[,sv_mean]))) else if(length(sv_mean)==1) log(sum(X[,sv_mean])) else 0
  logn_var = if(length(sv_var)>1) sum(log(colSums(X[,sv_var]))) else if(length(sv_var)==1) log(sum(X[,sv_var])) else 0
  logn0 = if(length(c(sv_mean,sv_var))==0) log(n) else if(length(c(sv_mean,sv_var))==1) log(n-sum(X[,c(sv_mean,sv_var)]>0)) else if (n-sum(rowSums(X[,c(sv_mean,sv_var)])>0)>0) log(n-sum(rowSums(X[,c(sv_mean,sv_var)])>0))  else 0
  mBIC = -2*loglik + (2*sum(beta!=0)+2*sum(gamma!=0)-1)*log(n) + logn_mean + logn_var + logn0
  pBIC = -2*loglik + (2*sum(beta!=0)+2*sum(gamma!=0))*log(2*n-3) + 2*log(n)+log(det(t(cbind(1,X[,sv_mean]))%*%InvSigma%*%cbind(1,X[,sv_mean])))
  #return(list('beta'=beta,'gamma'=gamma,'sigma2'=sigma2_0,'b0'=b0,'loglik' = loglik,'BIC'=BIC,'mBIC'=mBIC))
  result = list('beta'=beta,'gamma'=gamma,'sigma2'=sigma2_0,
              'b0'=b0,'loglik' = loglik,
              'BIC'=BIC,'mBIC' = mBIC,'pBIC' = pBIC,
              'fitted.values'=X %*% beta+b0,
              'Sigma'=Sigma)

  return(result)
}

#' @title backward_correction
#' @description This function conduct backward selection on given combination of shifts in optimal value and variance
#' @param tree a phylogenetic tree
#' @param Y a vector of trait values
#' @param alpha the parameter for selection strength
#' @param sv_mean shifts in optimal value
#' @param sv_var shifts in variance
#' @param criterion measurement criterion
#' @param original_model original OU model; optional, default NULL
#' @return
#' \item{OModel}{Fitted OU model}
#' @export

backward_correction = function(tree,Y,alpha,sv_mean,sv_var,criterion='BIC',original_model=NULL){

  if(is.null(original_model)){
    OModel = fit_OU_mean_var(tree,Y,alpha,sv_mean,sv_var)
  }else{
    original_model = OModel
  }
  original_score = OModel[[criterion]][[1]]
  best_score = original_score
  n1 = length(sv_mean)
  n2 = length(sv_var)
  if(n1+n2==0){
    return(OModel)
  }
  after_score_list = rep(0,n1+n2)
  for(i in 1:(n1+n2)){
    if(i <= length(sv_mean)){
      new_Model = fit_OU_mean_var(tree,Y,alpha,setdiff(sv_mean,c(sv_mean[i])),sv_var)
      after_score_list[i] = new_Model[[criterion]][[1]]
    } else{
      new_Model = fit_OU_mean_var(tree,Y,alpha,sv_mean,setdiff(sv_var,c(sv_var[i-n1])))
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
                                setdiff(sv_mean,sv_mean[c(i,remove_list)[c(i,remove_list)<=n1]]),
                                setdiff(sv_var,sv_var[c(i,remove_list)[c(i,remove_list)>n1]-n1])
                                )
   if(new_Model[[criterion]][[1]]<best_score){
      OModel = new_Model
      best_score = new_Model[[criterion]][[1]]
      remove_list = c(remove_list,i)
    }
  }
  return(OModel)
  }

#' @title get_mean_var_shifts_model_selection
#' @description This function provide estimation of shift sizes in optimal value and shifts in variance (including model selection)
#' @param Y a vector of trait values
#' @param tree a phylogenetic tree
#' @param alpha the parameter for selection strength
#' @param lambda1_list a list of degree of shrinkage parameters for shifts in optimal value
#' @param lambda2_list a list of degree of shrinkage parameters for shifts in variance
#' @param criterion the criterion for model selection, default BIC
#' @param max.steps maximum number of steps, default 100
#' @param nfolds number of folds for cross validation for tunning lambda1, default 8
#' @param top_k conduct backward correction on k original models with the best criterion scores, default 10
#' @return
#' \item{best_model}{Fitted OU model}
#' \item{score_summary}{score summary table}
#' @export
#' @import ape
#' @importFrom glmnet glmnet 
#' @importFrom glmnet cv.glmnet
#' @importFrom PIGShift OU.vcv
#' @examples
#' require(ape)
#' require(MASS)
#' require(PIGShift)
#' alpha = 1
#' sigma2_0 = alpha*2
#' true_shifts_var = c(8) #shift configuration in variance
#' size_var = 15 #shift size in variance
#' true_shifts_mean = 0 #shift configuration in mean
#' size_mean = 0 #shift size in mean
#' set.seed(345)
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
#' result = get_mean_var_shifts_model_selection(Y,tree,alpha,NULL,exp(1:10-6),nfolds=3)

get_mean_var_shifts_model_selection = function(Y,tree,alpha,lambda1_list=NULL,lambda2_list,criterion="BIC",max.steps=100,nfolds = 8,top_k=10){

  top_k = top_k
  X = generate_design_matrix(tree)
  Y = as.vector(Y)
  if(alpha == 0){
    V = vcv(tree)
  }else{
    V = OU.vcv(tree,alpha)
  }
  n = nrow(X)
  p = ncol(X)
  tb = node.depth.edgelength(tree)[tree$edge[,1]]
  q = exp(-2*alpha*max(node.depth.edgelength(tree)))/(2*alpha)*(1-exp(2*alpha*tb))
  score_summary = data.frame(matrix(nrow=0,ncol=8))
  names(score_summary) = c('lambda1','lambda2','sv_mean','sv_var',
                            'loglik','BIC','mBIC','pBIC')
  Y1 = Y
  #Y1 = (Y-mean(Y))/sd(Y)
  best_score = Inf
  for(lambda2 in lambda2_list){
    ret_pre = get_mean_var_shifts(Y,tree,alpha,Inf,lambda2,max.steps=max.steps)
    Sigma = V*ret_pre$sigma2 +  (X%*%diag(ret_pre$gamma)%*%t(X))*V+ X%*%diag(ret_pre$gamma*q)%*%t(X)
    InvSigma = solve(Sigma,tol=exp(-100))
    svd_result = svd(InvSigma)
    SqrtInvSigma = svd_result$u %*% diag(sqrt(svd_result$d))  %*% t(svd_result$v)
    YY = SqrtInvSigma %*% Y
    XX = SqrtInvSigma %*% X
    lambda1 = cv.glmnet(XX,YY,lambda=lambda1_list,intercept=FALSE,nfolds = nfolds)$lambda.1se
    ret = get_mean_var_shifts(Y,tree,alpha,lambda1,lambda2,max.steps = max.steps)
    OModel = fit_OU_mean_var(tree,Y1,alpha,ret$sv_mean,ret$sv_var)
    score_summary[nrow(score_summary)+1, ] = c(lambda1,lambda2,
                                              paste(ret$sv_mean,collapse=';'),
                                              paste(ret$sv_var,collapse=';'),
                                              OModel$loglik,
                                              OModel$BIC,
                                              OModel$mBIC,
                                              OModel$pBIC)
    if(OModel[[criterion]][[1]] < best_score){
      best_score = OModel[[criterion]][[1]]
      best_Model = OModel
      best_Model$lambda1 = as.numeric(lambda1)
      best_Model$lambda2 = as.numeric(lambda2)
    }
    }

  score_summary = score_summary[!duplicated(score_summary[,c('sv_mean','sv_var')]),]
  score_summary = cbind(score_summary,NA,NA,NA,NA,NA,NA)
  names(score_summary)[9:14] = c('sv_mean_corrected','sv_var_corrected',
                            'loglik_corrected','BIC_corrected',
                            'mBIC_corrected','pBIC_corrected')

  for(i in 1:min(top_k,nrow(score_summary))){
    ind = order(as.numeric(score_summary[[criterion]]))[i]
    sv_mean = as.numeric(strsplit(score_summary$sv_mean[ind],';')[[1]])
    sv_var = as.numeric(strsplit(score_summary$sv_var[ind],';')[[1]])
    OModel = backward_correction(tree,Y1,alpha,sv_mean,sv_var)
    score_summary[ind, 9:14] = c(paste((1:ncol(X))[OModel$beta!=0],collapse=';'),
                                        paste((1:ncol(X))[OModel$gamma!=0],collapse=';'),
                                        OModel$loglik,
                                        OModel$BIC,
                                        OModel$mBIC,
                                        OModel$pBIC)
    if(OModel[[criterion]][[1]] < best_score){
      best_score = OModel[[criterion]][[1]]
      best_Model = OModel
      best_Model$lambda1 = score_summary$lambda1[ind]
      best_Model$lambda2 = score_summary$lambda2[ind]
    }
  }


  return(list('best_model' = best_Model,'score_summary'=score_summary))
}




