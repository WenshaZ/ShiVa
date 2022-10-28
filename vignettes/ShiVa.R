### R code from vignette source 'ShiVa.Rnw'

###################################################
### code chunk number 1: ShiVa.Rnw:19-48
###################################################
library(ape)
library(MASS)
library(PIGShift)
library(ShiVa)
alpha = 1
sigma2_0 = alpha*2
true_shifts_var = c(2) #shift configuration in variance
size_var = 10 #shift size in variance
true_shifts_mean = 0 #shift configuration in mean
size_mean = 0 #shift size in mean
set.seed(123)
# generate a tree
tree = rcoal(80)
tree$edge.length = tree$edge.length/max(node.depth.edgelength(tree))
# generate design matrix
X = generate_design_matrix(tree,type='simpX')
V = OU.vcv(tree,alpha)
tb = node.depth.edgelength(tree)[tree$edge[,1]]
q = exp(-2*alpha*max(node.depth.edgelength(tree)))/(2*alpha)*(1-exp(2*alpha*tb))
# simulate shifts in variance
gamma = rep(0,ncol(X))
gamma[true_shifts_var]  = size_var
Sigma = V*sigma2_0 +  X%*%diag(gamma)%*%t(X)*V+ X%*%diag(gamma*q)%*%t(X)
# simulate shifts in mean
beta = rep(0,ncol(X))
beta[true_shifts_mean] = size_mean
# generate simulated data
eps = mvrnorm(n = 1, mu = rep(0,nrow(X)), Sigma = Sigma)
Y = X%*%beta+eps


###################################################
### code chunk number 2: ShiVa.Rnw:59-64
###################################################
result = get_mean_var_shifts(Y, tree, 1, 1, 0.02)
sv_mean = (1:ncol(X))[result$beta!=0]
sv_var = (1:ncol(X))[result$gamma!=0]
sv_mean # shifts in optimal values
sv_var # shifts in variance


###################################################
### code chunk number 3: ShiVa.Rnw:75-77
###################################################
OModel = fit_OU_mean_var(tree, Y ,1, sv_mean, sv_var)
OModel$gamma[sv_var]


###################################################
### code chunk number 4: ShiVa.Rnw:87-95
###################################################
result = get_mean_var_shifts_model_selection(Y,tree,1,NULL,exp(1:10-6))
sv_mean = (1:ncol(X))[result$beta!=0]
sv_var = (1:ncol(X))[result$gamma!=0]
sv_mean # shifts in optimal values
sv_var # shifts in variance
result$gamma[sv_var] # shift sizes
result$lambda1
result$lambda2


