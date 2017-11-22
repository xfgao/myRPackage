#########################################################
## Stat 202A - Homework 6
## Author: Xiaofeng Gao
## Date: 11/15/2017 
## Description: This script implements logistic regression
## using iterated reweighted least squares using the code 
## we have written for linear regression based on QR 
## decomposition
#########################################################

#############################################################
## INSTRUCTIONS: Please fill in the missing lines of code
## only where specified. Do not change function names, 
## function inputs or outputs. You can add examples at the
## end of the script (in the "Optional examples" section) to 
## double-check your work, but MAKE SURE TO COMMENT OUT ALL 
## OF YOUR EXAMPLES BEFORE SUBMITTING.
##
## Very important: Do not use the function "setwd" anywhere
## in your code. If you do, I will be unable to grade your 
## work since R will attempt to change my working directory
## to one that does not exist.
#############################################################


################ Optional ################ 
## If you are using functions from Rcpp, 
## uncomment and adjust the following two 
## lines accordingly. Make sure that your 
## cpp file is in the current working
## directory, so you do not have to specify
## any directories in sourceCpp(). For 
## example, do not write
## sourceCpp('John_Doe/Homework/Stats202A/QR.cpp')
## since I will not be able to run this.

# library(Rcpp)
# Rcpp::sourceCpp('QR_Logit.cpp')

##################################
## Function 1: QR decomposition ##
##################################
#' QR decomposition
#'
#' @param A, an n x m matrix
#' @return a list with Q.transpose and R, Q is an orthogonal n x n matrix,
#' R is an upper triangular n x m matrix, A = Q %*% R
#' @export
#' @examples 
#' myQR(A)
myQR <- function(A){
  
  n <- nrow(A)
  m <- ncol(A)
  Q <- diag(n)
  R <- A
  
  for(k in 1:(m - 1)){
    x      <- rep(0, n)
    x[k:n] <- R[k:n, k]
    s      <- -1 * sign(x[k])
    v      <- x
    v[k]   <- x[k] - s * norm(x, type = "2")
    u      <- v / norm(v, type = "2")
    
    R <- R - 2 * u %*% t(u) %*% R
    Q <- Q - 2 * u %*% t(u) %*% Q
    
  }
  return(list("Q" = t(Q), "R" = R))
  
}
#################################
## Function 2: Sweep operation ##
#################################
#' Perform a SWEEP operation on A with the pivot element A[m,m].
#'
#' @param A, an square matrix; m, the pivot element index.
#' @return a swept matrix
#' @export
#' @examples 
#' mySweep(A,m)
mySweep <- function(A, m){
  n <- dim(A)[1]
  for (k in 1:m)
  { 
    for (i in 1:n) 
      for (j in 1:n)
        if (i!=k  & j!=k)
          A[i,j] <- A[i,j] - A[i,k]*A[k,j]/A[k,k]
        for (i in 1:n)
          if (i!=k)
            A[i,k] <- A[i,k]/A[k,k]
          for (j in 1:n)
            if (j!=k)
              A[k,j] <- A[k,j]/A[k,k]
            A[k,k] <- - 1/A[k,k]
  }
  return(A)
}

###############################################
## Function 3: Linear regression based on QR ##
###############################################
#' Perform the linear regression of Y on X
#'
#' @param X is an n x p matrix of explanatory variables
#' @param Y is an n dimensional vector of responses
#' @param pFlag is whether to consider intercept (default 0)
#' @return least squares solution vector beta "beta" and standard error vector "std"
#' @export
#' @examples 
#' myLM(X,Y,1)
myLM <- function(X, Y, pFlag=0){
  n = dim(X)[1]
  p = dim(X)[2]
  if (pFlag == 1) {    
    p = p+1
    Z = cbind(rep(1,n),X,Y)
  }
  else  {
    Z = cbind(X,Y)
  }
  QR = myQR(Z)
  Q = QR[["Q"]]
  R = QR[["R"]]
  R1 = R[1:(p),1:(p)]
  Y1 = R[1:(p),(p+1)]
  beta_ls = solve(R1)%*%Y1  
  rss = sum(R[(p+1):n,p+1]^2)
  std = sqrt(diag(rss/(n-p) * solve(t(X)%*%X)))
  
  return(list("beta" = beta_ls, "std" = std))
}

##################################
## Function 4: PCA based on QR  ##
##################################
#' Perform PCA on A
#'
#' @param A is a square matrix, numIter (default 1000)
#' @return a list with D and V
## D is a vector of eigenvalues of A
## V is the matrix of eigenvectors of A
#' @export
#' @examples 
#' myEigen_QR(A,1000)
myEigen_QR <- function(A, numIter = 1000){
  
  n = dim(A)[1]
  V = matrix(rnorm(n*n),nrow=n)
  for (i in 1:numIter) {
    QR = myQR(V)
    Q = QR[["Q"]]
    R = QR[["R"]]
    V = A%*%Q
  }

  return(list("D" = diag(R), "V" = Q))
}

######################################
## Function 5: sigmoid function  ##
######################################
#' Expit/sigmoid function
#' @param x is a number
#' @return (1 / (1 + exp(-x)))
#' @export
expit <- function(x){
  return (1 / (1 + exp(-x)))
}

######################################
## Function 6: Logistic regression  ##
######################################
#' Perform logistic regression of Y on X
#'
#' @param X is an n x p matrix of explanatory variables
#' @param Y is an n dimensional vector of responses
#' @return logistic regression solution vector beta "beta" and standard error vector "std"
#' @export
#' @examples 
#' myLogistic(X,Y)
myLogistic <- function(X, Y){
  
  n = dim(X)[1]
  p = dim(X)[2]
  beta <- matrix(rep(0,p),nrow=p)
  eps = 1e-6
  err = 1
  while (err>eps) {
    eta = X%*%beta
    pr = expit(eta)
    w = pr*(1-pr)
    z = eta+(Y-pr)/w
    sw = sqrt(w)
    mw = matrix(rep(sw,p),nrow=n)
    xwork = mw*X
    ywork = sw*z
    lm = myLM(xwork,ywork)
    beta_new = lm$beta
    # std = lm$std
    # beta_new = t(myLinearRegressionC(xwork,ywork))
    err = sum(abs(beta_new-beta))
    # print(err)
    beta = beta_new
  }
  std = matrix(rep(0,p*p),nrow=p)
  for (i in (1:n)) {
    std = std+w[i]*X[i,]%*%t(X[i,])
  }
  std = sqrt(diag(solve(std)))
  ## Function returns the logistic regression solution vector
  return(list("beta" = beta, "std" = std))  
}

##################################
## Function 7: Ridge regression ##
##################################
#' Perform ridge regression of Y on X using lambda
#'
#' @param X is an n x p matrix of explanatory variables
#' @param Y is an n dimensional vector of dependent variables
#' @param lambda: regularization parameter (lambda >= 0)
#' @return beta, the (p+1) ridge regression solution
#' @export
#' @examples 
#' myRidge(X,Y,lambda)
myRidge <- function(X, Y, lambda){
  n = dim(X)[1]
  p = dim(X)[2]
  Z = cbind(rep(1, n), X, Y)
  A = t(Z) %*% Z
  D = diag(rep(lambda, p+2))
  D[1, 1] = 0
  D[p+2, p+2] = 0
  A = A + D
  S = mySweep(A, p+1)
  beta_ridge = S[1:(p+1), p+2]
  
  return(beta_ridge)
  
}

####################################################
## Function 8: Piecewise linear spline regression ##
####################################################
#' Perform spline regression of Y on X using lambda and p.
#'
#' @param x: An n x 1 vector or matrix of explanatory variables.
#' @param Y: An n x 1 vector of dependent variables. Y can also be a 
#' matrix, as long as the function works.
#' @param lambda: regularization parameter (lambda >= 0)
#' @param p: Number of cuts to make to the x-axis.
#' @return 
#' a list containing two elements:
#' The first element of the list is the spline regression
#' beta vector, which should be p + 1 dimensional.
#' The second element is y.hat, the predicted Y values
#' using the spline regression beta vector. This 
#' can be a numeric vector or matrix.
#' @export
#' @examples 
#' mySpline(x,Y,lambda,p)
mySpline <- function(x, Y, lambda, p = 100){
  
  n = length(x)
  xs = sort(x)
  kmin = min(xs)
  kmax = max(xs)
  delta = (kmax-kmin)/p
  k = seq((kmin),(kmax-delta),delta)
  X = c()
  for (i in 1:p) {
    X = cbind(X,(xs>k[i])*(xs-k[i]))
  }
  beta_spline = myRidge(X,Y,lambda)
  Yhat = cbind(rep(1,n),X)%*%beta_spline
  
  output <- list(beta_spline = beta_spline, predicted_y = Yhat)
  return(output)
  
}

#####################################
## Function 9: Lasso solution path ##
#####################################
#' Find the lasso solution path for various values of 
#' the regularization parameter lambda.
#' @param X: n x p matrix of explanatory variables.
#' @param Y: n dimensional response vector
#' @param lambda_all: Vector of regularization parameters.
#' 
#' @return 
#' Returns a matrix containing the lasso solution vector 
#' beta of size (p+1) x length(lambda_all) for each 
#' regularization parameter.
#' 
#' @export
#' @examples 
#' myLasso(X, Y, lambda_all)
myLasso <- function(X, Y, lambda_all){
  n = dim(X)[1]
  p = dim(X)[2]
  X = cbind(rep(1,n),X)
  p = p+1
  L = length(lambda_all)
  lambda_all = sort(lambda_all,decreasing=TRUE)
  T1 = 10
  beta = matrix(rep(0, p), nrow = p)
  beta_all = matrix(rep(0, p*L), nrow = p)
  R = Y
  ss = rep(0, p) 
  for (j in 1:p)
    ss[j] = sum(X[, j]^2)
  err = rep(0, L)
  for (l in 1:L)
  {
    lambda = lambda_all[l]
    for (t in 1:T1)
    { 
      for (j in 1:p)
      { 
        db = sum(R*X[, j])/ss[j]
        b = beta[j]+db
        if (j == 1)
          b = sign(b)*max(0, abs(b))
        else
          b = sign(b)*max(0, abs(b)-lambda/ss[j]) 
        db = b - beta[j]
        R = R - X[, j]*db
        beta[j] = b
      }
    }
    beta_all[,l] = beta
  }
  return(beta_all)
}