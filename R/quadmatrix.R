#' @title Solving Quadratic Matrix Equations
#'
#' @description Given inputs A,B and C, this package solves the matrix equation A*X^2 - B*X - C = 0.
#'
#' @param A  A Matrix
#' @param B  A Matrix
#' @param C  A Matrix
#'
#'
#'
#'
#' @return Matrix X
#'
#' @examples
#' B = matrix(c(1,3,3,1),nrow = 2,ncol = 2)
#' A = matrix(c(1.5,2,2,1.5),nrow = 2,ncol = 2)
#' C = matrix(c(1.5,1,1,1.5),nrow = 2,ncol = 2)
#' quadmatrix(A,B,C)
#'
#'
#'
#' @export quadmatrix


#options(warn=-1)
#options("getSymbols.warning4.0"=FALSE)
quadmatrix = function(A, B, C) {
  if (!all(dim(A) == dim(B)) &&
      !all(dim(C) == dim(B))) {
    stop("A,B and C does not have the same dimensions")
  }

  if (! matrixcalc::is.square.matrix(A)){
    stop("A has to be square matrix")
  }

  if (! matrixcalc::is.square.matrix(B)){
    stop("B has to be square matrix")
  }

  if (! matrixcalc::is.square.matrix(C)){
    stop("C has to be square matrix")
  }

  p = dim(A)[1]

  E1 = matrix(0, 2 * p, 2 * p)
  E2 = matrix(0, 2 * p, 2 * p)

  E1[1:p, 1:p] = B
  E1[1:p, (p + 1):(2 * p)] = C
  E1[(p + 1):(2 * p), (1:p)]  = diag(p)

  E2[1:p, 1:p] = A
  E2[(p + 1):(2 * p), (p + 1):(2 * p)] = diag(p)

  gen_eigen = geigen::geigen(E1,E2)

  ordered_values = sort(abs(gen_eigen$values), index.return = T)$ix

  eigenmatrix = gen_eigen$vectors[, ordered_values]

  X12 = eigenmatrix[(p + 1):(2 * p), 1:p]
  Delta_mat = diag(gen_eigen$values[ordered_values][1:p], p)


  cov_estimate = X12 %*% Delta_mat %*% solve(X12)
  return(cov_estimate)
}
