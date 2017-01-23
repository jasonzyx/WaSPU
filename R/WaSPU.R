#' A weighted Adaptive Sum of Powered Correlation Test (WaSPU) with Individual level data
#'
#' @param Y, an n x 1 response vector
#' @param X, an n x p design matrix
#' @param Z, an n x p nuisance covariace matrix
#' @param pow, power integer candidates, default c(1:6, Inf)
#' @param n.perm, number of permutations to calculate a P-value
#' @return p-values of aSPU, SSU, Sum and UminP tests.
#' @references Xu Z., Xu Gong, Pan W. Integrating genomic and imaging endophenotypes in GWAS
#' @examples
#' library(aSPU)
#' Y = sample(c(1,0), 100, replace = T)
#' X = matrix(sample(c(0,1,2), 1000, replace = T), nrow = 100)
#' Z = matrix(rnorm(200), nrow = 100)
#' weights = runif(10)
#' WaSPU(Y, X, Z, weights, pow = c(1:6, Inf), n.perm = 1e3)
#'
#'
#' @export
#' @importFrom aSPU aSPU


WaSPU = function(Y, X, Z, weights, pow = pow, n.perm = n.perm){
  # test parameter
  # df1 = X; df2 = Y; pow = c(1:8, Inf); B=2
  Xw = X
  for(i2 in 1:ncol(X)){
    Xw[,i2] <- X[,i2]*weights[i2]
  }

  sd = apply(X, 2, sd)
  Xw = as.matrix(Xw[,sd!=0])

  if(ncol(Xw) == 0) stop("after weighting, the design matrix has no variables on all columns")

  pTests3V1 = Tests3V1(Y=Y, X=Xw, Z=Z )

  pSSU = pTests3V1[1]

  pSum = pTests3V1[2]

  pUminP = pTests3V1[3]

  paSPU = aSPU::aSPU(Y,Xw,cov=Z,resample="perm",model="binomial",n.perm=n.perm)$pvs

  return(list(paSPU = paSPU,  pSSU = pSSU, pSum = pSum, pUminP = pUminP))
}
