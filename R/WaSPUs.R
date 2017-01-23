#' A weighted Adaptive Sum of Powered Correlation Test (WaSPU) with summary statistics
#'
#' @param Zstat, a p x 1 Z statistics corresponding to a gene
#' @param corSNP, a p x p correlation matrix of SNPs based on a reference panel
#' @param weights, a p x 1 weight vector corresponding to the SNPs
#' @param pow, power integer candidates, default c(1:6, Inf)
#' @param n.perm, number of permutations to calculate a P-value
#' @return p-values of aSPU, SSU, Sum and UminP tests.
#' @references Xu Z., Xu Gong, Pan W. Integrating genomic and imaging endophenotypes in GWAS
#' @examples
#' library(aSPU)
#' set.seed(123)
#' Zstat = runif(10,-2,2)
#' corSNP = 0.1 + diag(0.9,10)
#' weights = runif(10)
#' WaSPUs(Zstat, corSNP, weights, pow = c(1:6, Inf), n.perm = 1e3)
#'
#'
#' @export
#' @importFrom aSPUs aSPU


WaSPUs = function(Zstat, corSNP, weights, pow = c(1:6, Inf), n.perm = 1e3){
  # test parameter
  # df1 = X; df2 = Y; pow = c(1:8, Inf); B=2
  weight_diag = diag(as.vector(weights))
  Zstat.w <- weight_diag %*% Zstat
  corSNP.w <- weight_diag %*% corSNP %*% t(weight_diag)

  pSum = Sum(U=Zstat.w, CovS=corSNP.w)
  pSSU = SumSqU(U=Zstat.w, CovS=corSNP.w)
  pUminP = UminP(U=Zstat.w, CovS=corSNP.w)


  paSPU <- aSPU::aSPUs(Zstat.w, corSNP.w,pow = pow, n.perm = n.perm)$pvs


  return(list(paSPU = paSPU,  pSSU = pSSU, pSum = pSum, pUminP = pUminP))
}
