#' Vovk based multiple testing correction.
#' 
#' Multiple testing combination method based on Vovk et al. 2018.
#' @param pval Numeric vector of p-values used for combination.
#' @param r Numeric, parameter for Vovk's methods.
#' @references Vovk, V., Wang., R., (2018), Combining p-values via averaging, Biometrika 107(4): 791–808.
#' @return Numeric, a single pvalue.
#'
#' @export

merging_function = function(pval, r){
  #pval is a vector of pvalues as input
  #r the parameter
  #as output a single combined pvalue
  K = length(pval)
  if (r != 0){
    M = (sum(pval^r)/K)^(1/r)
  } else {
    M = prod(pval)^(1/K)
  }
  if (r == Inf){
    p = max(pval)
  }else if (r == -Inf){
    p = min(pval)*K
  }else if (r == -1){
    p = exp(1)*log(K)*M
  }else if (r == 0){
    p = exp(1)*M
  }else if (r < -1){
    p = K^(1+1/r)*M*r/(r+1)
  }else if (r > (K-1)){
    p = K^(1/r)*M
  }else{
    p = (r+1)^(1/r)*M
  }
  return(p)
}