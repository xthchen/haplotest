#' Simulation of estimated haplotype frequencies with some sample size
#'
#' Noisy estimated haplotype frequencies are simulated assuming multinomial distribution with some sample size.
#' @param freq Numeric matrix, with each i,j element being the frequency of haplotype i at sequenced time point j.
#' @param samsize Numeric vector or numeric. If numeric vector, a samplesize vector the same length as the ncol(freq). If numeric, the mean sample size, a poisson distribution will be used to determine the samplesize at different timepoints.
#' @return A list, with the first item being the noisy frequency matrix, the second item the samplesize matrix.
#' @seealso [haplotest()]
#' @export


hap_err_sim = function(freq, samsize){
  if (length(samsize) == 1){
    samsize = rpois(ncol(freq), lambda = samsize)
  }
  freq_err = c()
  for (i in 1:ncol(freq)){
    freq_err = cbind(freq_err, rmultinom(n=1, size = samsize[i],prob = freq[,i])/samsize[i])
    #freq_err = cbind(freq_err, rowSums(rmultinom(n=ncol(freq), size = samsize[i],prob = freq[,i])/samsize[i])/ncol(freq))

  }
  return(list(freq_err, samsize))
}
