#' Adapted CMH test if all frequencies are known.
#'
#' This function performs the adapted CMH test (Spitzer et al. 2019), but in cases where all frequencies are known, and the only variance of data comes from genetic drift. Drift variance is the only variance included in this test.
#' @param freq Numeric matrix of frequencies, with the row being the haplotype/SNP, and column being the sequenced time points.
#' @param Ne Numeric vector with length as number of replicates, containing information of Ne (effective population size) at each replicated population. If Ne changes over time, take as input a numeric matrix, with the column being the replicate position, row being the Ne at each sequenced time points.
#' @param gen Numeric vector of sequenced time points.
#' @param repli Numeric, specifying the number of replicated populations.
#' @return Numeric vector, p-values from the hypothesis test
#' @references Spitzer, K., Pelizzola, M., Futschik, A., (2019), Modifying the Chi-square and the CMH test for population genetic inference: adapting to over-dispersion, The Annals of Applied Statistics 14(1): 202-220.
#' @seealso [haplotest()]
#' @export

adapted.cmh.test.true = function(freq, Ne, gen,repli){
  k = length(gen)
  #no. timepoints
  repli_gen = seq(1,repli*k,k)
  #starting col of each repli
  tchi_num = c()
  tchi_den = c()
  if (length(Ne) < repli){
    return("Ne is not available for all replicates")
  }
  for (i in 1:repli){
    f_0 = freq[,repli_gen[i]]
    #starting frequencies
    t_T = repli_gen[i]+k-1
    #last timepoint
    f_T = freq[,t_T]
    #ending frequencies
    if (!is.matrix(Ne)){
      drift = rowSums(freq[,(repli_gen[i]:(t_T-1))]*(1-freq[,(repli_gen[i]:(t_T-1))])*
                        matrix(rep(1-(1-1/Ne[i])^(gen[-1]-gen[-k]),nrow(freq)), nrow = nrow(freq), byrow = TRUE))
      tchi_num = cbind(tchi_num,(Ne[i]^4*(f_0*(1-f_T)-(1-f_0)*f_T)^2))
      tchi_den = cbind(tchi_den,(Ne[i]^3*(Ne[i]-1)*drift))  
    }else{
      hm_Ne = floor(2/(1/Ne[-length(Ne[,i]),i]+1/Ne[-1,i]))
      drift = rowSums(freq[,(repli_gen[i]:(t_T-1))]*(1-freq[,(repli_gen[i]:(t_T-1))])*
                        matrix(rep(1-(1-1/hm_Ne)^(gen[-1]-gen[-k]),nrow(freq)), nrow = nrow(freq), byrow = TRUE))
      tchi_num = cbind(tchi_num,(Ne[1,i]^2*Ne[nrow(Ne),i]^2*(f_0*(1-f_T)-(1-f_0)*f_T)^2))
      tchi_den = cbind(tchi_den,(Ne[1,i]^2*Ne[nrow(Ne),i]*(Ne[nrow(Ne),i]-1)*drift))
    }
  }
  tchi = rowSums(tchi_num)/rowSums(tchi_den)
  res = pchisq(tchi, df = 1, lower.tail = FALSE)
  return(res)
}