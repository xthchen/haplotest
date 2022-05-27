#' Adapted CMH test if all frequencies estimated with sampling variances.
#'
#' This function performs the adapted CMH test (Spitzer et al. 2019), but in cases where all frequencies estimated with sampling variances. Drift variance and sampling variance are taken into account for all frequencies in this test.
#' @param freq Numeric matrix of frequencies, with the row being the haplotype/SNP, and column being the sequenced time points.
#' @param Ne Numeric vector with length as number of replicates, containing information of Ne (effective population size) at each replicated population. If Ne changes over time, take as input a numeric matrix, with the column being the replicate position, row being the Ne at each sequenced time points.
#' @param gen Numeric vector of sequenced time points.
#' @param cov Numeric matrix of sample size, with the row being the haplotype/SNP, and column being the sequenced time points.
#' @param repli Numeric, specifying the number of replicated populations.
#' @return Numeric vector, p-values from the hypothesis test
#' @references Spitzer, K., Pelizzola, M., Futschik, A., (2019), Modifying the Chi-square and the CMH test for population genetic inference: adapting to over-dispersion, The Annals of Applied Statistics 14(1): 202-220.
#' @seealso [haplotest()]
#' @export

adapted.cmh.test.est.pool = function(freq, Ne, gen, cov, pool, repli){
  k = length(gen)
  #no. timepoints
  repli_gen = seq(1, repli*k, k)
  #starting col of each repli
  tchi_num = c()
  tchi_den = c()
  for (i in 1:repli){
    f_0 = freq[,repli_gen[i]]
    #starting frequencies
    t_T = repli_gen[i]+k-1
    #last timepoint position of current replicate
    f_T = freq[,t_T]
    #ending frequencies
    n_0 = cov[,repli_gen[i]]
    n_T = cov[,t_T]
    #sample sizes
    s_0 = pool[,repli_gen[i]]
    s_T = pool[,t_T]
    #pooled sample sizes

    if (!is.matrix(Ne)){
      drift = rowSums(freq[,(repli_gen[i]:(t_T-1))]*(1-freq[,(repli_gen[i]:(t_T-1))])*
                        matrix(rep(1-(1-1/Ne[i])^(gen[-1]-gen[-k]),nrow(freq)), nrow = nrow(freq), byrow = TRUE))
      tchi_num = cbind(tchi_num, s_0^2*s_T^2*((f_0*(1-f_T)-(1-f_0)*f_T)^2))
      var_0 = s_0*f_0*(1-f_0)*(1+(s_0-1)/n_0)
      var_T = s_T*(f_T*(1-f_T)*(1+(s_T-1)/n_T)+(s_T-1)*((n_T-1)/n_T)*drift)
      tchi_den = cbind(tchi_den, (s_0^2*var_T + s_T^2*var_0))
    }else{
      hm_Ne = floor(2/(1/Ne[-length(Ne[,i]),i]+1/Ne[-1,i]))
      drift = rowSums(freq[,(repli_gen[i]:(t_T-1))]*(1-freq[,(repli_gen[i]:(t_T-1))])*
                        matrix(rep(1-(1-1/hm_Ne)^(gen[-1]-gen[-k]),nrow(freq)), nrow = nrow(freq), byrow = TRUE))
      tchi_num = cbind(tchi_num, s_0^2*s_T^2*((f_0*(1-f_T)-(1-f_0)*f_T)^2))
      var_0 = s_0*f_0*(1-f_0)*(1+(s_0-1)/n_0)
      var_T = s_T*(f_T*(1-f_T)*(1+(s_T-1)/n_T)+(s_T-1)*((n_T-1)/n_T)*drift)
      tchi_den = cbind(tchi_den, (s_0^2*var_T + s_T^2*var_0))
    }
  }
  tchi = rowSums(tchi_num)/rowSums(tchi_den)
  res = pchisq(tchi, df = 1, lower.tail = FALSE)
  return(res)
}

