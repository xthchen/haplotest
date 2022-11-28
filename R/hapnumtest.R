#' Iterative testing for the number of selected haplotypes.
#'
#' This is the post-hoc test for testing for the number of selected haplotypes. The test should only be performed if haplotype based test have shown evidence of selection.
#' @param frequency_matrix Numeric matrix, with each i,j element being the frequency of haplotype i at sequenced time point j.
#' @param deltat Numeric, number of generations between each pair of time points of interests
#' @param Ne Numeric vector with length as number of replicates, containing information of Ne (effective population size) at each replicated population. If Ne changes over time, take as input a numeric matrix, with the column being the replicate position, row being the Ne at each sequenced time points.
#' @param repli Numeric, specifying the number of replicated populations.
#' @param p_combine_method Factor, method of pvalue combination, can be "omnibus" (Futschik, A. et al. 2019), "harmonic" (Wilson, D.J. (2019)), "vovk" (Vovk, V. et al. 2018), "bonferroni", "BH" (Benjamini & Hochberg 1995).
#' @param seed setting seed of the run
#' @return A number specifying the number of selected haplotypes.
#' @seealso [haplotest()]
#' @export
#'
#' @import omnibus
#' @import harmonicmeanp
#'
#' @examples
#' #We show here an example for the test for the number of selected haplotypes
#' #Suppose we have the haplotype frequency matrix hap_freq, sequenced at every 10 generations, with effective population size 1000 and 3 replicate populations, using omnibus test as multiple testing correction, the number of selected haplotypes is computed by:
#' hap_num = hapnumtest(hap_freq, p_combine_method = "omnibus", deltat = 10, Ne = rep(1000,3), repli = 3)



hapnumtest = function(frequency_matrix, p_combine_method = "omnibus", deltat = 10, Ne = 1000,
                      repli = 1, seed = 2022){
  set.seed(seed)
  t = (ncol(frequency_matrix)/repli)
  nrm_fq = frequency_matrix
  freq_mod = frequency_matrix
  if (!is.matrix(Ne)){
    new_ne = matrix(rep(Ne, ncol(frequency_matrix)), ncol = repli)
  }else{
    new_ne = Ne
  }
  gen = seq(0,(t-1)*deltat,deltat)
  repli_gen = seq(1,repli*t,t)
  pval = 0
  curr_sel = 1
  while (pval < 0.05){
    if (curr_sel < (nrow(frequency_matrix)-1)){
      max_diff_repli = rep(0,nrow(freq_mod))
      for(z in 1:repli){
        f_0 = nrm_fq[,repli_gen[z]]
        #starting frequencies
        t_T = repli_gen[z]+t-1
        #last timepoint
        f_T = nrm_fq[,t_T]
        #harmonic means
        hm_Ne = floor(2/(1/new_ne[-length(new_ne[,z]),z]+1/new_ne[-1,z]))
        hm_Ne = hm_Ne[repli_gen[z]:(repli_gen[z]+t-2)]

        diff = f_T - f_0
        drift = rowSums(nrm_fq[,(repli_gen[z]:(t_T-1))]*(1-nrm_fq[,(repli_gen[z]:(t_T-1))])*
                          matrix(rep(1-(1-1/hm_Ne)^(gen[-1]-gen[-t]),nrow(nrm_fq)), nrow = nrow(nrm_fq), byrow = TRUE))
        max_diff_repli = max_diff_repli + diff/sqrt(drift)
      }
      freq_mod = nrm_fq[-which.max(max_diff_repli),]

      tot_ne = c()
      for(z in 1:repli){
        ne_hold = c()
        for (v in 1:t){
          ne_hold = c(ne_hold, floor(sum(freq_mod[,t*(z-1)+v])*new_ne[v,z]))
        }
        tot_ne = cbind(tot_ne, ne_hold)
      }
      new_ne = tot_ne
      #new ne after frequency removal


      for (x in 1:ncol(freq_mod)){
        if (sum(freq_mod[,x]) == 0){
          freq_mod[,x] = rep(1/nrow(freq_mod), nrow(freq_mod))
        }else{
          freq_mod[,x] = (1/sum(freq_mod[,x]))*freq_mod[,x]
        }
      }
      #normalisation of frequencies after removal

      nrm_fq = freq_mod
      for(z in 1:repli){
        for(v in 1:t)
          if(new_ne[v,z] == 0){
            new_ne[v,z] = 1
          }
      }
      hap_re = adapted.cmh.test.true(nrm_fq, Ne = new_ne, gen = seq(0,(t-1)*deltat,deltat), repli = repli)
      if (any(is.nan(hap_re))){
        hap_re = rep(1,length(hap_re))
      }
      if(p_combine_method == "omnibus"){
        pval = omnibus.test(hap_re,method = "log.p", N.sim = 1000)
      }else if (p_combine_method == "harmonic"){
        pval = hmp.stat(hap_re)
      }
      if (pval>0.05){
        return(curr_sel)
      }else{
        curr_sel = curr_sel+1
      }
    }else{
      curr_sel = nrow(frequency_matrix)-1
      pval = 1
    }
  }
}
