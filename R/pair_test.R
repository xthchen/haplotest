#' Apply pairwise test to all pairs of haplotypes.
#' 
#' Pairwise post hoc test. Should only be used if haplotype based test have shown to be significant.
#' @param freq Numeric matrix, with each i,j element being the frequency of haplotype i at sequenced time point j.
#' @param Ne Numeric vector with length as number of replicates, containing information of Ne (effective population size) at each replicated population. If Ne changes over time, take as input a numeric matrix, with the column being the replicate position, row being the Ne at each sequenced time points.
#' @param repli Numeric, specifying the number of replicated populations.
#' @param tdelta Numeric, number of generations between each pair of time points of interests
#' @return A numeric vector of p-values from all the pairwise test after B&H corrections. The pairwise test is sorted such that its between haplotypes :1-2, 1-3, 1-4,...,2-3,2-4,...,3-4,...
#'
#' @export




pair_test = function(freq, Ne, repli, tdelta){
  nhap = nrow(freq)
  t = ncol(freq)/repli
  
  pval = c()
  if (!is.matrix(Ne)){
    new_ne = matrix(rep(Ne, t), ncol = repli)
  }else{
    new_ne = Ne
  }
  #Put Ne into matrix if not
  for(i in 1:(nhap-1)){
    for(j in (i+1):nhap){
      freq_int = freq[c(i,j),]
      #the haplotype frequency pair
      ne_all = c()
      for (v in 1:repli){
        ne_norm = floor(new_ne[,v]*colSums(freq_int)[(t*(v-1)+1):(t*v)])
        ne_norm[ne_norm == 0] = 2
        ne_all = cbind(ne_all, ne_norm)
      }
      if (length(which(colSums(freq_int) == 0)) > 0){
        freq_int[,which(colSums(freq_int) == 0)] = 0.5
        freq_int[freq_int == 0] = 10^(-4)
      }
      freq_norm = freq_int[1,]*(1/colSums(freq_int))
      #normalising freq
      freq_norm[is.na(freq_norm)] = 0.5
      #replace Na elements
      freq_norm = rbind(freq_norm, 1-freq_norm)
      #new haplotype matrix for the pair
      pval_par = adapted.cmh.test.true(freq = freq_norm, Ne = ne_all, gen = seq(0,(t-1)*tdelta,tdelta), repli = repli)
      if (is.na(pval_par)){
        print(list(freq_norm,ne_all))
      }
      pval = c(pval,pval_par[1])
    }
  }
  return (p.adjust(pval, method = "BH"))
}

