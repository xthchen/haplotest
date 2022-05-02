#' Selection of non-overlapping haplotype block from filtered haplotype block set
#'
#' The algorithm selects non-overlapping haplotype block from the filtered haplotype block candidate output from the MIG algorithm, by favouring larger haplotype blocks first.
#' @param filt_hap Numeric matrix, the filtered haplotype block set, output from MIG()
#' 
#' @return Numeric matrix, set of non-overlapping haplotype block location, with the first column being the starting SNP position of the candidate, and second column ending SNP position of the candidate. Each row will corresponds to a single candidate.
#'
#' @export

hap_block = function(filt_hap){
  ordered_filt = order(filt_hap[,2]-filt_hap[,1], decreasing = TRUE)
  result_block = c()
  dummy_block = c()
  for (i in ordered_filt){
    if (any(seq(filt_hap[i,1], filt_hap[i,2], 1) %in% dummy_block) == FALSE){
      result_block = rbind(result_block, filt_hap[i,])
      dummy_block = c(dummy_block, seq(filt_hap[i,1], filt_hap[i,2], 1))
    }
  }
  return(result_block)
}