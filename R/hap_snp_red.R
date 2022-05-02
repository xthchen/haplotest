#' Combination of haplotypes
#' 
#' Combine haplotype frequencies according to the reduced haplotype structural matrix created through SNP based test/haplotype blocks.
#' @param red_snp Binary matrix, the reduced haplotype structural matrix with reduced SNPs.
#' @param frequency_matrix Numeric matrix, with each i,j element being the frequency of haplotype i at sequenced time point j.
#' 
#' @return A reduced haplotype frequency matrix (if haplotypes are combined due to reduced SNPs). Remains unchanged from the input frequency_matrix if no combination of haplotypes occurred.
#'
#' @export

hap_snp_red = function(red_snp, frequency_matrix){
  #red_snp, haplotype matrix with only chosen snps
  #frequency_matrix, the original frequency matrix
  all_hap = list()
  #each element of the list contain a vector of identical haplotype position
  passed_hap = c(0)
  #haplotypes that have already been considered
  
  for(k in 1:ncol(red_snp)){
    if (k %in% passed_hap){
      #do nth, skip to next k
    } else {
      int_hap = c()
      #int_hap records all j haplotype that are identical with haplotype k
      for(j in (1:ncol(red_snp))[-k]){
        if (identical(red_snp[,k], red_snp[,j])){
          int_hap = c(int_hap, k, j)
        }
      }
      int_hap = unique(int_hap)
      if (length(int_hap) > 1){
        #if there exists hap to combine
        passed_hap = c(passed_hap, int_hap)
        passed_hap = unique(passed_hap)
        #record hap pos that have already been processed
        all_hap = append(all_hap, list(int_hap))
        #new item for the list
      }
    }
  }
  if(length(all_hap)>0){
    #more than 1 item in combination list
    rm_hap = c()
    #position of removed haplotype
    red_frequency_matrix = frequency_matrix
    for (k in 1:length(all_hap)){
      red_frequency_matrix[all_hap[[k]][1],] = colSums(red_frequency_matrix[all_hap[[k]],])
      #combine frequency of all other haplotype to first one
      rm_hap = c(rm_hap, all_hap[[k]][-1])
      #record all haplotype position that are being combined
    }
    red_frequency_matrix = red_frequency_matrix[-rm_hap,]
  } else {
    red_frequency_matrix = frequency_matrix
  }
  return(red_frequency_matrix)
}