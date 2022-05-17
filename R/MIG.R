#' MIG algorithm for computation of haplotype blocks
#'
#' Algorithm for computation of haplotype blocks, from (Taliun, D. et al. 2014). 
#' @param haplotype Binary numeric matrix containing information of haplotype structure. The columns corresponds to the haplotypes and the rows to the SNPs.
#' @param frequency_matrix Numeric matrix, with each i,j element being the frequency of haplotype i at sequenced time point j.
#' @param allele_frequency Numeric matrix, with each i,j element being the frequency of SNP i at sequenced time point j.
#' @param LDU Numeric, parameter compared to upper 90 % confidence interval for normalised LD variance, chosen in line with Gabriel, S. B. et al. (2002)
#' @param LDL Numeric, parameter compared to lower 90 % confidence interval for normalised LD variance, chosen in line with Gabriel, S. B. et al. (2002)
#' @param EHR Numeric, parameter for showing strong evidence of historical recombination, chosen in line with Gabriel, S. B. et al. (2002)
#' @param d Numeric, parameter regarding haplotype block, chosen in line with Gabriel, S. B. et al. (2002)
#' @references Gabriel, S. B., Schaffner, S. F., Nguyen, H., Moore, J. M., Roy, J., Blumenstiel, B., Higgins, J., DeFelice, M., Lochner, A., Faggart, M., Liu-Cordero, S. N., Rotimi, C., Adeyemo, A., Cooper, R., Ward, R., Lander, E. S., Daly, M. J., Altshuler, D., (2002), The structure of haplotype blocks in the human genome, Science 296(5576): 2225–2229.
#' @references Taliun, D., Gamper, J., Pattaro, C., (2014), Efficient haplotype block recognition of very long and dense genetic sequences, BMC bioinformatics 15: 10.
#' @return Numeric matrix of filtered haplotype block candidates, with the first column being the starting SNP position of the candidate, and second column ending SNP position of the candidate. Each row will corresponds to a single candidate.
#' @seealso [haplotest()]
#' @export

MIG = function(haplotype, frequency_matrix, allele_frequency, LDU = 0.98, LDL = 0.7, EHR = 0.9, d = 0.95){
  w_matrix = matrix(0,nrow(haplotype),nrow(haplotype))
  diag(w_matrix) = 1-d
  H = c()
  new_b = 1
  n = nrow(haplotype)
  W = rep(0,nrow(haplotype))
  for (j in 2:n){
    s = 0
    b = new_b
    new_b = j
    if ((j-1) < b){
      k = b
    }else{
      k = j
    }
    for (i in ((k-1):b)){
      if (w_matrix[i,j] == 0){
        D = normalised_LD(i,j,haplotype, frequency_matrix, allele_frequency)
        CU = D[1]+qnorm(0.95,0,1)*sqrt(D[2])/ncol(haplotype)
        #Upper 90% confidence interval of D'
        CL = D[1]-qnorm(0.95,0,1)*sqrt(D[2])/ncol(haplotype)
        #Lower 90% confidence interval of D'
        if ((CU > LDU) & (CL > LDL)){
          w = 1-d
        }else if (CU < EHR){
          w = -d
        }else{
          w = 0
        }
      }
      w_matrix[i,j] = w
      w_matrix[j,i] = w
      s = s+w
      W[i] = W[i]+s
      if ((w == 1-d) & W[i] >= 0){
        H = rbind(H,c(i,j))
      }
      w_sum = sum(w_matrix[i:j,i:j]*upper.tri(w_matrix[i:j,i:j], diag = TRUE)) - 1 + d
      w_max = w_sum + (1-d)*((k-i+1)+(n-i))*(n-k)/2
      if (w_max >= 0){
        new_b = i
      }
    }
  }
  return(H)
}