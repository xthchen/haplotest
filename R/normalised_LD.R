#' Computation of normalised linkage disequillibrium.
#' 
#' Compute Linkage Disequillibrium(LD) between 2 SNPs, return normalised LD (Lewontin, R. C. 1964) and estimated variance of normalised LD (Zapata, C. et al. 1997). 
#' @param SNP_1 Numeric, SNP position of SNP 1
#' @param SNP_2 Numeric, SNP position of SNP 2
#' @param haplotype Binary numeric matrix containing information of haplotype structure. The columns corresponds to the haplotypes and the rows to the SNPs.
#' @param frequency_matrix Numeric matrix, with each i,j element being the frequency of haplotype i at sequenced time point j.
#' @param allele_frequency Numeric matrix, with each i,j element being the frequency of SNP i at sequenced time point j.
#' @references Lewontin, R. C., (1964), The Interaction of Selection and Linkage. I. General Considerations; Heterotic Models, Genetics 49(1): 49–67.
#' @references Zapata, C., Alvarez, G., Carollo, C., (1997), Approximate variance of the standardized measure of gametic disequilibrium D’, American journal of human genetics 61(3): 771–774.
#' @return A numeric vector with its first element being the normalised LD, second element the estimated variance of normalised LD.
#' @seealso [haplotest()]
#' @export

normalised_LD = function(SNP_1, SNP_2, haplotype, frequency_matrix, allele_frequency){
  
  int_time = ncol(frequency_matrix)
  #position of last timepoint
  p1 = allele_frequency[SNP_1, int_time]
  #allele frequency at position SNP_1
  q1 = allele_frequency[SNP_2, int_time]
  #allele frequency at position SNP_2
  int_SNP = haplotype[c(SNP_1,SNP_2), ]
  #all haplotype only at designated SNP position
  hap_record = c()
  for (i in 1:ncol(haplotype)){
    if(setequal(int_SNP[,i], c(1,1))){
      hap_record = c(hap_record,i)
    }
  }
  #find out which haplotype gives c(1,1)
  if (length(hap_record) == 0){
    D = - p1*q1
  } else {
    D = sum(frequency_matrix[hap_record, int_time]) - p1*q1
  }
  #linkage disequillibrium
  if (D<0){
    Dmax = max(c(-p1*q1), -(1-p1)*(1-q1))
  } else {
    Dmax = min(p1*(1-q1), (1-p1)*q1)
  }
  if (abs(D) < 0.0001){
    D_norm = 0
  } else if (Dmax == 0){
    Dmax = 0.00001
    D_norm = as.numeric(D/Dmax)
  } else {
    D_norm = as.numeric(D/Dmax)
  }
  #normalised LD
  D_var = (p1*q1*(1-p1)*(1-q1)+D*((1-p1)-p1)*((1-q1)-q1)-D^2)/ncol(haplotype)
  
  hap_record2 = c()
  for (i in 1:ncol(haplotype)){
    if(setequal(int_SNP[,i], c(1,0))){
      hap_record2 = c(hap_record2,i)
    }
  }
  
  hap_record3 = c()
  for (i in 1:ncol(haplotype)){
    if(setequal(int_SNP[,i], c(0,1))){
      hap_record3 = c(hap_record3,i)
    }
  }
  
  hap_record4 = c()
  for (i in 1:ncol(haplotype)){
    if(setequal(int_SNP[,i], c(0,0))){
      hap_record4 = c(hap_record4,i)
    }
  }
  
  if (Dmax == (-p1*q1)){
    xi = sum(frequency_matrix[hap_record, int_time])
  } else if (Dmax == (p1*(1-q1))){
    xi = sum(frequency_matrix[hap_record2, int_time])
  } else if (Dmax == ((1-p1)*q1)){
    xi = sum(frequency_matrix[hap_record3, int_time])
  } else {
    xi = sum(frequency_matrix[hap_record4, int_time])
  }
  
  if (abs(D_norm) == 1){
    D_norm_var = 0
  } else if (D_norm == 0){
    D_norm_var = 0
  } else if (D_norm < 0){
    D_norm_var = (1/(ncol(haplotype)*Dmax^2))*((1-abs(D_norm))*(ncol(haplotype)*D_var-abs(D_norm)*Dmax*(p1*(1-q1)+(1-p1)*q1-2*abs(D)))+abs(D_norm)*xi*(1-xi))
  } else {
    D_norm_var = (1/(ncol(haplotype)*Dmax^2))*((1-abs(D_norm))*(ncol(haplotype)*D_var-abs(D_norm)*Dmax*(p1*q1+(1-p1)*(1-q1)-2*abs(D)))+abs(D_norm)*xi*(1-xi))
  }
  return(c(D_norm, abs(D_norm_var)))
}