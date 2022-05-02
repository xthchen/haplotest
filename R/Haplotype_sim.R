#' Simulate a haplotype structure matrix
#'
#' Simulate a haplotype structure matrix through random sampling of binary argutments, where the rows correspond to the SNPs and the columns to the haplotypes.
#' @param nsnp Integer, number of SNPs
#' @param nhap Integer, number of haplotypes
#'
#' @return Returns a numeric matrix, with all entries being binary, where the rows correspond to the SNPs and the columns to the haplotypes.
#'
#' @export
#' 
Haplotype_sim = function(nsnp = 500,nhap = 5){
  if(nsnp == round(nsnp) & nhap == round(nhap) & nsnp>0 & nhap>0){
    return(matrix(round(runif(nsnp*nhap)), nrow = nsnp, ncol = nhap))
  } else {
    return("Both the number of SNPs (nhap) and the number of haplotypes (nhap) need to be positive integer numbers.")
  }
  
}