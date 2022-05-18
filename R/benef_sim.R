#' Generate the position of benefitial allele and selection strength.
#'
#'
#' @param haplotype Binary numeric matrix containing information of haplotype structure. The columns corresponds to the haplotypes and the rows to the SNPs.
#' @param n_benef Non-negative integer, containing information of number of benefitial alleles. 0 corresponds to the scenario of no benifitial alleles.
#' @param min Numeric, minimum selection strength of the benefitial allele, must be between 0 and max.
#' @param max Numeric, maximum selection strength of the benefitial allele, must be between min and 1.
#' @param fix_sel Numeric or numeric vector, whether to fix the number of selected haplotypes to some number. 0 means no fixing, any integer greater than 0 fixes the number of selected haplotypes to the number, any integer lower than 0 means no such number of selected haplotypes (i.e. -1 means no cases of 1 selected haplotypes). Can take numeric vector input, for example c(1,2) will fix the number of selected haplotypes to 1 and 2, c(-1,-2) will disallow number of selected haplotypes as 1 and 2.
#' @param repli Numeric, number of replicated populations.
#' @param diff_sel_str Logical statement, given there are replicates, whether to allow difference in selection strength between replicates.
#' @return A list, with the first item being a numeric vector of selective strength for all SNPs, the second item being the position of benefitial alleles, the third item being the number of selected haplotypes.
#' @seealso [haplotest()]
#' @export

benef_sim = function(haplotype, n_benef = 1, min = 0.05, max = 0.05, fix_sel = 0, repli = 1,
                     diff_sel_str = FALSE){
  if(is.null(haplotype)){
    warning("The haplotype structure matrix is missing. A random haplotype configuration over 1000 SNPs and 5 haplotypes is simulated")
  }
  if (n_benef != round(n_benef)|repli != round(repli)|n_benef<0|repli<0){
    return("The number of benefitial alleles (n_benef), the number of replicates (repli) needs to be positive integers or 0")
  }
  if (!is.numeric(min)|!is.numeric(max)|min<0|max<0|min>1|max>1|min>max){
    return("The minimum (min) and maximum (max) selection strength needs to be numeric values between 0 and 1. The minimum value must be smaller or equals to the maximum value.")
  }
  if(fix_sel != round(fix_sel)){
    return("The fixed/discarded number of selected haplotypes must be integers")
  }
  if(length(fix_sel) !=1){
    if(diff(sign(fix_sel)) != 0){
      return("When input a vector as the fixed/discarded number of selected haplotypes, all element in the vector must have the same sign")
    }
  }

  nhap = ncol(haplotype)
  nsnp = nrow(haplotype)
  #one selected snp case
  if (n_benef == 1){
    s_final = list()
    benef_all_final = list()
    num_sel_hap = c()
    for (i in 1:repli){
      s = rep(0, nsnp)
      samp_route = sample(seq(1,nsnp,1),nsnp)
      if (all(fix_sel > 0)){
        #choosing options that have correct number of selected haplotypes
        i = 1
        while (all(fix_sel != sum(haplotype[samp_route[i],]))){
          i = i+1
        }
        benef_all = samp_route[i]
      } else if (all(fix_sel < 0)){
        i = 1
        while (any(c(-fix_sel,0,nhap) == sum(haplotype[samp_route[i],]))){
          i = i+1
        }
        benef_all = samp_route[i]
      } else {
        #excluding options where no haplotypes are selected or all haplotypes are selected
        i = 1
        while (any(c(0,nhap) == sum(haplotype[samp_route[i],]))){
          i = i+1
        }
        benef_all = samp_route[i]
      }
      num_sel_hap = sum(haplotype[benef_all,])
      s[benef_all] = runif(n_benef,min,max)
      s_final = append(s_final,list(s))
      benef_all_final = append(benef_all_final, list(benef_all))
    }

    if (diff_sel_str == FALSE){
      # exact replicate
      s_final = rep(list(s_final[[1]]), repli)
      benef_all_final = rep(list(benef_all_final[[1]]), repli)
    } else {
      s_int = rep(list(rep(0,nsnp)), repli)
      for (i in 1:repli){
        s_int[[i]][benef_all_final[[1]]] = s_final[[i]][benef_all_final[[i]]]
      }
      s_final = s_int
      benef_all_final = rep(list(benef_all_final[[1]]), repli)
    }
  }
  if (n_benef == 0){
    benef_all_final = NULL
    s_final = NULL
    num_sel_hap = 0
  }

  return(list(s_final,benef_all_final,num_sel_hap))
}
