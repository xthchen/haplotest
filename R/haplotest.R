#' Testing for selection using temporal haplotype/allele frequencies.
#'
#' Selection can be tested using either haplotype or allele frequencies. For haplotype frequencies, 3 methods are available, basic, haplotype filtering using SNPs, haplotype block based. All tests uses the adapted chi-squared test to correct for variance during testing, dependent on whether the frequencies are estimated. For haplotype filtering and haplotype block based method, both haplotype and allele frequencies, as well as the haplotype structure will be required.
#' @param haplotype Binary numeric matrix containing information of haplotype structure. The columns corresponds to the haplotypes and the rows to the SNPs.
#' @param haplotype_frequency Numeric matrix, with each i,j element being the frequency of haplotype i at sequenced time point j.
#' @param allele_frequency Numeric matrix, with each i,j element being the frequency of SNP i at sequenced time point j.
#' @param coverage_matrix Numeric matrix, required if frequencies are estimated. Sample size matrix of either haplotype/allele frequencies. If test_method == "hap_filt" for haplotype based test, sample size matrix for both haplotype and allele needs to be present. It should  be in the form of list(haplotype_sample_size, allele_sample_size).
#' @param deltat Numeric, number of generations between each pair of time points of interests
#' @param Ne Numeric vector with length as number of replicates, containing information of Ne (effective population size) at each replicated population. If Ne changes over time, take as input a numeric matrix, with the column being the replicate position, row being the Ne at each sequenced time points.
#' @param repli Numeric, specifying the number of replicated populations.
#' @param test_type Factor, can be either "haplotype" or "allele", denoting whether to use haplotype based test or snp based test
#' @param test_method Factor, can be "base", "hap_filt" or "hap_block" for haplotype based, denoting basic haplotype based test, haplotype filtering using snp based test, or haplotype block based test. Can only take argument "base" for snp based.
#' @param p_combine_method Factor, method of pvalue combination, can be "omnibus" (Futschik, A. et al. 2019), "harmonic" (Wilson, D.J. (2019)), "vovk" (Vovk, V. et al. 2018), "bonferroni", "BH" (Benjamini & Hochberg 1995).
#' @param freq_est Factor, whether the frequencies are known or estimated. Input "known" if all frequencies are known, "half" if starting frequencies are known and others estimated, "est" if all are estimated.
#' @param significance Numeric, must be between 0 and 1, level of significance for snp based test if "hap_filt" method is used.
#' @param seed setting seed of the run
#' @references Spitzer, K., Pelizzola, M., Futschik, A., (2019), Modifying the Chi-square and the CMH test for population genetic inference: adapting to over-dispersion, The Annals of Applied Statistics 14(1): 202-220.
#' @references Futschik, A., Taus, T., Zehetmayer, S., (2019), An omnibus test for the global null hypothesis, Statistical Methods in Medical Research 28(8): 2292-2304.
#' @references Wilson, D. J., (2019), The harmonic mean p-value for combining dependent tests, Proceedings of the National Academy of Sciences of the United States of America 116(4): 1195-1200.
#' @references Benjamini, Y., Hochberg, Y., (1995), Controlling the false discovery rate: a practical and powerful approach to multiple testing, Journal of the Royal Statistical Society Series B 57(1): 289-300.
#' @references Vovk, V., Wang., R., (2018), Combining p-values via averaging, Biometrika 107(4): 791–808.
#' @return A numeric vector of p-values from the test after multiple testing corrections.
#'
#'@import omnibus
#'@import harmonicmeanp
#' @examples
#' #We show here examples from simulated data to p-value results of a haplotype based test.
#' #The data simulation part can be ignored if a real data set is used.
#' ##########Known frequency simulated data using basic haplotype based test##########
#'
#' #Simulate haplotype structure with 500 SNPs and 5 haplotypes:
#' hap_struc = Haplotype_sim(nsnp = 500,nhap = 5)
#'
#' #Simulation of selected SNP, with 1 selected SNP private to 1 haplotype, 1 population and 0.05 selective strength:
#' SNP_sel = benef_sim(haplotype = hap_struc, n_benef = 1, min = 0.05, max = 0.05, fix_sel = 1)
#'
#' #Simulation of haplotype trajectory, 60 generations, sequenced at every 10 timepoints, Ne of 500:
#' hap_freq = Frequency_sim(haplotype = hap_struc, t = 6, tdelta = 10, benef_sim = SNP_sel,Ne = 500)[[1]]
#'
#' #Computation of p-value using basic haplotype based test and omnibus p-value combination method:
#' pval = haplotest(haplotype = hap_struc, haplotype_frequency = hap_freq, deltat = 10, Ne = 500, test_type = "haplotype", test_method = "base", p_combine_method = "omnibus", freq_est = "known")
#'
#' ##########Known frequency simulated data using HapSNP##########
#'
#' #In addition to the previous simulations, we need to compute the allele frequencies:
#' allele_freq = hap_struc %*% hap_freq
#'
#' #Computation of p-value using HapSNP and omnibus p-value combination method:
#' pval = haplotest(haplotype = hap_struc, haplotype_frequency = hap_freq, allele_frequency = allele_freq, deltat = 10, Ne = 500, test_type = "haplotype", test_method = "hap_filt", p_combine_method = "omnibus", freq_est = "known")
#'
#' ##########Estimated frequency simulated data using basic haplotype based test##########
#'
#' #Simulation of noisy haplotype frequencies with mean samplesize 80:
#' hap_freq_est = hap_err_sim(hap_freq, 80)
#'
#' #Computation of p-value using basic haplotype based test and omnibus p-value combination method:
#' pval = haplotest(haplotype = hap_struc, haplotype_frequency = hap_freq_est[[1]], coverage_matrix = hap_freq_est[[2]], deltat = 10, Ne = 500, test_type = "haplotype", test_method = "base", p_combine_method = "omnibus", freq_est = "est")
#'
#' @export



haplotest = function(haplotype = NULL, haplotype_frequency = NULL, allele_frequency = NULL, coverage_matrix = NULL, deltat = 10, Ne = 1000, repli = 1, test_type = "haplotype",
                         test_method = "base", p_combine_method = "omnibus", freq_est = "known",significance = 0.05, seed = 2021){

  if (p_combine_method != "omnibus" & p_combine_method != "harmonic" & p_combine_method != "vovk" & p_combine_method != "bonferroni" & p_combine_method != "BH"){
    return ("The given multiple testing method is not available for this function")
  }else if (test_method == "base" & freq_est == "half"){
    if (test_type == "haplotype"){
      if (ncol(haplotype_frequency) != ncol(coverage_matrix)){
        return("The haplotype frequency matrix have different dimensions then the sample size matrix")
      }else if (nrow(haplotype_frequency) != nrow(coverage_matrix)){
        return("The haplotype frequency matrix have different dimensions then the sample size matrix")
      }
    }
  }else if (test_method == "base" & freq_est == "est"){
    if (test_type == "haplotype"){
      if (ncol(haplotype_frequency) != ncol(coverage_matrix)){
        return("The haplotype frequency matrix have different dimensions then the sample size matrix")
      }else if (nrow(haplotype_frequency) != nrow(coverage_matrix)){
        return("The haplotype frequency matrix have different dimensions then the sample size matrix")
      }
    }
  }else if (test_type == "haplotype"){
    if (is.null(haplotype_frequency)){
      return("The frequency matrix of the chosen test_type is missing")
    } else if (test_method != "base" & test_method != "hap_filt" & test_method != "hap_block"){
      return("The chosen test_method is not available")
    }
  }else if (test_type == "allele"){
    if (is.null(allele_frequency)){
      return("The frequency matrix of the chosen test_type is missing")
    }else if (test_method != "base"){
      return ("The chosen test_method is not available")
    }
  }else if (test_type != "haplotype" & test_type != "allele"){
    return("The chosen test_type is not available")
  }

  set.seed(seed)
  if (!is.null(allele_frequency) & !is.null(haplotype)){
    allele_frequency = allele_frequency[rowSums(allele_frequency)>0,]
    haplotype = haplotype[rowSums(haplotype) > 0,]
    allele_frequency = allele_frequency[rowSums(allele_frequency)<ncol(allele_frequency),]
    haplotype = haplotype[rowSums(haplotype) < ncol(haplotype),]
  }
  if (!is.null(allele_frequency)){
    #allele_frequency = allele_frequency[rowSums(haplotype)>0,]
    allele_frequency = allele_frequency[rowSums(allele_frequency)>0,]

    #allele_frequency = allele_frequency[rowSums(haplotype)<ncol(haplotype),]
    allele_frequency = allele_frequency[rowSums(allele_frequency)<ncol(allele_frequency),]
  }
  if (test_method == "base"){
    if (test_type == "haplotype"){
      freq = haplotype_frequency
    } else {
      freq = allele_frequency
    }
    t = (ncol(freq)/repli)-1
    if (freq_est == "known"){
      pval = adapted.cmh.test.true(freq = freq, repli = repli, Ne = Ne, gen = seq(0,t*deltat,deltat))
    }else if (freq_est == "half"){
      pval = adapted.cmh.test.str.true(freq = freq, repli = repli, Ne = Ne, cov = coverage_matrix,gen = seq(0,t*deltat,deltat))
    }else if (freq_est == "est"){
      pval = adapted.cmh.test.est(freq = freq, repli = repli, Ne = Ne, cov = coverage_matrix,gen = seq(0,t*deltat,deltat))
    }
  } else if (test_method == "hap_filt"){
    t = ncol(allele_frequency)/repli-1
    if (freq_est == "known"){
      pval_allele = adapted.cmh.test.true(freq = allele_frequency, Ne = Ne, gen = seq(0,t*deltat,deltat), repli = repli)
    }else if (freq_est == "half"){
      pval_allele = adapted.cmh.test.str.true(freq = allele_frequency, Ne = Ne, cov = coverage_matrix[[2]], gen = seq(0,t*deltat,deltat), repli = repli)
    }else if (freq_est == "est"){
      pval_allele = adapted.cmh.test.est(freq = allele_frequency, Ne = Ne, cov = coverage_matrix[[2]], gen = seq(0,t*deltat,deltat), repli = repli)
    }

    allele_pos = which(as.numeric(p.adjust(pval_allele, method = "BH")) < significance)

    if (length(allele_pos) > 0){
      if (length(allele_pos) == 1){
        #if only one snp satisfy the threshold, highly unlikely
        red_haplotype = t(as.matrix(haplotype[allele_pos,]))
      }else{
        red_haplotype = haplotype[allele_pos, ]
        #matrix with only SNPs satisfying threshold
      }
      #combining haplotype frequencies
      red_haplotype_frequency = hap_snp_red(red_haplotype, haplotype_frequency)

      red_coverage = coverage_matrix[[1]][1:nrow(red_haplotype_frequency),]
      if (freq_est == "known"){
        pval = adapted.cmh.test.true(freq = red_haplotype_frequency, Ne = Ne, gen = seq(0,t*deltat, deltat), repli = repli)
      }else if (freq_est == "half"){
        pval = adapted.cmh.test.str.true(freq = red_haplotype_frequency, Ne = Ne, cov = red_coverage, gen = seq(0,t*deltat, deltat), repli = repli)
      }else if (freq_est == "est"){
        pval = adapted.cmh.test.est(freq = red_haplotype_frequency, Ne = Ne, cov = red_coverage, gen = seq(0,t*deltat, deltat), repli = repli)
      }

    } else {
      if (freq_est == "known"){
        pval = adapted.cmh.test.true(freq = haplotype_frequency, Ne = Ne, gen = seq(0,t*deltat,deltat), repli = repli)
      } else if (freq_est == "half"){
        pval = adapted.cmh.test.str.true(freq = haplotype_frequency, Ne = Ne, cov = coverage_matrix[[1]], gen = seq(0,t*deltat,deltat), repli = repli)
      } else if (freq_est == "est"){
        pval = adapted.cmh.test.est(freq = haplotype_frequency, Ne = Ne, cov = coverage_matrix[[1]], gen = seq(0,t*deltat,deltat), repli = repli)
      }
      warning("SNP based test did not find any significance, will use basic haplotype based test instead")
    }
  } else {
    t = ncol(allele_frequency)-1
    hap_block_cand = MIG(haplotype, haplotype_frequency, allele_frequency, LDU = 0.98, LDL = 0.7, EHR = 0.9, d = 0.95)
    #compute all potential haplotype candidates
    hap_block_filt = hap_block(hap_block_cand)
    #filtered by favouring larger haplotype block first
    hap_block_filt = hap_block_filt[((hap_block_filt[,2] - hap_block_filt[,1]) > 1),]
    #remove haplotype blocks that are only 2 snps in length
    pval = c()
    for (num in 1:nrow(hap_block_filt)){
      red_haplotype = haplotype[seq(hap_block_filt[num,1],hap_block_filt[num,2],1),]
      red_haplotype_frequency = hap_snp_red(red_haplotype, haplotype_frequency)
      red_coverage = coverage_matrix[1:nrow(red_haplotype_frequency),]
      if (freq_est == "known"){
        red_pvalue = adapted.cmh.test.true(freq = red_haplotype_frequency, Ne = Ne, gen = seq(0,t*deltat, deltat), repli = repli)
      }else if (freq_est == "half"){
        red_pvalue = adapted.cmh.test.str.true(freq = red_haplotype_frequency, Ne = Ne, cov = red_coverage, gen = seq(0,t*deltat, deltat), repli = repli)
      }else if (freq_est == "est"){
        red_pvalue = adapted.cmh.test.est(freq = red_haplotype_frequency, Ne = Ne, cov = red_coverage, gen = seq(0,t*deltat, deltat), repli = repli)
      }
      pval = c(pval, p.hmp(red_pvalue+10^(-7), L = length(red_pvalue)))

    }
  }
  if (p_combine_method == "omnibus"){
    comb_pval = omnibus.test(pval, method = "log.p", N.sim = 1000)
  } else if (p_combine_method == "harmonic"){
    comb_pval = p.hmp(pval+10^(-7), L = length(pval))
  } else if (p_combine_method == "bonferroni"){
    comb_pval = min(as.numeric(p.adjust(pval, method = "bonferroni")))
  } else if (p_combine_method == "BH"){
    comb_pval = min(as.numeric(p.adjust(pval, method = "BH")))
  } else {
    comb_pval = merging_function(pval, -1)
    #maybe allow other parameter input for vovk method at start
  }
  return(comb_pval)
}
