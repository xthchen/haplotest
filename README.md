# haplotest
This R package is used for haplotype based selection testing in the context of evolve and resequence. The main function haplotest() takes as input temporal haplotype or allele frequencies, and gives p-value of selection testing as output. Haplotype combination methods can be used if the input contains many haplotypes. Post-hoc test in the form of testing for the number of selected haplotypes and pair-wise tests are also available in this package.

In order to run this package, both the [omnibus](https://github.com/ThomasTaus/omnibus) and [harmonicmeanp](https://github.com/danny-wilson/harmonicmeanp) R packages are needed. 
## Example of haplotest usage
_p_value = haplotest(haplotype_frequency = haplo_freq_mat, delta = 10, Ne = rep(1000,3), repli = 3, test_type = "haplotype", test_method = "base", p_combine_method = "omnibus", freq_est = "known")_

Here haplo_freq_mat is the haplotype frequency matrix, with frequencies available at every 10 generations. There are 3 replicate populations, each with a effective population size of 1000. The haplotype frequencies are assumed to be known, with no noise coming from estimations. If haplotype frequencies are all estimated with some sample size, a sample size vector is needed to be provided to the argument "coverage_matrix", and the argument "freq_est" needs to be set to "est".

SNP based test is also available by changing the "haplotype_frequency" argument to "allele_frequency", and setting test_type = "allele".

For haplotype frequency combination methods in case of large haplotypes, both haplotype and allele frequency needs to be provided, as well as the haplotype structural matrix as argument "haplotype". Two combination methods are available, test_method = "hap_filt" and test_method = "hapblock".

##Example of post hoc tests
Using the same notations, the test for the number of selected haplotypes can be tested using:

_hap_num = hapnumtest(haplo_freq_mat, p_combine_method = "omnibus", deltat = 10, Ne = 1000, repli = 3)_

The pairwise test for fitness difference can be tested using:

_pair_pval = pair_test(haplo_freq_mat, Ne = rep(1000,3), repli = 3, deltat = 10)_

For further information and examples refer to the manual typing ?haplotest in R.
