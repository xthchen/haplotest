% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MIG.R
\name{MIG}
\alias{MIG}
\title{MIG algorithm for computation of haplotype blocks}
\usage{
MIG(
  haplotype,
  frequency_matrix,
  allele_frequency,
  LDU = 0.98,
  LDL = 0.7,
  EHR = 0.9,
  d = 0.95
)
}
\arguments{
\item{haplotype}{Binary numeric matrix containing information of haplotype structure. The columns corresponds to the haplotypes and the rows to the SNPs.}

\item{frequency_matrix}{Numeric matrix, with each i,j element being the frequency of haplotype i at sequenced time point j.}

\item{allele_frequency}{Numeric matrix, with each i,j element being the frequency of SNP i at sequenced time point j.}

\item{LDU}{Numeric, parameter compared to upper 90 \% confidence interval for normalised LD variance, chosen in line with Gabriel, S. B. et al. (2002)}

\item{LDL}{Numeric, parameter compared to lower 90 \% confidence interval for normalised LD variance, chosen in line with Gabriel, S. B. et al. (2002)}

\item{EHR}{Numeric, parameter for showing strong evidence of historical recombination, chosen in line with Gabriel, S. B. et al. (2002)}

\item{d}{Numeric, parameter regarding haplotype block, chosen in line with Gabriel, S. B. et al. (2002)}
}
\value{
Numeric matrix of filtered haplotype block candidates, with the first column being the starting SNP position of the candidate, and second column ending SNP position of the candidate. Each row will corresponds to a single candidate.
}
\description{
Algorithm for computation of haplotype blocks, from (Taliun, D. et al. 2014).
}
\references{
Gabriel, S. B., Schaffner, S. F., Nguyen, H., Moore, J. M., Roy, J., Blumenstiel, B., Higgins, J., DeFelice, M., Lochner, A., Faggart, M., Liu-Cordero, S. N., Rotimi, C., Adeyemo, A., Cooper, R., Ward, R., Lander, E. S., Daly, M. J., Altshuler, D., (2002), The structure of haplotype blocks in the human genome, Science 296(5576): 2225–2229.

Taliun, D., Gamper, J., Pattaro, C., (2014), Efficient haplotype block recognition of very long and dense genetic sequences, BMC bioinformatics 15: 10.
}
\seealso{
\code{\link[=haplotest]{haplotest()}}
}
