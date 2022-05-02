#' Simulation of haploid haplotype frequencies over time.
#' 
#' Simulation of haplotype frequencies over time using a multinomial distribution. Simulations can be run with or without replicated populations, and for constant or changing selection over time. Frequencies are recorded at time points provided in the input data.
#' @param haplotype Binary numeric matrix containing information of haplotype structure. The columns corresponds to the haplotypes and the rows to the SNPs.
#' @param t Numeric, number of time points of interest, excluding the initial time point
#' @param deltat Numeric, number of generations between each pair of time points of interests
#' @param s Numeric vector of selection coefficients. s = 0 corresponds to neutral simulations. Can use first element of list from output of benef_sim() as input.
#' @param benef_all Numeric vector of position of beneficial alleles. Can use second element of list from output of benef_sim() as input.
#' @param rand_start Logical. rand_start = TRUE corresponds to random sampling of the starting haplotype frequency from a uniform distribution. If rand_start = FALSE the starting haplotype frequency is set to equal for all haplotypes.
#' @param Ne Positive integer, the effective population size.
#' @param sel_dec Logical. sel_dec = TRUE corresponds to decreasing selective strength (fitness of selected haplotype) over time
#' @param dec_ne Integer, strength by which Ne decreases over time. dec_ne = number of Ne decreased per generation. Contant Ne if equals to 0.
#' @param dec_strength Numeric, strength by which selection decreases if sel_dec=TRUE. If sel_dec = FALSE this parameter is ignored.
#' @param repli Numeric, specifying the number of replicated populations.
#' @return Returns a list. The first item is the haplotype frequency matrix with row being the haplotypes and columns the chosen timepoints. If replicates are present, the column will be the time points of first replicate, followed by the timepoints of second replicates etc. Second item is the fitness of haplotypes.
#' @export


Frequency_sim = function(haplotype, t = 6, tdelta = 10, benef_sim, rand_start = FALSE,
                         Ne = 500, sel_dec = FALSE, dec_strength = 0.7, dec_ne = 0, repli = 1){
  nhap = ncol(haplotype)
  s = benef_sim[[1]]
  benef_all = benef_sim[[2]]
  freq_matrix = matrix(NA, ncol= (t+1)*repli, nrow = nhap)
  #frequency matrix initialisation
  if (rand_start == TRUE){
    time_0 = runif(nhap,0,1)
    time_0 = time_0/sum(time_0)
  }else{
    time_0 = rep(1/nhap, nhap)
  }
  #initialisation of haplotype frequencies at time 0
  for (i in 1:repli){
    freq_matrix[,1 + (i-1)*(t+1)] = time_0
  }
  #record frequency at time 0 into matrix
  timepoints = seq(0,t*tdelta,tdelta)
  #timepoints of interests
  total_fitness = c()
  for (k in 1:repli){
    fitness = rep(1,nhap)
    #initialisation for fitness, set all to 1
    if (length(s) == 0){
      fitness = fitness
    }
    else {
      #apply selection if there are more than one selection required
      for (i in 1:nhap){
        #loop through all haplotypes
        if (sum(haplotype[benef_all[[k]], i]) > 0){
          #if not all benef_all rows at column i are 0
          fitness[i] = fitness[i]+sum(s[[k]][benef_all[[k]]][as.logical(haplotype[benef_all[[k]],i])])
          #fitness of haplotype i = sum(wA[benef_all]) only if hp_str[benef_all,i] = 1
          #else = 0
          #basically a sum of selections for each haplotype
        }
      }
    }
    freq = time_0
    #frequency initialisation
    ti = 1
    #current generation
    while(ti <= t*tdelta){
      if (length(s) > 0){
        freq = freq*fitness
      }
      # apply drift
      if (!is.na(Ne))
        #as long as effective population size is present
        freq = rmultinom(1, Ne, freq)/Ne
      #frequency follows a multinomial distribution
      
      if (ti %in% timepoints){
        freq_matrix[,(ti/tdelta+1)+(t+1)*(k-1)] = freq
        if (sel_dec == TRUE){
          fitness = 1 + (fitness - 1)*dec_strength
        }
      }
      ti = ti + 1
      Ne = Ne - dec_ne
    }
    total_fitness = rbind(total_fitness, fitness)
  }
  
  return(list(freq_matrix, total_fitness))
}
