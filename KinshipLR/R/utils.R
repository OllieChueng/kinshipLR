#' Create a random person with genotypes
#'
#' @param freq A data frame containing allele frequencies
#' @return A data frame representing a person's genotype
#' @export
randomperson <- function(freq) {
  randMan <- data.frame(matrix(nrow = 2, ncol = length(freq)))
  randMan[,1] <- "R"
  for (i in 2:length(freq)) {
    sample_genotypes <- sample(freq[,1], 2, replace = TRUE, prob = freq[,i])
    randMan[,i] <- sample_genotypes
  }
  return(randMan)
}

#' Create offspring from two parents
#'
#' @param father A data frame representing father's genotype
#' @param mother A data frame representing mother's genotype
#' @return A data frame representing offspring's genotype
#' @export
produce_offspring <- function(father, mother) {
  offspring <- father
  offspring[1,1] <- "F"
  offspring[2,1] <- "M"

  for (i in 1:(ncol(father)-1)) {
    offspring[1,(i+1)] <- father[sample(1:2,1),(i+1)]
    offspring[2,(i+1)] <- mother[sample(1:2,1),(i+1)]
  }
  return(offspring)
}
