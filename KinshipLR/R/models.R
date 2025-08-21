#' Calculate likelihood ratio for Parent-Child relationship
#'
#' @param freq Allele frequency data frame
#' @param n Number of simulations (default: 10000)
#' @return A matrix with likelihood ratios
#' @export
LR_incl_PC <- function(freq, n = 10000) {

  CLR <- matrix(0, nrow = n, ncol = ncol(freq))
  colnames(CLR) <- colnames(freq)

  for (k in 1:n) {
    personA <- randomperson(freq)
    personB <- randomperson(freq)
    personD <- produce_offspring(personA, personB)
    LR <- LR_calculation(personD, personA, freq, r = 1/2)
    CLR[k,] <- LR
    if (k %% 1000 == 0) print(paste("Processing", k, "of", n))
  }
  return(CLR)
}

#' Calculate likelihood ratio for Full Sibling relationship
#'
#' @param freq Allele frequency data frame
#' @param n Number of simulations (default: 10000)
#' @return A matrix with likelihood ratios
#' @export
LR_incl_FS <- function(freq, n = 10000) {

  CLR <- matrix(0, nrow = n, ncol = ncol(freq))
  colnames(CLR) <- colnames(freq)

  for (k in 1:n) {
    personA <- randomperson(freq)
    personB <- randomperson(freq)
    personD <- produce_offspring(personA, personB)
    personE <- produce_offspring(personA, personB)
    LR <- LR_calculation_FS(personD, personE, freq, r = 1/4)
    CLR[k,] <- LR
    if (k %% 1000 == 0) print(paste("Processing", k, "of", n))
  }
  return(CLR)
}

#' Calculate likelihood ratio for Grandparent-Grandchild relationship
#'
#' @param freq Allele frequency data frame
#' @param n Number of simulations (default: 10000)
#' @return A matrix with likelihood ratios
#' @export
LR_incl_GP <- function(freq, n = 10000) {

  CLR <- matrix(0, nrow = n, ncol = ncol(freq))
  colnames(CLR) <- colnames(freq)

  for (k in 1:n) {
    personA <- randomperson(freq)
    personB <- randomperson(freq)
    personC <- randomperson(freq)
    personD <- produce_offspring(personA, personB)
    personI <- produce_offspring(personC, personD)
    LR <- LR_calculation(personA, personI, freq, r = 1/4)
    CLR[k,] <- LR
    if (k %% 1000 == 0) print(paste("Processing", k, "of", n))
  }
  return(CLR)
}

#' Calculate likelihood ratio for Uncle-Nephew relationship
#'
#' @param freq Allele frequency data frame
#' @param n Number of simulations (default: 10000)
#' @return A matrix with likelihood ratios
#' @export
LR_incl_UA <- function(freq, n = 10000) {

  CLR <- matrix(0, nrow = n, ncol = ncol(freq))
  colnames(CLR) <- colnames(freq)

  for (k in 1:n) {
    personA <- randomperson(freq)
    personB <- randomperson(freq)
    personC <- randomperson(freq)
    personD <- produce_offspring(personA, personB)
    personE <- produce_offspring(personA, personB)
    personI <- produce_offspring(personC, personD)
    LR <- LR_calculation(personE, personI, freq, r = 1/4)
    CLR[k,] <- LR
    if (k %% 1000 == 0) print(paste("Processing", k, "of", n))
  }
  return(CLR)
}

#' Calculate likelihood ratio for Half-Sibling relationship
#'
#' @param freq Allele frequency data frame
#' @param n Number of simulations (default: 10000)
#' @return A matrix with likelihood ratios
#' @export
LR_incl_HS <- function(freq, n = 10000) {

  CLR <- matrix(0, nrow = n, ncol = ncol(freq))
  colnames(CLR) <- colnames(freq)

  for (k in 1:n) {
    personG <- randomperson(freq)
    personH <- randomperson(freq)
    personS <- randomperson(freq)
    personO <- produce_offspring(personG, personH)
    personU <- produce_offspring(personH, personS)
    LR <- LR_calculation(personO, personU, freq, r = 1/4)
    CLR[k,] <- LR
    if (k %% 1000 == 0) print(paste("Processing", k, "of", n))
  }
  return(CLR)
}

#' Calculate likelihood ratio for unrelated individuals (r=1/2)
#'
#' @param freq Allele frequency data frame
#' @param n Number of simulations (default: 10000)
#' @return A matrix with likelihood ratios
#' @export
LR_excl_FK <- function(freq, n = 10000) {

  CLR <- matrix(0, nrow = n, ncol = ncol(freq))
  colnames(CLR) <- colnames(freq)

  for (k in 1:n) {
    personX <- randomperson(freq)
    personY <- randomperson(freq)
    LR <- LR_calculation(personX, personY, freq, r = 1/2)
    CLR[k,] <- LR
    if (k %% 1000 == 0) print(paste("Processing", k, "of", n))
  }
  return(CLR)
}

#' Calculate likelihood ratio for unrelated individuals (r=1/4)
#'
#' @param freq Allele frequency data frame
#' @param n Number of simulations (default: 10000)
#' @return A matrix with likelihood ratios
#' @export
LR_excl_SK <- function(freq, n = 10000) {

  CLR <- matrix(0, nrow = n, ncol = ncol(freq))
  colnames(CLR) <- colnames(freq)

  for (k in 1:n) {
    personX <- randomperson(freq)
    personY <- randomperson(freq)
    LR <- LR_calculation(personX, personY, freq, r = 1/4)
    CLR[k,] <- LR
    if (k %% 1000 == 0) print(paste("Processing", k, "of", n))
  }
  return(CLR)
}

#' Calculate likelihood ratio for unrelated individuals using FS method
#'
#' @param freq Allele frequency data frame
#' @param n Number of simulations (default: 10000)
#' @return A matrix with likelihood ratios
#' @export
LR_excl_FS <- function(freq, n = 10000) {

  CLR <- matrix(0, nrow = n, ncol = ncol(freq))
  colnames(CLR) <- colnames(freq)

  for (k in 1:n) {
    personX <- randomperson(freq)
    personY <- randomperson(freq)
    LR <- LR_calculation_FS(personX, personY, freq, r = 1/4)
    CLR[k,] <- LR
    if (k %% 1000 == 0) print(paste("Processing", k, "of", n))
  }
  return(CLR)
}
