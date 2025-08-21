#' Calculate likelihood ratio
#'
#' @param personA First person's genotype
#' @param personB Second person's genotype
#' @param allele_fre Allele frequency data frame
#' @param r Relationship coefficient
#' @return A matrix with likelihood ratios
#' @export
LR_calculation <- function(personA, personB, allele_fre, r) {

  LRdat <- matrix(0, nrow = 1, ncol = ncol(allele_fre))
  colnames(LRdat) <- colnames(allele_fre)

  CLR <- 1
  for (l in 2:ncol(allele_fre)) {
    A <- sort(personA[,l])
    B <- sort(personB[,l])
    G1 = 1/(allele_fre[which(allele_fre[,1]==B[1]),l]) + 1
    G2 = 1/(allele_fre[which(allele_fre[,1]==B[2]),l]) + 1
    m = length(intersect(intersect(A,B),B[1]))
    n = length(intersect(intersect(A,B),B[2]))

    if((A[1]==A[2])||(B[1]==B[2])) {
      LR = r*(G1^m+G2^n) + (1-4*r)
    } else {
      LR = (r/2)*(G1^m+G2^n) + (1-3*r)
    }
    CLR <- CLR * LR
    LRdat[1,l] <- LR
  }
  LRdat[1,1] <- CLR
  return(LRdat)
}

#' Calculate likelihood ratio for Full Sibling
#'
#' @param personA First person's genotype
#' @param personB Second person's genotype
#' @param allele_fre Allele frequency data frame
#' @param r Relationship coefficient
#' @return A matrix with likelihood ratios
#' @export
LR_calculation_FS <- function(personA, personB, allele_fre, r) {

  LRdat <- matrix(0, nrow = 1, ncol = ncol(allele_fre))
  colnames(LRdat) <- colnames(allele_fre)

  CLR <- 1
  for (l in 2:ncol(allele_fre)) {
    A <- sort(personA[,l])
    B <- sort(personB[,l])
    G1 = 1/(allele_fre[which(allele_fre[,1]==B[1]),l]) + 1
    G2 = 1/(allele_fre[which(allele_fre[,1]==B[2]),l]) + 1
    m = length(intersect(intersect(A,B),B[1]))
    n = length(intersect(intersect(A,B),B[2]))

    if((A[1]==A[2])||(B[1]==B[2])) {
      LR = r*(G1^m+G2^n) + (1-4*r)
    } else {
      LR = (r/2)*(G1^m+G2^n) + (1-3*r) - (1/8)
    }
    CLR <- CLR * LR
    LRdat[1,l] <- LR
  }
  LRdat[1,1] <- CLR
  return(LRdat)
}
