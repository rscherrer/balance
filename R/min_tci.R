#' Calculate minimum total cophenetic index
#' 
#' This function returns the minimum TCI value (Mir et al. 2013) for a binary tree of size \code{n}.
#' @param n The number of tips in the tree
#' @return An integer, the minimum value of the TCI
#' @keywords Total Cophenetic Index, tree imbalance
#' @export
#' @example 
#' min_tci(5) #should give 5
#' @references Mir, A., Rosselló, F., & Rotger, L. a. (2013). A new balance index for phylogenetic trees. Mathematical Biosciences, 241(1), 125–136. https://doi.org/10.1016/j.mbs.2012.10.005

min_tci <- function(n) {
  
  if(!is.integer(n) | n <= 1) stop("n must be an integer > 1")
  
  # For each integer k from 0 to n-1, calculate the highest power of 2 that divides k! (factorial k) (sequence A011371 in Sloane's Encyclopedia of Integer Sequences)
  a <- sapply(0:(n-1), function(k) {

    # First turn k into binary
    vecBinary <- binary(k)
    
    # Then sum all ones
    sumOnes <- sum(vecBinary)
    
    # Substract this number to k
    return(k - sumOnes)
    
  })
  
  # Sum these numbers over all k values
  return(sum(a))
  
}