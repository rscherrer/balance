#' Calculate maximum total cophenetic index
#' 
#' This function returns the maximum TCI value (Mir et al. 2013) for a binary tree of size \code{n}.
#' @param n The number of tips in the tree
#' @return An integer, the maximum value of the TCI
#' @keywords Total Cophenetic Index, tree imbalance
#' @export
#' @example 
#' max_tci(5) #should give 10
#' @references Mir, A., Rosselló, F., & Rotger, L. a. (2013). A new balance index for phylogenetic trees. Mathematical Biosciences, 241(1), 125–136. https://doi.org/10.1016/j.mbs.2012.10.005

max_tci <- function(n) {
  
  if(!is.integer(n) | n <= 1) stop("n must be an integer > 1")
  
  print(1)
  
  # The max TCI is the number of combinations of n elements taken 3 at a time (number of possible triplets among n elements)
  return(choose(n, 3))
  
}