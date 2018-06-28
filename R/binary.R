#' Convert to binary
#' 
#' This function turns the input number into a binary sequence.
#' @param x A number
#' @return The binary version of \code{x} as a vector of 0 and 1
#' @author Michael J Crawley
#' @export
#' @example 
#' binary(3) #should give 11
#' @references Crawley, Michael J. The R book. John Wiley & Sons, 2012.

binary <- function(x) {
  i <- 0
  string <- numeric(32)
  while(x > 0) {
    string[32 - i] <- x %% 2
    x <- x %/% 2
    i <- i + 1 
  }
  first <- match(1, string)
  string[first:32] 
}
