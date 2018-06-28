#' Create balance datasets from sets of trees
#' 
#' This function turns the input number into a binary sequence.
#' @param path_to_data A string, path to where the RDS files containing the trees are
#' @param path_to_output A string, path to where the RDS files containing the balance data are (to be saved)
#' @param nLim The number of replicate trees to look for per parameter set
#' @return As many zeros as files in \code{path_to_data} if everything went fine
#' @author Raphael Scherrer
#' @export

make_balance_files <- function(path_to_data, path_to_output, nLim) {
  
  # Load packages
  library(ape)
  library(TotalCopheneticIndex)
  library(balance)
  library(pbapply)
  
  # Read through all tree files
  treeFiles <- list.files(path_to_data)
  
  # Read through all balance files
  balanceFiles <- list.files(path_to_output)
  
  # For each file...
  pbapply::pbSapply(treeFiles, function(curr.file) {
    
    # Is there a balance file with that name?
    outputfile <- gsub("simul_", "balance_", curr.file)
    i <- grep(outputfile, balanceFiles)
    
    # If yes, load it
    if(length(i) > 0) {
      
      balance <- readRDS(paste0(path_to_output, "/", balanceFiles[i]))
      
      # Is the balance file full? If yes, move to the next file
      if(nrow(balance) >= nLim) return(0)
      
    } else {
      
      # Otherwise, create a balance table
      balance <- data.frame()
      
    }
    
    # Load the trees
    trees <- readRDS(paste0(path_to_data, "/", curr.file))
    
    # Set the row to start with
    j <- nrow(balance) + 1
    a <- 1
    
    # Initialize vectors
    nTrees <- tciTrees <- mintciTrees <- maxtciTrees <- ntciTrees <- NULL
    nStrees <- tciStrees <- mintciStrees <- maxtciStrees <- ntciStrees <- NULL
    
    # For each tree from j, record the number of tips and the balance of both the good and the incipient species tree
    while(j <= nLim) {
      
      # Extract trees
      curr.tree <- trees$trees[[j]]
      curr.stree <- trees$strees[[j]]
      
      # Count numbers of tips
      nTrees[a] <- ape::Ntip(curr.tree)
      nStrees[a] <- ape::Ntip(curr.stree)
      
      # Calculate imbalance (total cophenetic index)
      tciTrees[a] <- TotalCopheneticIndex::tci(curr.tree)
      tciStrees[a] <- TotalCopheneticIndex::tci(curr.stree)
      
      # Calculate minimum and maximum TCI for each tree
      mintciTrees[a] <- balance::min_tci(nTrees[a])
      mintciStrees[a] <- balance::min_tci(nStrees[a])
      maxtciTrees[a] <- balance::max_tci(nTrees[a])
      maxtciStrees[a] <- balance::max_tci(nStrees[a])
      
      # Calculate normalized imbalance (corrected for tree size)
      ntciTrees[a] <- (tciTrees[a] - mintciTrees[a]) / (maxtciTrees[a] - mintciTrees[a])
      ntciStrees[a] <- (tciStrees[a] - mintciStrees[a]) / (maxtciStrees[a] - mintciStrees[a])
      
      # Update j and a
      j <- j + 1
      a <- a + 1
      
    }
    
    # Append the results to the balance data frame
    newBalance <- data.frame(nTrees, tciTrees, mintciTrees, maxtciTrees, ntciTrees, nStrees, tciStrees, mintciStrees, maxtciStrees, ntciStrees)
    balance <- rbind(balance, newBalance)
    
    # Save
    saveRDS(balance, paste0(path_to_output, "/", outputfile))
    
    return(0)
    
  })
  
}


