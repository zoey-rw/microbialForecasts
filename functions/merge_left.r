# Merge left function, for NEON data in particular
#
# Merges two tables, keeping all columns/rows from the left one 
# and only the unique columns/overlapping rows from the right one
merge_left <- function(left, right, by.col = NULL, all.x = FALSE) {
  
  if (is.null(by.col)) {  # Choose merge column
    common_col <- intersect(colnames(left), colnames(right))
    merge_cols <- intersect(common_col, c("dnaSampleID", "sampleID", "geneticSampleID", "deprecatedVialID")) # preferred order
    by.col <- merge_cols[[1]] 
    if (is.null(by.col)) {
      cat("No merge column found, please specify using the 'by.col' argument.")
      return()
    }
  }
  
  cat(paste0(by.col, " used as merge column."))
  right_to_merge <- right[,!colnames(right) %in% colnames(left)]
  cat(c("\nDropped the non-unique columns: ", colnames(right)[colnames(right) %in% colnames(left)]))
  right_to_merge[,by.col] <- right[,by.col]
  merged <- merge(left, right_to_merge, all.x = all.x)
  cat(c("\nStarting dimensions for left table: ", dim(left), 
             "\nStarting dimensions for right table: ", dim(right),
             "\nMerged table dimensions:", dim(merged)))
  return(merged)
}
