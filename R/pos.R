#' Pseudotime-temporal ordering score
#'
#' Pseudotime-temporal ordering (POS) is a metric introduced in the TSCAN paper, used to compare the estimated ordering of cells to the true ordering.
#'
#' @param pseudotime An n-vector of pseudotime estimates
#' @param truetime An n-vector of true timepoint values
#' @return The positional ordering score
#' @export
pos <- function(pseudotime, truetime)
{
    pseudotime_order <- order(pseudotime)
    truetime_order <- order(truetime)
  if (cor(pseudotime_order, truetime_order, method="kendall") < 0)
    pseudotime_order <- rev(pseudotime_order)

  n <- length(truetime_order)
  gscores <- c()
  for (i in 1:(n - 1))
  {
    for (j in (i+1):n)
    {
      lower_cell <- pseudotime_order[i]
      upper_cell <- pseudotime_order[j]
      gscores <- c(gscores, truetime_order[upper_cell] - truetime_order[lower_cell])
    }
  }
  return(sum(gscores))
}

