#' Compute the mean squared error between pseudotime estimates and true timepoints.
#'
#' @param pseudotime An n-vector of pseudotime estimates.
#' @param truetime An n-vector of true timepoint values.
#' @return The MSE between pseudotime and truetime.
#' @export
mse <- function(pseudotime, truetime)
{
  if (cor(pseudotime, truetime, method="kendall") < 0) pseudotime <- rev(pseudotime)
  sum((pseudotime - truetime)^2) / length(truetime)
}


#' Compute the mean absolute error between pseudotime estimates and true timepoints.
#'
#' @param pseudotime An n-vector of pseudotime estimates.
#' @param truetime An n-vector of true timepoint values.
#' @return The MAE between pseudotime and truetime.
#' @export

mae <- function(pseudotime, truetime)
{
  if (cor(pseudotime, truetime, method="kendall") < 0) pseudotime <- rev(pseudotime)
  sum(abs(pseudotime - truetime)) / length(truetime)

}
