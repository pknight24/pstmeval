#' Flips the pseudotime estimates over their midpoint
#'
#' @param pseudotime A vector of pseudotime estimates
#' @return A vector of flipped pseudotime estimates
#' @export
flip <- function(pseudotime)
{
    midpoint <- (max(pseudotime) + min(pseudotime)) / 2
    pseudotime - 2 * (pseudotime - midpoint)
}
