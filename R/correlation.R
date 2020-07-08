#' Compute the correlation between pairwise distances
#'
#' This computes the pairwise distances between cells according to the pseudotime
#' estimates and truetime points, and then returns the Spearman correlation between these
#' distances.
#' It is recommended to normalize the pseudotime and truetime values before calling this function.
#' @param pseudotime An n-vector of pseudotime estimates.
#' @param truetime An n-vector of true time values.
#' @return The spearman correlation between computed distances.
#'
#' @importFrom stats cor dist
#' @importFrom magrittr %>%
#' @export
distance_correlation <- function(pseudotime, truetime)
{
  if (cor(pseudotime, truetime, method="kendall") < 0) pseudotime <- rev(pseudotime)

  ps_distances <- dist(pseudotime) %>% as.vector
  true_distances <- dist(truetime) %>% as.vector
  cor(ps_distances, true_distances, method = "spearman")
}

#' Correlation between gene importances computed by random forests.
#'
#' This metric is introduced in the review paper by Saelens et al. We first choose trajectory "milestones" according to our knowledge of the true trajectory. Then, for each milestone, we compute the geodesic distance according our estimated trajectory from each cell to the milestone, and use the normalized gene expression data as a covariate matrx to predict these distances with a Random Forests model. We then extract the Gini impurity index for each gene (as a measure of gene "importance"). This process is then repeated using the geodesic distances according to the true trajectory. After computing the Gini importances for each gene at all of the milestones, we take the average importance of each gene, and compute the Pearson correlation of the estimated importances and the true importances.
#'
#'
#' High correlation values are supposed to indicate that the estimated trajectory is informed by the same genes as the true trajectory, which we can take as a sort of surrogate for accuracy.
#'
#' @param expression_matrix A p x n normalized gene expression matrix.
#' @param pseudotime An n-vector of pseudotime estimates.
#' @param truetime An n-vector of true timpoints.
#' @param milestones A subset of the true time points, used to compute Gini impurty indices.
#' @param ... Additional parameters passed to \code{ranger}.
#'
#' @return A list containing the gene importances according to pseudotime, true time, and the Pearson correlation between the two.
#' @importFrom ranger ranger
#' @export
rf_importance_correlation <- function(expression_matrix, pseudotime, truetime, milestones, ...)
{
  if (cor(pseudotime, truetime, method="kendall") < 0) pseudotime <- rev(pseudotime)

  ps_importances <- matrix(0, nrow=nrow(expression_matrix), ncol=length(milestones))
  true_importances <- matrix(0, nrow=nrow(expression_matrix), ncol=length(milestones))
  colnames(ps_importances) <- as.character(milestones)
  colnames(true_importances) <- as.character(milestones)
  rownames(ps_importances) <- rownames(expression_matrix)
  rownames(true_importances) <- rownames(expression_matrix)

  for (i in 1:length(milestones))
  {
    ps_dist <- abs(milestones[i] - pseudotime)
    true_dist <- abs(milestones[i] - truetime)
    ps_rf <- ranger(y = ps_dist,x = t(expression_matrix),
                    importance = "impurity", ...)
    true_rf <- ranger(y = true_dist,x = t(expression_matrix),
                    importance = "impurity", ...)
    ps_importances[,i] <- ps_rf$variable.importance
    true_importances[,i] <- true_rf$variable.importance
  }
  correlation <- cor(rowMeans(ps_importances), rowMeans(true_importances))
  list(ps_importances = ps_importances, true_importances = true_importances, correlation = correlation)
}

#' Correlation between LASSO coefficients.
#'
#' Similar to the random forest correlation metric, but rather than Gini impurity indices we compute the coefficients corresponding to each gene with a LASSO regression model.
#'
#' @param expression_matrix A p x n normalized gene expression matrix.
#' @param pseudotime An n-vector of pseudotime estimates.
#' @param truetime An n-vector of true timpoints.
#' @param milestones A subset of the true time points, used to compute LASSO coefficients.
#' @param ... Additional parameters passed to \code{cv.glmnet}.
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
#' @export
lasso_coef_correlation <- function(expression_matrix, pseudotime, truetime, milestones, ...)
{
  if (cor(pseudotime, truetime, method="kendall") < 0) pseudotime <- rev(pseudotime)

  ps_betahat <- matrix(0, nrow=nrow(expression_matrix), ncol=length(milestones))
  true_betahat <- matrix(0, nrow=nrow(expression_matrix), ncol=length(milestones))
  colnames(ps_betahat) <- as.character(milestones)
  colnames(true_betahat) <- as.character(milestones)
  rownames(ps_betahat) <- rownames(expression_matrix)
  rownames(true_betahat) <- rownames(expression_matrix)

  for (i in 1:length(milestones))
  {
    ps_dist <- abs(milestones[i] - pseudotime)
    true_dist <- abs(milestones[i] - truetime)
    ps_cv <- cv.glmnet(y = ps_dist,x = t(expression_matrix),alpha = 1, intercept = FALSE, ...)
    true_cv <- cv.glmnet(y = true_dist,x = t(expression_matrix), alpha = 1, intercept=FALSE, ...)
    ps_betahat[,i] <- coef(ps_cv, s = "lambda.min")[-1,1]
    true_betahat[,i] <- coef(true_cv, s = "lambda.min")[-1,1]
  }
  correlation <- cor(rowMeans(ps_betahat), rowMeans(true_betahat))
  list(ps_betahat = ps_betahat, true_betahat = true_betahat, correlation = correlation)
}
