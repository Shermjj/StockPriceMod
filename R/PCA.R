#' PCA Function
#'
#' @param data_matrix
#' @param scale PCA scaling
#' @param center PCA centering
#' @return Result from PCA
#' @export
pca <- function(data_matrix, variance_threshold = 0.9){
  pca_result <- prcomp(data_matrix, scale = FALSE, center= TRUE)

  explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cumulative_variance <- cumsum(explained_variance)

  # Find the number of components that explain at least the threshold variance
  num_components <- which(cumulative_variance >= variance_threshold)[1]

  # Return the principal components
  return(list(pcs = pca_result$x[, 1:num_components],
              variance_explained = cumulative_variance[num_components]))
}
