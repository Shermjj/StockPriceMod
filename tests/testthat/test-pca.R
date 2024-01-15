library(testthat)

# Sample data for testing
data_matrix <- matrix(rnorm(100), ncol = 10)

# Test 1: Check if the function returns the correct number of principal components
test_that("Correct number of principal components are returned", {
  result <- pca(data_matrix, variance_threshold = 0.9)
  expected_num_components <- length(which(cumsum(prcomp(data_matrix, scale = FALSE, center= TRUE)$sdev^2 /
                                                   sum(prcomp(data_matrix, scale = FALSE, center= TRUE)$sdev^2)) < 0.9)) + 1
  expect_equal(ncol(result$pcs), expected_num_components)
})

# Test 2: Check if the function handles different cumulative variance thresholds correctly
test_that("Handles different variance thresholds correctly", {
  for (threshold in seq(0.7, 0.95, by = 0.05)) {
    result <- pca(data_matrix, variance_threshold = threshold)
    cumulative_variance <- sum(result$variance_explained)
    expect_true(cumulative_variance >= threshold)
  }
})

