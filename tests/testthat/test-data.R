test_that("example_data loads correctly", {
  data(example_data, package = "SafeBRMS")
  
  expect_s3_class(example_data, "data.frame")
  expect_equal(nrow(example_data), 100)
  expect_equal(ncol(example_data), 4)
  expect_named(example_data, c("y", "x1", "x2", "x3"))
})

test_that("example_data has expected properties", {
  data(example_data, package = "SafeBRMS")
  
  # All columns should be numeric
  expect_true(all(sapply(example_data, is.numeric)))
  
  # No missing values
  expect_false(any(is.na(example_data)))
  
  # Reasonable ranges (should be roughly standard normal)
  expect_true(all(abs(example_data$x1) < 5))
  expect_true(all(abs(example_data$x2) < 5))
  expect_true(all(abs(example_data$x3) < 5))
})
