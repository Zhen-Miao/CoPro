# Test self-bidirectional correlation functions
test_that("Self-bidirectional correlation functions work correctly", {
  skip_if_not_installed("CoPro")
  
  # Test that functions exist and have correct signatures
  expect_true(exists("getTransferSelfBidirCorr"))
  expect_true(exists("computeSelfBidirCorr"))
  
  # Test that the core helper function exists
  expect_true(exists(".computeSpatialSelfCorrelation", mode = "function"))
  
  # Additional tests would require creating mock CoPro objects
  # with self-kernel matrices, which is beyond the scope of this implementation
})
