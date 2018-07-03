context('Full model test run')

test_that('The full model runs without error',{
  setwd('../../sandbox')
  expect_error(source('DIPM_Model.R'), NA)
})