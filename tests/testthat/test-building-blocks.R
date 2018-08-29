context('Default parameter builders')

test_that('create.kernel.params() works as expected',{
  kp <- create.kernel.params()
  expect_length(kp, 6)
  expect_setequal(names(kp), c("U", "L", "n.bins", "y", "h", "b"))
})

test_that('create.vital.params() works as expected',{
  kp <- mock.vital.params()
  expect_length(kp, 20)
})

context('Vital rate functions working')

test_that('survival kernel runs with defaults',{
  S <- test.vital.fun('S')
  expect(all(diag(S)>0))
  delta <- row(S)-col(S) # creating matrix to check off-diag are 0
  expect(all(S[delta>0]==0))
})

test_that('reproduction kernel runs with defaults',{
  R <- test.vital.fun('R')
  delta <- row(R)-col(R) # creating matrix to check off-diag are 0
  expect(all(R[delta>0]==0))
})

test_that('movement kernel runs with defaults',{
  expect_error(test.vital.fun('G'), NA)
})