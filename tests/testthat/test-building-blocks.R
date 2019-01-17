context('Default parameter builders')

test_that('create.kernel.params() works as expected',{
  kp <- create.kernel.params(0, 10, 10)
  expect_length(kp, 6)
  expect_setequal(names(kp), c("U", "L", "n.bins", "y", "h", "b"))
  
  expect_equal(kp$b, 0:10)
  expect_equal(kp$h, 1)
  expect_equal(kp$y, 0:9+0.5)
  
  kp <- create.kernel.params(-12, 20.4, 57)
  expect_length(kp, 6)
  expect_length(kp$y, 57)
  expect_length(kp$b, 57+1)
  expect_equal(kp$y[1], -12 + (20.4-(-12))/(2*57))
  expect_equal(kp$h, (20.4-(-12))/(57))
})

