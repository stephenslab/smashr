test_that("smash correctly dispatches to smash.pois",{
  expect_equal(smash(1:8),smash.poiss(1:8))
})