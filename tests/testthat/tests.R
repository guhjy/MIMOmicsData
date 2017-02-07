library(MIMOmicsData)
context("MIMOmicsData main function checking")

test_that("Normal input goes without error", {
  expect_error(GenData(), NA)
})

# test_that("Examples run", {
#   expect_error(example("o2m",package = "O2PLS"), NA)
#   expect_error(example("crossval_o2m",package = "O2PLS"), NA)
#   expect_error(example("crossval_o2m_adjR2",package = "O2PLS"), NA)
#   expect_error(example("orth",package = "O2PLS"), NA)
#   expect_error(example("summary.o2m",package = "O2PLS"), NA)
#   expect_error(example("loadings.o2m",package = "O2PLS"), NA)
# })
#
# test_that("Errors in o2m are thrown", {
#   expect_error(o2m(diag(3),diag(4),1,0,0),"rows")
# })
#
# test_that("size, ratios and names are correct", {
#   expect_equal(o2m(diag(4),diag(4),1,1,0)$R2X,      0.5)
# })
#
# test_that("S3 Methods are working OK", {
#   fit = o2m(data.frame(a=1:10,b=2:11,c=3:12),data.frame(d=1:10,e=2:11,f=3:12),1,0,0)
#   expect_error(print(fit), NA)
# })
