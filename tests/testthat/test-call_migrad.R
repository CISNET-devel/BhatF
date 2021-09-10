test_that("migrad works", {
    data = matrix(data=c(0,0,0,
                         10.,10.,0.01,
                         10.,10.,0.01,
                         0.,0.,0.,
                         100, 20, 0.1),
                  nrow=3, ncol=5,
                  dimnames=list(c("A","B","C")))
    result = call_migrad(data)
    expect_equal(result$fmin,-2.579809e+04)
    expect_equal(result$est, c(2.093722e+01, 4.424903e-01, 4.627753e-02))
    expect_equal(result$label, dimnames(data)[[1]])
    expect_equal(result$status, 0)
    expect_equal(result$nfcn, 0)
})
