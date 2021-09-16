test_that("migrad with the built-in function", {
    x = list(label=c("A","B","C"),
             est=c(10.,10.,.01),
             low=c(0,0,0),
             upp=c(100,20,.1))
    f = function(x) { 0 }
    result = call_migrad(x,f)
    
    expect_equal(result$fmin,-2.579809e+04)
    expect_equal(result$est, c(2.093722e+01, 4.424903e-01, 4.627753e-02), tolerance=1e-5)
    expect_equal(result$label, x$label)
    expect_equal(result$status, 0)
    expect_equal(result$nfcn, 164)
})

## FIXME:TODO: test for "fixed" variables (logicals? integers?)
## FIXME:TODO: test for inconsistent array lengths
## FIXME:TODO: test for missing/inconsistent-array-length labels

