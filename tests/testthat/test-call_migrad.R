test_that("migrad with the built-in C++ function", {
    x = list(label=c("A","B","C"),
             est=c(10.,10.,.01),
             low=c(0,0,0),
             upp=c(100,20,.1))

    result = call_migrad(x,"test_fn")
    
    expect_equal(result$fmin,-2.579809e+04)
    expect_equal(result$est, c(2.093722e+01, 4.424903e-01, 4.627753e-02), tolerance=1e-5)
    expect_equal(result$label, x$label)
    expect_equal(result$status, 0)
    expect_equal(result$nfcn, 164)
})

## FIXME:TODO: test for "fixed" variables (logicals? integers?)
## FIXME:TODO: test for inconsistent array lengths
## FIXME:TODO: test for missing/inconsistent-array-length labels

test_that("migrad with R function", {
    f = function(x) { (x[1]-1)^2 + (x[2]-2)^2 + 3 } # min=3 at (1,2)
    
    x = list(label=c("x","y"),
             est=c(0., 0.),
             low=c(-1,-1),
             upp=c(5,5))

    result = call_migrad(x,f)
    
    expect_equal(result$fmin, 3)
    expect_equal(result$est, c(1, 2), tolerance=1e-5)
    expect_equal(result$label, x$label)
    expect_equal(result$status, 0)
    expect_equal(result$nfcn, 55)
})

test_that("migrad with R function, fixed vars", {
    f = function(x) { (x[1]-1)^2 + (x[2]-2)^2 + x[3]^2 + 3 } # min=3 at (1,2,0); min=5.25 at (1, fixed(0.5), 0)
    
    x = list(label=c("x","y", "z"),
             fixed=c(FALSE, TRUE, FALSE),
             est=c(0.5, 0.5, 0.5),
             low=c(-1,-1,-1),
             upp=c(5,5,5))

    result = call_migrad(x,f)
    
    expect_equal(result$fmin, 5.25)
    expect_equal(result$est, c(1.0, 0.5, 0.0), tolerance=1e-5)
    expect_equal(result$label, x$label)
    expect_equal(result$status, 0)
})

