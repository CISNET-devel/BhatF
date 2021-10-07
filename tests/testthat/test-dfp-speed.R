test_that("DFP benchmark on R function", {
    if (!getOption("BhatF.run_benchmarks", default=FALSE)) {
        skip("To run benchmarks, set option `BhatF.run_benchmarks` to TRUE")
        return()
    }

    ## nvars = 80
    nvars = 30
    nparams = 5000
    
    xmin = -5
    xmax = 5
    aa = runif(nvars+nparams, min=1, max=10)
    bb = runif(nvars+nparams, min=0.5*xmin, max=0.5*xmax)
    cc = 5
    parms = runif(nparams, min=0, max=1)
    est = runif(nvars, min=0.1*xmin, max=0.1*xmax)

    fn = function(x) {
        for (i in 1:nvars) {
            cc = cc + aa[i]*(x[i]-bb[i])^2
        }
        for (i in 1:nparams) {
            cc = cc + aa[i+nvars]*(parms[i]-bb[i+nvars])^2
        }
        return(cc)
    }
    expected_min = fn(bb[1:nvars])

    x = list(label=paste0("x", 1:nvars),
             est=est,
             low=rep(xmin, nvars),
             upp=rep(xmax, nvars))
    
    #suppressMessages(
        result<-dfp(x, fn)
    #)

    expect_equal(result$fmin, expected_min)
    expect_equal(result$est, bb[1:nvars], tolerance=1e-5)
    expect_equal(result$label, x$label)
    expect_equal(result$status, 0)
})
