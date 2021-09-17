test_that("DFP benchmark on R function", {
    if (!getOption("BhatF.run_benchmarks", default=FALSE)) {
        skip("To run benchmarks, set option `BhatF.run_benchmarks` to TRUE")
        return()
    }

    nvars = 30
    xmin = -5
    xmax = 5
    aa = runif(nvars, min=1, max=10)
    bb = runif(nvars, min=0.5*xmin, max=0.5*xmax)
    cc = 5
    est = runif(nvars, min=0.1*xmin, max=0.1*xmax)

    cat("aa=",aa,"\n",sep=" ")
    cat("bb=",bb,"\n",sep=" ")
    cat("est=",est,"\n",sep=" ")

    
    fn = function(x) {
        for (i in 1:nvars) {
            cc = cc + aa[i]*(x[i]-bb[i])^2
        }
        return(cc)
    }

    x = list(label=paste0("x", 1:nvars),
             est=est,
             low=rep(xmin, nvars),
             upp=rep(xmax, nvars))
    
    #suppressMessages(
        result<-dfp(x, fn)
    #)

    cat("Convergence:",result$status,"\n")
    cat("Rurun...\n")
    
    x$est = result$est
    #suppressMessages(
        result<-dfp(x, fn)
    #)

    expect_equal(result$fmin, cc)
    expect_equal(result$est, bb, tolerance=1e-5)
    expect_equal(result$label, x$label)
    expect_equal(result$status, 0)
})
