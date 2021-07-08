library(Rcpp)

xx = matrix(byrow=TRUE,
            data=c(0, 10, 0, 0, 100,
                   0, 10, 0, 0, 20,
                   0, 0.01, 0, 0, 0.1),
            nrow=3, ncol=5,
            dimnames=list(c("A","B","C")))



## requires a crafted makefile
Sys.setenv(PKG_LDFLAGS="-L/home/galexv/Work/UMich/CISNet/LCSModel/Bhat/Bhat_f95 -lBhat_f95")

## compiles, but does not load successfully:
sourceCpp(file="../call_migrad.cpp", verbose=TRUE, cacheDir="./tmp")

## still needs a manual compilation:
##    LD_RUN_PATH=/home/galexv/Work/UMich/CISNet/LCSModel/Bhat/Bhat_f95/
##    g++ ..... -lBhat_f95

