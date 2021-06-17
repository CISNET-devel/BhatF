
y _ read.table('mcmc.out',header=TRUE,quote="")

npar <-  ncol(y)-2
par(mfrow=c(npar+1,1),bg="grey")

par(mar=c(0, 5, 0.1, 4) + 0.1)
cat('90% Credibility Regions','\n')
plot(y[,2], type='l', xlab = " ", ylab = "-log-density",col=2)
for (i in 1:(npar-1)) {
  par(mar=c(0, 5, 0, 4) + 0.1)
  plot(y[,2+i], type='l', xlab = " ", ylab = names(y)[i+2], col=2)
  tmp <- quantile(y[,2+i],probs=c(.025,.5,.975))
  cat(names(y)[i+2],tmp,'\n')
}
par(mar=c(0.1, 5, 0, 4) + 0.1);
plot(y[,2+npar], type='l', xlab = " ", xaxt='n', ylab = names(y)[npar+2], col=2);
tmp <- quantile(y[,2+npar],probs=c(.025,.5,.975));
cat(names(y)[npar+2],tmp,'\n')
