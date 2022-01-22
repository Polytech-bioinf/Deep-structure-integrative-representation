
par(mfrow=c(1,1))

Algorithm <- c("SNF","CC","SNFCC","PINS","LRAcluster","MCCA","K-Means", "Spectral", "NEMO","iCluster","DSIR")
# parameter for LSCC
surv_pval <- c(1.9, 1.4,2.1,1.5, 1.7, 2.4, 1.5,1.7,2.4,2.2,2.7)
num_enrich <- c(1.5,1.5,1.8,1.1,1.8,1.3,1.5,1.3,1.7,1.6,2.1)


alg_cols <- c( "cyan","gray","cadetblue1", "magenta", "blue", "green", "pink","royalblue","tomato", "azure3","red")
# print(alg_cols)
pch = c( 16, 10,12,17, 18, 2, 3,4,6,5,15)

available.indices = !is.na(surv_pval)
surv_pval = surv_pval[available.indices]
num_enrich = num_enrich[available.indices]
current.cols = alg_cols[available.indices]
current.pch = pch[available.indices]

surv_significance <- -log10(0.05)


ylab = 'Mean -log10 of logrank p-value'
xlab = 'Mean number of enriched clinical parameters'
 

plot(num_enrich,surv_pval, xlab=xlab, ylab=ylab, ylim=c(0.5, max(surv_pval, surv_significance)+0.2), 
    xlim=c(0.5, max(num_enrich+0.2)), col = alg_cols, pch=pch, cex.lab = 1, cex=1, cex.axis=1, cex.main=1.6, lwd=2,xaxt="n",yaxt="n")

#abline(v=surv_significance, col="red")

#legend("bottomleft", Algorithm, pch=pch, col=alg_cols)
legend("bottomleft",Algorithm, pch=pch, col=alg_cols,cex=0.5)
grid(30,30, lwd = 0.5,col = "lightgray")
axis(1,at=c(0.5,1.0,1.5,2.0,2.5))
axis(2,at=c(0.5,1.0,1.5,2.0,2.5))
abline(v=2.1, col="red",lty=3)
abline(h=2.7, col="red",lty=3)





