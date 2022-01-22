DSIR.num.clusters <- function(W, NUMC=2:15) {
  if (min(NUMC) == 1) {
    warning("Note that we always assume there are more than one cluster.")
    NUMC = NUMC[NUMC > 1]
  }
  W = (W + t(W))/2
  diag(W) = 0
  if (length(NUMC) > 0) {
    degs = rowSums(W)
    degs[degs == 0] = .Machine$double.eps
    D = diag(degs)
    L = D - W
    Di = diag(1/sqrt(degs))
    L = Di %*% L %*% Di
    print(dim(L))
    eigs = eigen(L)
    eigs_order = sort(eigs$values, index.return = T)$ix
    eigs$values = eigs$values[eigs_order]
    eigs$vectors = eigs$vectors[, eigs_order]
    eigengap = abs(diff(eigs$values))
    eigengap = (1:length(eigengap)) * eigengap
    
    t1 <- sort(eigengap[NUMC], decreasing = TRUE, index.return = T)$ix
    return(NUMC[t1[1]])
  }
}