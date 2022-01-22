## ----message=FALSE------------------------------------------------------------
library('iClusterPlus')
library(CancerSubtypes)
library(survival)
biocLite("iClusterPlus")
# load the data
# --------------------------------------------
# change the cacer name
# 1 "BREAST"
# 2 "COLON"
# 3 "GLIO"
# 4 "KIDNEY"
# 5 "LUNG"
# 6 "OVARY"

aa="LSCC"
cancername <- "LUNG"
classi <- 4
parentpath <- paste0("/Users/yangyan/Desktop/data/",aa,"/")
genefile <- paste(parentpath, cancername, "_Gene_Expression.txt", sep="")
methyfile <- paste(parentpath, cancername, "_Methy_Expression.txt", sep="")
mirnafile <- paste(parentpath, cancername, "_Mirna_Expression.txt", sep="")
survfile <- paste(parentpath, cancername, "_Survival.txt", sep="")
# read mRNA
gene <- read.table(genefile, header = T)
dim(gene)
# read Methylatio
methy <- read.table(methyfile, header = T)
dim(methy)
# read miRNA
mirna <- read.table(mirnafile, header = T)
dim(mirna)
# read survival
surv <- read.table(survfile, header = T)
dim(surv)

# -----------------------------------
# ExecuteCC
# names is not same, check by names(gene) == names(methy)
# change to same first
names(methy) <- names(gene)
names(mirna) <- names(gene)

gene = as.matrix(gene)
print(nrow(gene))
methy = as.matrix(methy)
mirna = as.matrix(mirna)


file_dir = paste0("/Users/yangyan/Desktop/enrichment_test/",aa,"/icluster/icluster_param_comb/Labels/")
save_dir_1 = paste0(file_dir, "label.txt")
save_dir_2 = paste0("/Users/yangyan/Desktop/enrichment_test/",aa,"/icluster/surv_pval.txt")


data=list(gene, methy, mirna)
# print(nrow(data[[1]]))


## -----------------------------------------------------------------------------
omics.list = data
MAX.NUM.CLUSTERS=5
icluster.rets = list()

if (length(omics.list) == 1) {
  icluster.ret = iClusterPlus::tune.iClusterBayes(cpus=(MAX.NUM.CLUSTERS - 1), t(omics.list[[1]]), 
                                                  K=1:(MAX.NUM.CLUSTERS - 1), type=c('gaussian'))$fit
} else {
  icluster.ret = iClusterPlus::tune.iClusterBayes(cpus=(MAX.NUM.CLUSTERS - 1), t(omics.list[[1]]), 
                                                  t(omics.list[[2]]), 
                                                  t(omics.list[[3]]), 
                                                  K=1:(MAX.NUM.CLUSTERS - 1), type=rep('gaussian', 3))$fit
}
dev.ratios = lapply(1:(MAX.NUM.CLUSTERS - 1), function(i) icluster.ret[[i]]$dev.ratio)


print('dev.ratios are:')
print(dev.ratios)

optimal.solution = icluster.ret[[which.max(dev.ratios)]]
Labels <- optimal.solution$clusters
#print(Labels)
write.table(Labels, save_dir_1, row.names = F, col.names = F)

# Labels_dir = paste0(file_dir, "label.txt")
# Labels <- read.table(Labels_dir)
# dim(Labels)
survresult<-survdiff(Surv(Survival,Death)~Labels, data=surv)
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
print(p.val)
write.table(p.val, save_dir_2, row.names = F, col.names = F)


