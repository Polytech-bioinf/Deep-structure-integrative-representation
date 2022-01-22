## ----message=FALSE------------------------------------------------------------

library(CancerSubtypes)
library(survival)

# load the data
# --------------------------------------------
# change the cacer name
# 1 "BREAST"
# 2 "COLON"
# 3 "GLIO"
# 4 "KIDNEY"
# 5 "LUNG"
# 6 "OVARY"

aa="AML"
cancername <- "AML"
classi <- 9
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


file_dir = paste0("/Users/yangyan/Desktop/enrichment_test/",aa,"/spectral/spectral_param_comb/Labels/")
save_dir_1 = paste0(file_dir, "label.txt")
save_dir_2 = paste0("/Users/yangyan/Desktop/enrichment_test/",aa,"/spectral/surv_pval.txt")


data=list(gene, methy, mirna)
# print(nrow(data[[1]]))


## -----------------------------------------------------------------------------

concat.omics = do.call(rbind, data)
similarity.data = affinityMatrix(dist2(as.matrix(t(concat.omics)),
                                       as.matrix(t(concat.omics))), 
                                 20, 0.5)
num.clusters = estimateNumberOfClustersGivenGraph(similarity.data, 
                                                  classi)[[3]]  
clustering = spectralClustering(similarity.data, num.clusters)
Labels <- clustering
#print(Labels)
write.table(Labels, save_dir_1, row.names = F, col.names = F)

# Labels_dir = paste0(file_dir, "label.txt")
# Labels <- read.table(Labels_dir)
# dim(Labels)
survresult<-survdiff(Surv(Survival,Death)~Labels, data=surv)
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
print(p.val)
write.table(p.val, save_dir_2, row.names = F, col.names = F)


