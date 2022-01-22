# install.packages("LRAcluster")
library("LRAcluster")
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
cancername <- "BREAST"
classi <- 5
dime = 5

parentpath <- "/Users/yangyan/Desktop/data/BIC/"
genefile <- paste(parentpath, cancername, "_Gene_Expression.txt", sep="")
methyfile <- paste(parentpath, cancername, "_Methy_Expression.txt", sep="")
mirnafile <- paste(parentpath, cancername, "_Mirna_Expression.txt", sep="")
survfile <- paste(parentpath, cancername, "_Survival.txt", sep="")

file_dir <- "/Users/yangyan/Desktop/enrichment_test/BIC/LRAcluster/LRAcluster_param_comb/Labels/"
save_dir_1 = paste0(file_dir, "label.txt")
save_dir_2 = "/Users/yangyan/Desktop/enrichment_test/BIC/LRAcluster/surv_pval.txt"


gene <- read.table(genefile, header = T)
methy <- read.table(methyfile, header = T)
mirna <- read.table(mirnafile, header = T)
surv <- read.table(survfile, header = T)

names(methy) <- names(gene)
names(mirna) <- names(gene)

gene = as.matrix(gene)
methy = as.matrix(methy)
mirna = as.matrix(mirna)



#data = list(m, mi, me)
data = list(gene, methy, mirna)
#types = list("gaussian", "gaussian")
#types = list("gaussian", "gaussian", "gaussian")
types = list("gaussian", "gaussian", "gaussian")

names = list("Gene", "Methy", "Mirna")

res = LRAcluster(data,types, dimension=dime, names=as.character(1:length(data)))

result = t(res$coordinate)
kc <- kmeans(result, classi)

Labels <- kc$cluster
# print(Labels)

write.table(Labels, save_dir_1, row.names = F, col.names = F)


survresult<-survdiff(Surv(Survival, Death)~Labels, data=surv)
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
print(p.val)
write.table(p.val, save_dir_2, row.names = F, col.names = F)
