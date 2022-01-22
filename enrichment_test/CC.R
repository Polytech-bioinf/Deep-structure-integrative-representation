# clear enviroment
rm(list = ls())
# -----------------------------------
# load library
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
cancername <- "AML"
classi <- 3
parentpath <- "/Users/yangyan/Desktop/data/AML/"
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

file_dir = "/Users/yangyan/Desktop/enrichment_test/AML/CC/CC_param_comb/Labels/"
save_dir_1 = paste0(file_dir, "label.txt")
save_dir_2 = "/Users/yangyan/Desktop/enrichment_test/AML/CC/surv_pval.txt"


datas=list(GeneExp=gene,MethyExp=methy,miRNAExp=mirna)
resultCC=ExecuteCC(clusterNum=classi, #cluster number
                   d=datas, # data to be clustered
                   maxK = 10, clusterAlg = "hc",
                   distance = "pearson", title = "ConsensusClusterResult", reps = 20,
                   pItem = 0.8, pFeature = 1, plot = "png", innerLinkage = "average",
                   finalLinkage = "average", writeTable = FALSE, weightsItem = NULL,
                   weightsFeature = NULL, verbose = FALSE, corUse = "everything")
# add group column to survival data
surv$Labels<-resultCC$group
write.table(surv$Labels, save_dir_1, row.names = F, col.names = F)


survresult<-survdiff(Surv(Survival,Death)~Labels, data=surv)
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)

print(p.val)
write.table(p.val, save_dir_2, row.names = F, col.names = F)


