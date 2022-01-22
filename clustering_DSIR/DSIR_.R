#source("/Users/yangyan/Desktop/DSIR/clustering_DSIR/DSIR.num.clusters.R")
library("xlsx")
library("survival")

file_dir = "/Users/yangyan/Desktop/test4/COAD/"
param_dir = paste0(file_dir, "_param_comb/")
enrich_dir = paste0(file_dir, "_enrich_surv/")
file_dir_1 = paste0(param_dir,"/param2_10")
file_name = paste0(file_dir_1, "/param3_20/epoch_2/matrix.txt")
save_label_dir=paste0(file_dir_1, "/param3_20/epoch_2/labels.xlsx")
p_value <- c()
C <- read.table(file_name, header = F)
C = as.matrix(C)
S=0.5*(C+t(C))
cluster_num=DSIR.num.clusters(C)

#cluster_num=6
clustering=spectralClustering(S,cluster_num)
write.xlsx(clustering, sheetName="Sheet1",save_label_dir, row.names = T)

surv_name = paste0(enrichdata_dir,"_Survival.txt")
print(file_name)

survdata<-read.table(file = surv_name, header = T)
survdata$labels=clustering
survdata$Survival=as.numeric(survdata$Survival)
survdata$Death=as.numeric(survdata$Death)
survdata$labels=as.numeric(survdata$labels)
survresult <- survdiff(Surv(Survival,Death)~labels, data=survdata)
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
print(survresult)
print(p.val)
write.table(computep, paste0(enrich_dir,"/param2_10/param3_20/surv_pval", ".txt"), row.names = F, col.names = F,append=TRUE)