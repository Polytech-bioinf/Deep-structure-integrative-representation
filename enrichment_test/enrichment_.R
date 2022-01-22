# install.packages("xlsx", repos="https://mran.microsoft.com/snapshot/2019-02-01/")
library("xlsx")


cal.age.enrichment <- function(age, cluster){
  # ??????????KW????
  # cluster <- read.xlsx("D:\\enrichment_test\\lung\\lung_clustering.xlsx", sheetIndex=1)
  # age <- read.xlsx("D:\\enrichment_test\\lung\\lung_age.xlsx", sheetIndex=1)
  age.test.result = kruskal.test(as.numeric(unlist(age)), as.numeric(unlist(cluster)))
  
  return(age.test.result)
}

cal.discrete.enrichment <- function(file_dir){
  # ??????ɢ?Ͳ????ĸ???????
  para_tb = table(as.data.frame(file_dir[,2:3]))
  print(para_tb)
  enrichment.result <- chisq.test(para_tb)
  return(enrichment.result)
}
aa="MELANOMA"
bb="MELANOMA"
#file_dir = "/Users/yangyan/Desktop/enrichment_test/"+a+"/SNF/SNF_param_comb/"
file_dir = paste0("/Users/yangyan/Desktop/enrichment_test/",aa,"/SNFCC/SNFCC_param_comb/")
#param_dir = paste0(file_dir, "_param_comb/model5/")
enrich_dir = paste0("/Users/yangyan/Desktop/enrichment_test/",aa,"/SNFCC/_enrich_surv/")
  
for (i in (1:126)){
  # age enrichment
  age = paste0("/Users/yangyan/Desktop/test/",aa,"/_param_comb/",bb,"_Age.xls")
  age_file <- read.xlsx(age, sheetIndex = 1)
  # print(age_file)
  
  cluster_name = paste0(file_dir, "Labels/label",i,".xls")
  #cluster_name = paste0(param_dir, "labels1.xlsx")
  print(cluster_name)
  cluster_file <- read.xlsx(cluster_name, sheetIndex = 1)
  a=cluster_file
  b=cluster_file
  c=cbind(a,b)
  cluster_file=c
  age_test_res = cal.age.enrichment(age_file, cluster_file)
  print(age_test_res)
  p_val = as.matrix(age_test_res)
  p_value = as.data.frame(p_val[3,])
  age_dir = paste0(enrich_dir,"age_res",".txt")
  write.table(p_value,age_dir,append = T, row.names = F, col.names = F)
  
  print("_______________________________________________________________________________")
  
  
  # gender enrichment analysis
  #gender_name = paste0(param_dir, "demo1.xlsx")
  gender_name = paste0(file_dir, "demo",i,".xlsx")
  print(gender_name)
  
  gender_file = read.xlsx(gender_name, sheetIndex = 1)
  gender.pval = cal.discrete.enrichment(gender_file)
  print(gender.pval)
  p_val = as.matrix(gender.pval)
  p_value = as.data.frame(p_val[3,])
  
  gender_dir = paste0(enrich_dir,"gender_res",".txt")
  write.table(p_value,gender_dir,append = T, row.names = F, col.names = F)
  

  print("_______________________________________________________________________________")
  
  #path_M enrichment analysis
  #path_M_name = paste0(param_dir, "demo1.xlsx")
  path_M_name = paste0(file_dir, "demo",i,".xlsx")
  print(path_M_name)
  
  pathM_file = read.xlsx(path_M_name, sheetIndex = 2)
  pathM.pval = cal.discrete.enrichment(pathM_file)
  print(pathM.pval)
  p_val = as.matrix(pathM.pval)
  p_value = as.data.frame(p_val[3,])
  
  pathM_dir = paste0(enrich_dir,"pathM_res.txt")
  write.table(p_value,pathM_dir,append = T, row.names = F, col.names = F)
  print("_______________________________________________________________________________")
  
  
  # path_N enrichment analysis
  #path_N_name = paste0(param_dir, "demo1.xlsx")
  path_N_name = paste0(file_dir, "demo",i,".xlsx")
  print(path_N_name)
  
  pathN_file = read.xlsx(path_N_name, sheetIndex = 3)
  pathN.pval = cal.discrete.enrichment(pathN_file)
  print(pathN.pval)
  p_val = as.matrix(pathN.pval)
  p_value = as.data.frame(p_val[3,])
  
  pathN_dir = paste0(enrich_dir,"pathN_res.txt")
  write.table(p_value, pathN_dir, append = T, row.names = F, col.names = F)
  
  print("_______________________________________________________________________________")
  
  
  # path_T enrichment analysis
  #path_T_name = paste0(param_dir, "demo1.xlsx")
  path_T_name = paste0(file_dir, "demo",i,".xlsx")
  print(path_T_name)
  
  pathT_file = read.xlsx(path_T_name, sheetIndex = 4)
  pathT.pval = cal.discrete.enrichment(pathT_file)
  print(pathT.pval)
  p_val = as.matrix(pathT.pval)
  p_value = as.data.frame(p_val[3,])
  
  pathT_dir = paste0(enrich_dir,"pathT_res.txt")
  write.table(p_value,pathT_dir,append = T, row.names = F, col.names = F)
  
  print("_______________________________________________________________________________")
  
  
  #path_stage_name = paste0(param_dir, "demo1.xlsx")
  path_stage_name = paste0(file_dir, "demo",i,".xlsx")
  print(path_stage_name)
  
  path_stage_file = read.xlsx(path_stage_name, sheetIndex = 5)
  path_stage.pval = cal.discrete.enrichment(path_stage_file)
  print(path_stage.pval)
  p_val = as.matrix(path_stage.pval)
  p_value = as.data.frame(p_val[3,])
  
  path_stage_dir = paste0(enrich_dir,"path_stage_res.txt")
  
  write.table(p_value,path_stage_dir,append = T, row.names = F, col.names = F)
  
}







