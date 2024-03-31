# DSIR
Code of Deep structure integrative representation of multi-omics data for cancer subtyping. DSIR is a cancer subtype framework based on multi-omics profiles. The input for the framework is mRNA, miRNA and DNA methylation. The output is the corresponding similarity matrix between samples. 

learn similarity matrix
```
#train.py is used to obtain global matrix to learn consensus similarity
#the models are saved in train_models
```

DSIR is mainly divided into two components: 1. networks is used to learn consensus similarity matrix between samples. 2. Consensus clustering determine the number of subtypes and the cluster label corresponding to each sample. 


DSIR's Consensus clustering module:  
```
clustering_DSIR/DSIR_.R
# record the corresponding class label for each sample and the output file is results_BIC
```

the program of the comparison methods :
 ```
 ./enrichment_test
 ```

DSIR is based on the Python program language and R language. The network's implementation was based on the open-source library Pytorch 1.9.0. 
