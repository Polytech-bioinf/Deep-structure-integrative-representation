# DSIR
Code of Deep structure integrative representation of multi-omics data for cancer subtyping. DSIR is a cancer subtype framework based on multi-omics profiles. The input for the framework is mRNA, miRNA and DNA methylation. The output is the corresponding similarity matrix between samples. 

The input for the framework
```
#the input of the data is finished by ./dataset.py and ./dataloader.py
```
compute local similarity matrix
```
#the compute of local similarity matrix is finished by ./snf
```
train global similarity matrix
```
#networks.py is used to obtain global matrix to learn consensus similarity
#the models are saved in ./train_models
```

DSIR is mainly divided into two components: 1. networks is used to learn consensus similarity matrix between samples. 2. Consensus clustering determine the number of subtypes and the cluster label corresponding to each sample. 
```
# the input of BIC omics data set is txt. We can use the following command to finish the subtyping process: 
python train_C.py -m BIC -i ./DATA/BIC/..txt  
# the out of the process is the consensus similarity
# the consensus similarity matries are stored in ./C/BIC
```

DSIR's Consensus clustering module is used as follows:  
```
R clustering_DSIR/DSIR_.R
# record the corresponding class label for each sample and the output file is ./results_BIC
```

the program of the comparison methods :
 ```
 ./enrichment_test
 ```

DSIR is based on the Python program language and R language. The network's implementation was based on the open-source library Pytorch 1.9.0. 


