---
title: "A Tutorial"
output: html_document
---



## A Quick Tutorial on Conducting WaSPU Test

---

### A Weighted Adaptive Sum of Powered Score Test (WaSPU)

This is an R package for performing association study by integrating Genomic and Imaging Endophenotypes with individual/summary level data. Please cite the following manuscript for WaSPU methods:

---

**Xu Z., Xu Gong, Pan W. Integrating genomic and imaging endophenotypes in GWAS**

---

### Outline
1. Installation
2. Association Testing with Individual Level Data
3. Association Testing with Summary Statistics
4. How to get a correlation matrix of a set of SNPs using Hapmap as reference panel
---

### Installation
1. Fire up R 

2. install devtools by:

```R
install.packages("devtools")
```

3. install WaSPU package from github:

```R
install.packages("aSPU")
library(devtools)
install_github("jasonzyx/WaSPU")
```


### Association Testing with Individual Level Data


```r
library(WaSPU)
set.seed(123)
Y = sample(c(1,0), 100, replace = T)
X = matrix(sample(c(0,1,2), 1000, replace = T), nrow = 100)
Z = matrix(rnorm(200), nrow = 100)
weights = runif(10)
WaSPU(Y, X, Z, weights, pow = c(1:6, Inf), n.perm = 1e3)
```

```
## $paSPU
##      SPU1      SPU2      SPU3      SPU4      SPU5      SPU6      SPU7 
## 0.4350000 0.5490000 0.5210000 0.4790000 0.3920000 0.4330000 0.3800000 
##      SPU8      SPU0      aSPU 
## 0.4050000 0.3820000 0.6203796 
## 
## $pSSU
## [1] 0.5044016
## 
## $pSum
## [1] 0.4167493
## 
## $pUminP
## [1] 0.1793441
```

### Association Testing with Summary Statistics


```r
library(WaSPU)
set.seed(123)
Zstat = runif(10,-2,2)
corSNP = 0.1 + diag(0.9,10)
weights = runif(10)
WaSPUs(Zstat, corSNP, weights, pow = c(1:6, Inf), n.perm = 1e3)
```

```
## $paSPU
##     SPUs1     SPUs2     SPUs3     SPUs4     SPUs5     SPUs6   SPUsInf 
## 0.6740000 0.3060000 0.3020000 0.2790000 0.2580000 0.2560000 0.2400000 
##     aSPUs 
## 0.3936064 
## 
## $pSSU
## [1] 0.3013146
## 
## $pSum
## [1] 0.6524539
## 
## $pUminP
## [1] 0.5028263
```

### Calculate a correlation matrix of a set of SNPs using Hapmap as reference panel

1. Download hapmap reference data, Plink software and an example file "ex_mysnps.txt"

```
wget https://github.com/jasonzyx/WaSPU_resource/archive/master.zip
```

2. uncompress files

```
unzip WaSPU_resource-master.zip
rm WaSPU_resource-master.zip
cd WaSPU_resource-master
```

3. some system may need to change accessbility of plink, in order to use.

```
chmod 777 plink
```

4. fire up R

```
source("impute.R")
system("mkdir temp")
mysnps <- read.table("ex_mysnps.txt", header = T)
mysnps$SNP_map = paste0(mysnps$CHR,":",mysnps$BP)
write.table(mysnps$SNP_map,"temp/mysnps.txt",row.names=F,col.names=F,quote=F,append=F)
write.table(mysnps[,c("SNP_map","refAllele")],"temp/refAllele.txt",row.names=F,col.names=F,quote=F,append=F)
system(paste0("./plink --bfile hapmap_CEU_r23a_hg19", 
              " --extract ", "temp/mysnps.txt --recodeA --recode-allele ",
              "temp/refAllele.txt --out ","temp/ex"))
dat = read.table("temp/ex.raw",header=T)
X0 = dat[,7:ncol(dat)]
X00 = as.matrix(X0)
X0 = impute(X0)
X0 = as.matrix(X0)
### remove SNPs w/o variation
X0_sd <- apply(X0,2,sd)
X0 <- X0[,X0_sd!=0]


snp.extract = colnames(X0)
for(l in 1:length(snp.extract)) snp.extract[l] <- substr(snp.extract[l],2,nchar(snp.extract[l])-2)   
snp.extract = gsub("[.]", ":", snp.extract)

colnames(X0) <- snp.extract

commonSNPs = intersect(snp.extract, mysnps$SNP_map)

mysnps <- mysnps[match(commonSNPs,mysnps$SNP_map),]  
X0 <- X0[,commonSNPs] 

corMat <- cor(X0)

corMat[1:5, 1:5]
```

The object "corMat" is the correlation matrix.


