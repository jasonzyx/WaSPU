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


