library(aSPU)
setwd("~/Google_Drive/MyPackages/WaSPU/data")


tissue="TW_Artery-Coronary_elasticNet0_0.5.db.txt";
method = "perm"
B = 1e3

genename =  "ANKRD34A"



db = readRDS(paste0(tissue,".rds")); db = subset(db,alpha==0.5)
allgenes = unique(db$gene)

system("mkdir -p temp")


mysnps = subset(db,gene==genename)
write.table(mysnps$rs,paste0("temp/mysnps_",genename,".txt"),row.names=F,col.names=F,quote=F,append=F)
write.table(mysnps[,c("rs","refAllele")],paste0("temp/refAllele_",genename,".txt"),row.names=F,col.names=F,quote=F,append=F)
system(paste0("./plink --bfile wgas_maf5 --extract ","temp/mysnps_", genename,".txt --recodeA --recode-allele ","temp/refAllele_",genename,".txt --out ","temp/genotype_",genename))

dat = read.table(paste0("temp/genotype_",genename,".raw"),header=T)

system(paste0("rm ","temp/genotype_",genename,"*"))
system(paste0("rm ","temp/mysnps_",genename,"*"))
system(paste0("rm ","temp/refAllele_",genename,"*"))

#### extract X0:the un-weighted SNPs and Y
y = readRDS("binary_phenotype.rds")# convert to 1,0 coding
y = rnorm(90)
X0 = as.matrix(dat[,7:ncol(dat)])

#### extract the weights from PrediXcan database
snp.extract = names(dat)[7:ncol(dat)]
snp.extract = sub("_.*", "", snp.extract) # replace everything from the start of the string to the "_" with ""


beta <- mysnps[match(snp.extract,mysnps$rs),"beta"]

#### get X: the weighted SNPs

X = X0
X.sum = rep(0,nrow(X))
for(i in 1:ncol(X0)){
  X[,i] <- X0[,i]*beta[i]
  X.sum = X.sum + X[,i]
}

### aSPU test based on weighted X
p.aSPU = aSPU(y,X,resample=method,pow = pow,model="binomial",n.perm=B)$pvs
p.aSPU

p.aSPU = aSPU(y,X,resample=method,pow = pow, model="gaussian",n.perm=B)$pvs

p.aSPU






