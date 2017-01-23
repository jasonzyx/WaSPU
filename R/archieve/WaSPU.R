source("~/Google_Drive/MyPackages/WaSPU/R/impute.R")

WaSPU <- function(work.dir, genename, GWAS.plink = "wgas_maf5", y, Weight.db, method = "perm", model = "binomial",B = 1e3, pow = c(1:8, Inf)){
  # test parameters
  # work.dir = "~/Google_Drive/MyPackages/WaSPU/data"
  # genename = "ANKRD34A"
  # GWAS.plink = "wgas_maf5";
  # Weight.db = "TW_WholeBlood_ElasticNet.0.5.db"
  # method = "perm"; model = "binomial"
  # B = 1e3; pow = c(1:8, Inf)
  # y = readRDS("binary_phenotype.rds")# convert to 1,0 coding

  setwd(work.dir)

  cat("loading weight database ... \n")
  db = read.table(paste0(Weight.db,".txt"),header = T); db = subset(db,alpha==0.5)
  allgenes = unique(db$gene)

  system("mkdir -p temp")


  mysnps = subset(db,gene==genename)
  write.table(mysnps$rsid,paste0("temp/mysnps_",genename,".txt"),row.names=F,col.names=F,quote=F,append=F)
  write.table(mysnps[,c("rsid","refAllele")],paste0("temp/refAllele_",genename,".txt"),row.names=F,col.names=F,quote=F,append=F)
  system(paste0("./plink --bfile wgas_maf5 --extract ","temp/mysnps_", genename,".txt --recodeA --recode-allele ","temp/refAllele_",genename,".txt --out ","temp/genotype_",genename))

  is.genotype.avail = system(paste0("ls temp/genotype_",genename,".raw"), intern = TRUE)
  if(length(is.genotype.avail)==0) stop("no matching SNPs in GWAS data")
  dat = read.table(paste0("temp/genotype_",genename,".raw"),header=T)

  system(paste0("rm ","temp/genotype_",genename,"*"))
  system(paste0("rm ","temp/mysnps_",genename,"*"))
  system(paste0("rm ","temp/refAllele_",genename,"*"))

  #### extract X0:the un-weighted SNPs and Y

  X0 = as.matrix(dat[,7:ncol(dat)])
  if(sum(is.na(X0))>0) cat("imputing missing genotype by mean ... \n")

  X0 = impute(X0)
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

  cat(paste0(ncol(X)," SNPs were included in final analysis \n"))
  cat("calculating P-values ... \n")
  ### aSPU test based on weighted X
  p.aSPU = aSPU(y,X,resample=method,pow = pow,model="binomial",n.perm=B)$pvs
  names(p.aSPU)[names(p.aSPU)=="SPU0"] <- "SPUInf"
  p.aSPU

}
