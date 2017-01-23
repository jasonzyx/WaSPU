# work.dir = "~/Google_Drive/MyPackages/WaSPU/data"
# genename = "ANKRD34A"
# GWAS.plink = "wgas_maf5";
# Weight.db = "TW_WholeBlood_ElasticNet.0.5.db"
# method = "perm"; model = "binomial"
# B = 1e3; pow = c(1:8, Inf)

y = readRDS("binary_phenotype.rds")# convert to 1,0 coding
WaSPU(work.dir = "~/Google_Drive/MyPackages/WaSPU/data",
      # genename= "ANKRD34A",
      genename= "APOE",
      GWAS.plink = "wgas_maf5",
      y=y,
      Weight.db = "TW_WholeBlood_ElasticNet.0.5.db",
      method = "perm", model = "binomial",
      B = 1e3, pow = c(1:8, Inf))
