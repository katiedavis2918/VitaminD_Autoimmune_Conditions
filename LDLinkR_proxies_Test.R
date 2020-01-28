
# ----------- LDLinkR package to find proxies for a list of SNPs ----------

# Populations: EUR = European; FIN = Finnish in Finland

my_proxies <- LDproxy(snp = "rs10741657", 
                      pop = "EUR", 
                      r2d = "r2",
                      token = "c7e34cfd1059")

# Subset of proxies which have R2>0.8
LD_proxies <-subset(my_proxies, my_proxies$R2>0.8)

# grep function to check for top snps in output 

#--------------------

# Test list of SNPs (batch)

snps_to_test<- c("rs10741657", "rs10745742", "rs12785878")

my_proxies_batch <- LDproxy_batch(snp = snps_to_test, 
                                  pop = "EUR", 
                                  r2d = "r2",
                                  token = "c7e34cfd1059")

#--------------------