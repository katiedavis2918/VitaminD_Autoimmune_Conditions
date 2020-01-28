#-------------------------------------------------------------------------------------------------
# "Generic script" 0.1:    ADD Effect Allele Frequencies (EAFs)
#-------------------------------------------------------------------------------------------------

# Pull Allele Frequencies using Biomart

devtools::install_github("JhuangLab/annovarR")
library(annovarR)

dat<-(ExposureX)
names(dat)<-c("")



# if ref allele == minor allele, EAF = MAF
# if ref allele == major allele, EAF = 1- MAF