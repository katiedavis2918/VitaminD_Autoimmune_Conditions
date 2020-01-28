#-----------------------------------------------------------------------------------------------------
# 6 November 2019
# Load Test Dataset and flip beta's in exposure dataframe:

# NB NOTE: need to create standard column order in all datasets!

# setwd()
#-----------------------------------------------------------------------------------------------------

# Jiang et al 6 SNPs for Vitamin D used as exposure:
Exposure_Test<- read.table(file="VITD_SAMPLE_6SNPs", header = FALSE)

#Adding names of col's 
##For vitamin D dataframe:
names(Exposure_Test) <- c("SNP", "EffectAllele", "OtherAllele", "Effect", "StdErr", "P.value", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "Position", "Chromosome", "Sample.size")

#Select col's 
Exposure_Test<- Exposure_Test[ ,c("SNP", "EffectAllele", "OtherAllele", "Effect", "StdErr", "P.value", "Position", "Chromosome", "Sample.size")]

# --------

# Add column of EAF for vitD GWAS signif SNPs

## Effect Allele Frequency values put into a vector

SNPs.vitD<- c("rs12785878", "rs3755967", "rs8018720", "rs10745742", "rs10741657", "rs17216707")
#eaf from 1000 genomes EUR population not this cohort:
eaf.vitD<- c(0.701, 0.248, 0.827, 0.418, 0.381, 0.773)

# add eaf column to the VITDdf_6SNPS dataframe 
df_eaf<- data.frame(SNPs.vitD,eaf.vitD)

#Changing names of col's 
names(df_eaf) <- c("SNP", "eaf.exposure")

#merge eaf df with VITDdf_6SNPS to add this as a column
Exposure_Test<- merge(Exposure_Test,df_eaf, by= "SNP")

#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

# MIP1B dataset used as outcome:
Outcome_Test<- read.table(file="MIP1Bdf_SAMPLE_6SNPs", header = FALSE)

#Adding names of col's 
##For mip1B dataframe:
names(Outcome_Test) <-c("SNP", "Chromosome", "Position", "OtherAllele", "EffectAllele", "Effect", "StdErr", "Direction", "P.value", "HetPVal")

# Select col's
Outcome_Test<- Outcome_Test[, c("SNP", "Chromosome", "Position", "OtherAllele", "EffectAllele", "Effect", "StdErr", "P.value")]

# --------

# Add column of EAF for vitD SNPs in MIP1B DATASET:

EAFs_Fin_6SNPs<- read.csv(file="VitD_6SNPs_FinEAF.csv", header=TRUE)

names(EAFs_Fin_6SNPs) <-c("SNP", "eaf.outcome")

Outcome_Test<-merge(Outcome_Test,EAFs_Fin_6SNPs, by="SNP")


#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

# Standardise beta's (all positive in exposure dataset)

#------------------
df<-Exposure_Test

df$EffectAllele <- as.character(df$EffectAllele)        # Change from factor to character otherwise "NA" values generated in next coloumn switching steps 
df$OtherAllele <- as.character(df$OtherAllele)

# df[i,4] = beta / effect
# df[i,2] = EffectAllele
# df[i,3] = OtherAllele
# df[i,11] = help.col (equal to EffectAllele)
# df[i,10] = EAF column
#------------------

f_out<-rep(0,6)

for(i in 1:6){
  
 if (df[i,4] <0)
  
 {f_out<- (df[i,4] = abs(df[i,4]))
  
  df$help.col = df[,"EffectAllele"]             # Adds (duplicated) 'EffectAllele' column for switch step
  
  df[i,2] = df[i,3]                             # switch: df[,"EffectAllele"] = df[,"OtherAllele"]
  
  df[i,3] = df[i,11]                            # switch: df[,"OtherAllele"] = help.col
  
  df[i,10] = 1- df[i,10]} }                     # subtract: EAF = 1- EAF

# remove help column
df<- subset(df, select = -c(help.col))

#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------




