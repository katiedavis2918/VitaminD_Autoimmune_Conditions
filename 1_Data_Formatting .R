#-------------------------------------------------------------------------------------------------
# "Generic script" 1:    DATA FORMATTING 
#-------------------------------------------------------------------------------------------------

# Output from previous script: 'ExposureX' and 'OutcomeX' tables containing selected instruments
# from exposure and outcome datasets.

library(MendelianRandomization)
library(TwoSampleMR)

#---------    Put dataframes into required format for TwoSampleMR    -----------------------------

#For test datasets used:
# ExposureX<-df
# OutcomeX<-Outcome_Test

Exposure_DataX <- format_data( ExposureX,
                               type= "exposure",
                               snps = NULL,
                               snp_col="SNP",
                               beta_col="Effect",
                               se_col="StdErr",
                               eaf_col="eaf.exposure",
                               effect_allele_col="EffectAllele",
                               other_allele_col="OtherAllele",
                               pval_col="P.value")


Outcome_DataX <- format_data(OutcomeX,
                             type= "outcome",
                             snps = NULL,
                             snp_col="SNP",
                             beta_col="Effect",
                             se_col="StdErr",
                             eaf_col="eaf.outcome",
                             effect_allele_col="EffectAllele",
                             other_allele_col="OtherAllele",
                             pval_col="P.value")


#-----------------------------------------------------------------------------------------------------

## Harmonise data: 

harmonised.file.preclump <- harmonise_data( exposure_dat = Exposure_DataX, 
                                   outcome_dat = Outcome_DataX)

#-------------------------------------------------------------------------------------------------

## Clump data:
harmonised.file<- clump_data(harmonised.file, clump_kb = 500, clump_r2 = 0.001)

# ?? Record SNPs removed ??
# SNPs_removed_by_clumping<-harmonised.file.preclump$SNP[!test$SNP ] # '!' doesnt work for factors like 'SNP'

#-------------------------------------------------------------------------------------------------------

