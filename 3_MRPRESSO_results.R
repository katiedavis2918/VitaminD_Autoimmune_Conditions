#-------------------------------------------------------------------------------------------------
# "Generic script" 2.1:    MR-PRESSO 
#-------------------------------------------------------------------------------------------------

library(MRPRESSO)

#-------------------------------------------------------------------------------------------------

#Inputs 
beta.exposure<-(harmonised.file$beta.exposure)
se.exposure<-(harmonised.file$se.exposure)
beta.outcome<-(harmonised.file$beta.outcome)
se.outcome<-(harmonised.file$se.outcome)
SNPcol<- (harmonised.file$SNP)

# ================================================================================================

MRPRESSO_input = data.frame(beta_exposure=beta.exposure, 
                            se_exposure=se.exposure, 
                            beta_outcome=beta.outcome, 
                            se_outcome=se.outcome, 
                            row.names = SNPcol)

#-------------------------------------------------------------------------------------------------

MRPRESSO_ResultsX <- mr_presso(BetaOutcome = "beta_outcome", 
                              BetaExposure = "beta_exposure", 
                              SdOutcome = "se_outcome", 
                              SdExposure = "se_exposure", 
                              OUTLIERtest = TRUE, 
                              DISTORTIONtest = TRUE, 
                              data = MRPRESSO_input, 
                              NbDistribution = 1000,  
                              SignifThreshold = 0.05)
#-------------------------------------------------------------------------------------------------