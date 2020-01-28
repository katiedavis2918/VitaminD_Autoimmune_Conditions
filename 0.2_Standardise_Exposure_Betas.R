#-------------------------------------------------------------------------------------------------
# "Generic script" 0.2:    STANDARDISE DIRECTION OF BETA'S IN EXPOSURE DATAFRAME
#-------------------------------------------------------------------------------------------------

# Standardise direction of beta's (all positive in exposure dataset)

#------------------
df<-ExposureX

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

