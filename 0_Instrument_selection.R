#-------------------------------------------------------------------------------------------------
# "Generic script" 0:    INSTRUMENT SELECTION
#-------------------------------------------------------------------------------------------------

# Select SNPs which pass a given significance level (done before importing into R ?)

# Remove SNPs with horizontal pleiotropy

# Import dataframe with selected SNPs
ExposureX<- read.csv(file="XXX", header = TRUE)

# or add column names if needed
# names(ExposureX) <-c("SNP", "...", ...)

#-------------------------------------------------------------------------------------------------

# Check 'ExposureX' column order correct for applying loop to standardise the direction of the beta's 

col_order <- c("SNP", "EffectAllele", "OtherAllele", "Effect", "StdErr", "P.value", 
               "Position", "Chromosome", "Sample.size", "EAF")

ExposureX <- ExposureX[, col_order]

# -------   Calculate R^2 and F statistics -------------------------------------------------------
# NB: check column order & number of repitions of loop needed 

df<-ExposureX

# Make all EAFs = MAFs for this calculation

f_out<-rep(0,6)

for(i in 1:6){
  
  if (df[i,10] >0.5)                   # column 10 contains EAF values 
    
  {f_out<- (df[i,10] = 1 -df[i,10])}}

#-----------

# Calculate R^2 Statistic:
# df[i,4] = beta
# df[i,10] = MAF
# df[i,5] = se(B) 
# df[i,9] = N (sample size) 

# Formula for R^2:
#f_out[i]<-(2*(beta^2)*(MAF)*(1-(MAF)))/((2*(beta^2)*(MAF)*(1-(MAF)))+((se^2)*(2*N)*(MAF)*(1-(MAF))))

#-----------

R2<-rep(0,6)
for(i in 1:6){
  R2[i]<-(2*df[i,4]^2*df[i,10]*(1-df[i,10]))/((2*df[i,4]^2*df[i,10]*(1-df[i,10]))+(df[i,5]^2*2*df[i,9]*df[i,10]*(1-df[i,10])))
}

# Add R^2 to table:
df$R2<-(R2)

# -----  Compute F.statistic  -----------------

# df[i,11] = R2 
# df[i,9] = N

#f_out[i]<-(R2(N-2))/(1-R2)

F.stat<-rep(0,6)
for(i in 1:6){
  F.stat[i]<-(df[i,11]*(df[i,9] -2))/(1- df[i,11])
}

# Add to table :
df$F.stat<-F.stat

# Save All F stats and R2 values:
Fstat_all_SNPs<-df

# Remove rows with F stats <10:
Instrument_SNPs<- df[df$F.stat>=10,]
# Save names of SNPs to include
Instrument_SNPs_names<- data.frame(Instrument_SNPs$SNP)
names(Instrument_SNPs_names) <-("SNP")
# Merge with 'ExposureX' (So that Allele Frequencies are EAFs and not MAFs as used for previous calculations)

#SAVE FINAL TABLE
ExposureX<-merge(ExposureX,Instrument_SNPs_names, by= "SNP")

#record names of removed SNPs?

#----------------------------------------

# Find proxy SNPs for those not found in the outcome dataset (Manually?) 

# import outcome data

OutcomeX<-read.csv(file = "XXX", header = TRUE)

# or add column names if needed
# names(OutcomeX) <-c("SNP", "...", ...)
#----------------------------------------

