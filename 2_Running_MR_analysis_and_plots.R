#-------------------------------------------------------------------------------------------------
# "Generic script" 2:    DATA FORMATTING 
#-------------------------------------------------------------------------------------------------

# Output from '1_Data_Formatting' script: "harmonised.file"

#=========================MR analysis: Exposure- VEGF & Outcome- VITAMIN D =======================

#INPUT
harmonised.file

#Define exposure and outcome
exp = ("VitaminD")
out = ("MIP1B")
#=================================================================================================

#Perform MR Analysis
mr_results<-mr(harmonised.file, method_list = c("mr_ivw","mr_egger_regression", "mr_weighted_median"))

# Create MRInput to do these with MendelianRandomisation Package (to get CI's)

MRInputObject<-mr_input(bx=harmonised.file$beta.exposure,
                        bxse = harmonised.file$se.exposure,
                        by= harmonised.file$beta.outcome,
                        byse= harmonised.file$se.outcome,
                        exposure = exp,
                        outcome = out,
                        snps = harmonised.file$SNP)

#Run Analyses (This will, by default, perform random effects IVW as there are >3 genetic variants)
MR_all<- mr_allmethods(MRInputObject, method="all")


# To perform fixed effects IVW
IVW<- MendelianRandomization::mr_ivw(MRInputObject,
                                     model = "fixed",
                                     robust = FALSE,
                                     penalized = FALSE,
                                     weights = "simple",
                                     distribution = "normal",
                                     alpha = 0.05)

# To get I^2 Statistic
Egger<- mr_egger(MRInputObject,
                 robust=FALSE,
                 penalized =FALSE,
                 distribution ="normal",
                 alpha= 0.05 )


## ----Sensitivity analysis---------------------------

#Obtain heterogeneity statistics 
mr_heterogeneity<- mr_heterogeneity(harmonised.file)

# Test for directional pleiotropy
mr_pleiotropy_test<- mr_pleiotropy_test(harmonised.file)  

#Obtain MR estimates for each of the selected IVs
res_single<-mr_singlesnp(harmonised.file, 
                         all_method = c("mr_ivw","mr_egger_regression", "mr_weighted_median"))

#Obtain MR estimates excluding one IV at a time 
res_loo<-mr_leaveoneout(harmonised.file)

##--------CREATING PLOTS --------------------------

#Create scatter plots of IV-outcome associations against IV-exposure
p1<- mr_scatter_plot(mr_results, harmonised.file)

#Create forest plot of causal effects calculated from each IV
p2<-mr_forest_plot(res_single)

#Create plot of causal effects derived from removing IVs sequentially
p3<-mr_leaveoneout_plot(res_loo)

#Create funnel plot of the reciprocal of the SE vs the MR causal estimate 
p4<- mr_funnel_plot(res_single)

#----------------------------------------------------

