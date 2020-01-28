#-------------------------------------------------------------------------------------------------
# "Generic script" 4:    COMPILE RESULTS 
#-------------------------------------------------------------------------------------------------

# Extract 4 plots to place in a grid:

library(devtools)
install_github("AntonioJBT/episcout")
library(episcout)

# Put 4 plots together:

my_plots <- list(p1[[1]], p2[[1]], p3[[1]], p4[[1]]) # subset 
sapply(my_plots, class)

grid_plots <- episcout::epi_plots_to_grid(my_plots,
                                          label_size = 12,
                                          align = 'v')

#----------------------------------------------------

# **** collect output values from these analyses ****

test.table<-mr_heterogeneity[,c("method","Q_pval")]

tester400<-write.table(MR_all[,c("Method", "Estimate")])

#----------------------------------------------------