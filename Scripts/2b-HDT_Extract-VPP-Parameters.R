#install.packages("rjags")
library("rjags")

# This script was written by Bari Fuchs and Nicole Roberts to process stanfit object (output from fitting Stan model in 2a_HDT_VPPmodel.R)

# set basedir to navigate to HDT directory (modify if on different machine)
basedir <- "~/OneDrive - The Pennsylvania State University/b-childfoodlab_Shared/Inactive_Studies/DMK_Study/AB_reprocess/HDT/HDT-intake-git"

# Set directory to run script from 
setwd(file.path(basedir,"Scripts/"))

# load VPP model Rdata file
load(file.path(basedir,"Data/GeneratedDatabases/VPPKids2019-10-11.Rdata"))

# Generate results table with individual parameters, posterior stats
source("postcalc.R")
codaSamples = mcmc.list(lapply(1:ncol(fit), function(x) {
  mcmc(as.array(fit)[,x, ])
}))

resulttable <- zcalc(codaSamples)


### Extract individual-level mean estimates for VPP model parameters from resulttable ###

# Make new dataframe with columns for each parameter
HDT_VPP_means <- data.frame(matrix(nrow = 70, ncol = 9))
colnames(HDT_VPP_means) <- c("subjID", "A_mean", "alpha_mean", "lambda_mean", "epP_mean", "epN_mean", "K_mean", "w_mean", "cons_mean")

# Extract means for each variable
HDT_VPP_means$A_mean <- resulttable[17:86, "mean"]
HDT_VPP_means$alpha_mean <- resulttable[87:156, "mean"]
HDT_VPP_means$cons_mean <- resulttable[157:226, "mean"]
HDT_VPP_means$lambda_mean <- resulttable[227:296, "mean"]
HDT_VPP_means$epP_mean <- resulttable[297:366, "mean"]
HDT_VPP_means$epN_mean <- resulttable[367:436, "mean"]
HDT_VPP_means$K_mean <- resulttable[437:506, "mean"]
HDT_VPP_means$w_mean <- resulttable[507:576, "mean"]

# Add subject IDs to table
setwd(file.path(basedir,"Data/GeneratedDatabases"))
pardat <- read.csv("HDT_participant_data.csv", stringsAsFactors = FALSE)

### Export tables ###
# Set directory to export file to
setwd(file.path(basedir,"Data/GeneratedDatabases"))

# write to csv files
write.csv(resulttable,'HDT_VPP_resulttable.csv')
write.csv(HDT_VPP_means,'HDT_VPP_means.csv')

