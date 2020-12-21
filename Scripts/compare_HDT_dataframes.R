# Script to verify that Bari_HDT_datalist.txt is the same as HDT_raw-task-data.txt
library(arsenal)

# set basedir to navigate to HDT directory (modify if on different machine)
basedir <- "~/OneDrive - The Pennsylvania State University/b-childfoodlab_Shared/Inactive_Studies/DMK_Study/AB_reprocess/HDT/HDT-intake-git"

# Set directory to run script from 
setwd(file.path(basedir,"Data/GeneratedDatabases"))

# Read in files
HDT_raw_task <-  read.delim("HDT_raw-task-data.txt", stringsAsFactors = FALSE)
Bari_datalist <-  read.delim("Bari_HDT_datalist.txt", stringsAsFactors = FALSE)

comparedf(HDT_raw_task, Bari_datalist)

