rm(list=ls())
gc()

## load required packages
library(STAAR)
library(STAARpipeline)
library(MetaSTAAR)
library(MetaSTAARlite)

###########################################################
#           User Input
###########################################################
## Directories of the study-specific summary statistics file folders
file.dir <- c("/path_to_JHS_individual_analysis/",
              "/path_to_MESA_individual_analysis/")
file.prefix <- c("JHS_individual_analysis","MESA_individual_analysis")
## Sample sizes of participating studies
sample.sizes <- c(2923,4791)

## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "individual_analysis"
## input array id from batch file
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat_",arrayid,".Rdata")
individual_analysis_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)

results_individual_analysis <- individual_analysis_MetaSTAARlite(sample.sizes=sample.sizes,
                                                                 sumstat.list=individual_analysis_sumstat_list,
                                                                 mac_cutoff=20,check_qc_label=TRUE)

save(results_individual_analysis,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

