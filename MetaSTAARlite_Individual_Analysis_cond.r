rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(MetaSTAAR)
library(MetaSTAARlite)

###########################################################
#           User Input
###########################################################
## Directories of the study-specific summary statistics file folders
file.dir <- c("/path_to_JHS_individual_analysis_cond/",
              "/path_to_MESA_individual_analysis_cond/")
file.prefix <- c("JHS_LDLR_individual_analysis","MESA_LDLR_individual_analysis")
## Sample sizes of participating studies
sample.sizes <- c(2923,4791)

individual_results <- read.csv("/path_to_the_file/individual_results.csv")
individual_results <- individual_results[individual_results$CHR == 19 & 
                                           individual_results$POS >= 10500001 & 
                                           individual_results$POS <= 11500000, c(1:4)]

## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "LDLR_individual_analysis_cond"

###########################################################
#           Main Function 
###########################################################
sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat.Rdata")
covcond.file.path <- paste0(file.dir,file.prefix,"_cov_cond.Rdata")
individual_analysis_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)
individual_analysis_cov_cond <- sapply(covcond.file.path, function(x) mget(load(x)), simplify = TRUE)

results_individual_analysis_cond <- individual_analysis_MetaSTAARlite_cond(individual_results=individual_results,sample.sizes=sample.sizes,
                                                                           sumstat.list=individual_analysis_sumstat_list,
                                                                           covcond.list=individual_analysis_cov_cond,
                                                                           mac_cutoff=20,check_qc_label=TRUE)

row.names(results_individual_analysis_cond) <- NULL
save(results_individual_analysis_cond,file=paste0(output_path,output_file_name,".Rdata"))

