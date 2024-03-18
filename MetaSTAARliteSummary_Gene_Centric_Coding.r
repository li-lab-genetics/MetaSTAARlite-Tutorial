rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(MetaSTAAR)
library(MetaSTAARlite)

###########################################################
#           User Input
###########################################################
## results path
input_path <- "/path_to_the_results_file/"
output_path <- input_path
## number of jobs
gene_centric_coding_jobs_num <- 381
## results name
gene_centric_results_name <- "coding"

## alpha level
alpha <- 5E-07

###########################################################
#           Main Function
###########################################################
Gene_Centric_Coding_Results_Summary_meta(gene_centric_coding_jobs_num=gene_centric_coding_jobs_num,
                                         input_path=input_path,output_path=output_path,
                                         gene_centric_results_name=gene_centric_results_name,
                                         alpha=alpha,manhattan_plot=TRUE,QQ_plot=TRUE)
