rm(list=ls())
gc()

## load required packages
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(MetaSTAAR)
library(MetaSTAARlite)

###########################################################
#           User Input
###########################################################
## Number of jobs for each chromosome
jobs_num <- read.csv("/path_to_the_file/jobs_num.csv")
## results path
input_path <- "/path_to_the_results_file/"
output_path <- input_path
## results name
individual_results_name <- "individual_analysis"

## alpha level
alpha <- 5E-08

###########################################################
#           Main Function 
###########################################################
Individual_Analysis_Results_Summary_meta(jobs_num=jobs_num,input_path=input_path,output_path=output_path,
                                         individual_results_name=individual_results_name,
                                         alpha=alpha,manhattan_plot=TRUE,QQ_plot=TRUE)
