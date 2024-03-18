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
gene_centric_noncoding_jobs_num <- 387
## results name
gene_centric_results_name <- "noncoding"

## alpha level
alpha <- 3.57E-07

## ncRNA
ncRNA_jobs_num <- 223
ncRNA_input_path <- "/path_to_the_results_file/"
ncRNA_output_path <- ncRNA_input_path
ncRNA_results_name <- "ncRNA"

###########################################################
#           Main Function
###########################################################
## gene info
Gene_Centric_Noncoding_Results_Summary_meta(gene_centric_noncoding_jobs_num=gene_centric_noncoding_jobs_num,
                                            input_path=input_path,output_path=output_path,
                                            gene_centric_results_name=gene_centric_results_name,
                                            ncRNA_jobs_num=ncRNA_jobs_num,ncRNA_input_path=ncRNA_input_path,
                                            ncRNA_output_path=ncRNA_output_path,ncRNA_results_name=ncRNA_results_name,
                                            alpha=alpha,manhattan_plot=TRUE,QQ_plot=TRUE)
