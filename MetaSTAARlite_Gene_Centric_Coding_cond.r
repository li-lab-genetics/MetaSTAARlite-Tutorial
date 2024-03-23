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
file.dir <- c("/path_to_JHS_coding_cond/",
              "/path_to_MESA_coding_cond/")
file.prefix <- c("JHS_LDLR_coding","MESA_LDLR_coding")
## Sample sizes of participating studies
sample.sizes <- c(2923,4791)

## variant_type
variant_type <- "SNV"
## cov_maf_cutoff
cov_maf_cutoff <- c(0.05,0.05)

## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "LDLR_coding_cond"

genes <- genes_info

results_coding_cond <- c()

sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat.Rdata")
cov.file.path <- paste0(file.dir,file.prefix,"_cov.Rdata")
covcond.file.path <- paste0(file.dir,file.prefix,"_cov_cond.Rdata")
coding_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)
coding_cov_list <- sapply(cov.file.path, function(x) mget(load(x)), simplify = TRUE)
coding_cov_cond_list <- sapply(covcond.file.path, function(x) mget(load(x)), simplify = TRUE)

chr <- 19
gene_name <- "LDLR"
coding_sumstat_gene_list <- lapply(sumstat.file.path, function(x) {
  coding_sumstat_list[[paste0(x,".coding_sumstat")]][[gene_name]]
})
coding_cov_gene_list <- lapply(cov.file.path, function(x) {
  coding_cov_list[[paste0(x,".coding_cov")]][[gene_name]]
})
coding_cov_cond_gene_list <- lapply(covcond.file.path, function(x) {
  coding_cov_cond_list[[paste0(x,".coding_cov_cond")]][[gene_name]]
})
results_cond <- coding_MetaSTAARlite_cond(chr=chr,gene_name=gene_name,genes=genes,
                                          sample.sizes=sample.sizes,coding_sumstat_gene_list=coding_sumstat_gene_list,
                                          coding_cov_gene_list=coding_cov_gene_list,
                                          coding_cov_cond_gene_list=coding_cov_cond_gene_list,
                                          cov_maf_cutoff=cov_maf_cutoff,
                                          rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                          check_qc_label=TRUE,variant_type=variant_type,
                                          Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  
results_coding_cond <- append(results_coding_cond,results_cond)

save(results_coding_cond,file=paste0(output_path,output_file_name,".Rdata"))

