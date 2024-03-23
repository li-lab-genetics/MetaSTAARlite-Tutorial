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
file.dir <- c("/path_to_JHS_noncoding_cond/",
              "/path_to_JHS_noncoding_cond/")
file.prefix <- c("JHS_LDLR_noncoding","MESA_LDLR_noncoding")
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
output_file_name <- "LDLR_noncoding_cond"

results_noncoding_cond <- c()

sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat.Rdata")
cov.file.path <- paste0(file.dir,file.prefix,"_cov.Rdata")
covcond.file.path <- paste0(file.dir,file.prefix,"_cov_cond.Rdata")
noncoding_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)
noncoding_cov_list <- sapply(cov.file.path, function(x) mget(load(x)), simplify = TRUE)
noncoding_cov_cond_list <- sapply(covcond.file.path, function(x) mget(load(x)), simplify = TRUE)

chr <- 19
gene_name <- "LDLR"
noncoding_sumstat_gene_list <- lapply(sumstat.file.path, function(x) {
  noncoding_sumstat_list[[paste0(x,".noncoding_sumstat")]][[gene_name]]
})
noncoding_cov_gene_list <- lapply(cov.file.path, function(x) {
  noncoding_cov_list[[paste0(x,".noncoding_cov")]][[gene_name]]
})
noncoding_cov_cond_gene_list <- lapply(covcond.file.path, function(x) {
  noncoding_cov_cond_list[[paste0(x,".noncoding_cov_cond")]][[gene_name]]
})
results_cond <- noncoding_MetaSTAARlite_cond(chr=chr,gene_name=gene_name,
                                             sample.sizes=sample.sizes,noncoding_sumstat_gene_list=noncoding_sumstat_gene_list,
                                             noncoding_cov_gene_list=noncoding_cov_gene_list,
                                             noncoding_cov_cond_gene_list=noncoding_cov_cond_gene_list,
                                             cov_maf_cutoff=cov_maf_cutoff,
                                             rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                             check_qc_label=TRUE,variant_type=variant_type,
                                             Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

results_noncoding_cond <- append(results_noncoding_cond,results_cond)

save(results_noncoding_cond,file=paste0(output_path,output_file_name,".Rdata"))

