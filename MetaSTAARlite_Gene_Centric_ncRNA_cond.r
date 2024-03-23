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
file.dir <- c("/path_to_JHS_ncRNA_cond/",
              "/path_to_MESA_ncRNA_cond/")
file.prefix <- c("JHS_MIR4497_ncRNA","MESA_MIR4497_ncRNA")
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
output_file_name <- "MIR4497_ncRNA_cond"

sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat.Rdata")
cov.file.path <- paste0(file.dir,file.prefix,"_cov.Rdata")
covcond.file.path <- paste0(file.dir,file.prefix,"_cov_cond.Rdata")
ncRNA_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)
ncRNA_cov_list <- sapply(cov.file.path, function(x) mget(load(x)), simplify = TRUE)
ncRNA_cov_cond_list <- sapply(covcond.file.path, function(x) mget(load(x)), simplify = TRUE)

chr <- 12
gene_name <- "MIR4497"
ncRNA_sumstat_gene_list <- lapply(sumstat.file.path, function(x) {
  ncRNA_sumstat_list[[paste0(x,".ncRNA_sumstat")]][[gene_name]]
})
ncRNA_cov_gene_list <- lapply(cov.file.path, function(x) {
  ncRNA_cov_list[[paste0(x,".ncRNA_cov")]][[gene_name]]
})
ncRNA_cov_cond_gene_list <- lapply(covcond.file.path, function(x) {
  ncRNA_cov_cond_list[[paste0(x,".ncRNA_cov_cond")]][[gene_name]]
})

results_ncRNA_cond <- c()
results_ncRNA_cond <- try(ncRNA_MetaSTAARlite_cond(chr=chr,gene_name=gene_name,
                                                   sample.sizes=sample.sizes,ncRNA_sumstat_gene_list=ncRNA_sumstat_gene_list,
                                                   ncRNA_cov_gene_list=ncRNA_cov_gene_list,
                                                   ncRNA_cov_cond_gene_list=ncRNA_cov_cond_gene_list,
                                                   cov_maf_cutoff=cov_maf_cutoff,
                                                   rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                   check_qc_label=TRUE,variant_type=variant_type,
                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))

save(results_ncRNA_cond,file=paste0(output_path,output_file_name,".Rdata"))

