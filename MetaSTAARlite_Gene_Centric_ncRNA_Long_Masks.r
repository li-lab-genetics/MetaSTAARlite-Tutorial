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
file.dir <- c("/path_to_JHS_ncRNA/",
              "/path_to_MESA_ncRNA/")
file.prefix <- c("JHS_ncRNA","MESA_ncRNA")
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
output_file_name <- "ncRNA"

###########################################################
#           Main Function 
###########################################################
## analyze large ncRNA masks
arrayid <- c(117,218,220,220,221,156,219)
sub_seq_id <- c(53,19,208,274,311,41,103)

region_spec <- data.frame(arrayid,sub_seq_id)

gene_num_in_array <- 100
group.num.allchr <- ceiling(table(ncRNA_gene[,1])/gene_num_in_array)
sum(group.num.allchr)

results_ncRNA <- c()

sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat_",223,".Rdata")
cov.file.path <- paste0(file.dir,file.prefix,"_cov_",223,".Rdata")
ncRNA_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)
ncRNA_cov_list <- sapply(cov.file.path, function(x) mget(load(x)), simplify = TRUE)

for(kk in 1:dim(region_spec)[1])
{
  arrayid <- region_spec$arrayid[kk]
  sub_seq_id <- region_spec$sub_seq_id[kk]
  
  chr <- which.max(arrayid <= cumsum(group.num.allchr))
  ncRNA_gene_chr <- ncRNA_gene[ncRNA_gene[,1]==chr,]
  gene_name <- ncRNA_gene_chr[sub_seq_id,2]
  ncRNA_sumstat_gene_list <- lapply(sumstat.file.path, function(x) {
    ncRNA_sumstat_list[[paste0(x,".ncRNA_sumstat")]][[gene_name]]
  })
  ncRNA_cov_gene_list <- lapply(cov.file.path, function(x) {
    ncRNA_cov_list[[paste0(x,".ncRNA_cov")]][[gene_name]]
  })
  results <- c()
  results <- try(ncRNA_MetaSTAARlite(chr=chr,gene_name=gene_name,
                                     sample.sizes=sample.sizes,ncRNA_sumstat_gene_list=ncRNA_sumstat_gene_list,
                                     ncRNA_cov_gene_list=ncRNA_cov_gene_list,
                                     cov_maf_cutoff=cov_maf_cutoff,
                                     rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                     check_qc_label=TRUE,variant_type=variant_type,
                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
  results_ncRNA <- rbind(results_ncRNA,results)
}

save(results_ncRNA,file=paste0(output_path,output_file_name,"_",223,".Rdata"))

