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
file.dir <- c("/path_to_JHS_noncoding/",
              "/path_to_MESA_noncoding/")
file.prefix <- c("JHS_noncoding","MESA_noncoding")
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
output_file_name <- "noncoding"
## input array id from batch file
arrayid_longmask <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
## gene number in job
gene_num_in_array <- 50
group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
sum(group.num.allchr)

## analyze large noncoding masks
arrayid <- c(21,39,44,45,46,53,55,83,88,103,114,127,135,150,154,155,163,164,166,180,189,195,200,233,280,285,295,313,318,319,324,327,363,44,45,54)
sub_seq_id <- c(1009,1929,182,214,270,626,741,894,83,51,611,385,771,493,671,702,238,297,388,352,13,303,600,170,554,207,724,755,1048,319,324,44,411,195,236,677)

region_spec <- data.frame(arrayid,sub_seq_id) 
sub_seq_id <- ((arrayid_longmask-1)*5+1):min(arrayid_longmask*5,length(arrayid))

results_noncoding <- c()

sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat_",arrayid_longmask+379,".Rdata")
cov.file.path <- paste0(file.dir,file.prefix,"_cov_",arrayid_longmask+379,".Rdata")
noncoding_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)
noncoding_cov_list <- sapply(cov.file.path, function(x) mget(load(x)), simplify = TRUE)

for(kk in sub_seq_id)
{
  print(kk)
  arrayid <- region_spec$arrayid[kk]
  sub_id <- region_spec$sub_seq_id[kk]
  
  chr <- which.max(arrayid <= cumsum(group.num.allchr))
  genes_info_chr <- genes_info[genes_info[,2]==chr,]
  gene_name <- genes_info_chr[sub_id,1]
  noncoding_sumstat_gene_list <- lapply(sumstat.file.path, function(x) {
    noncoding_sumstat_list[[paste0(x,".noncoding_sumstat")]][[gene_name]]
  })
  noncoding_cov_gene_list <- lapply(cov.file.path, function(x) {
    noncoding_cov_list[[paste0(x,".noncoding_cov")]][[gene_name]]
  })
  results <- noncoding_MetaSTAARlite(chr=chr,gene_name=gene_name,
                                     sample.sizes=sample.sizes,noncoding_sumstat_gene_list=noncoding_sumstat_gene_list,
                                     noncoding_cov_gene_list=noncoding_cov_gene_list,
                                     cov_maf_cutoff=cov_maf_cutoff,
                                     rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                     check_qc_label=TRUE,variant_type=variant_type,
                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  results_noncoding <- append(results_noncoding,results)
}

save(results_noncoding,file=paste0(output_path,output_file_name,"_",arrayid_longmask+379,".Rdata"))

