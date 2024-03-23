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
file.dir <- c("/path_to_JHS_coding/",
              "/path_to_MESA_coding/")
file.prefix <- c("JHS_coding","MESA_coding")
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
output_file_name <- "coding"
## input array id from batch file
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
## gene number in job
gene_num_in_array <- 50
group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
sum(group.num.allchr)

chr <- as.integer(which.max(arrayid <= cumsum(group.num.allchr)))
group.num <- group.num.allchr[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(group.num.allchr)[chr-1]
}

genes_info_chr <- genes_info[genes_info[,2]==chr,]
sub_seq_num <- dim(genes_info_chr)[1]

if(groupid < group.num)
{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
}else
{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
}

## exclude large coding masks
if(arrayid==57)
{
  sub_seq_id <- setdiff(sub_seq_id,840)
}

if(arrayid==112)
{
  sub_seq_id <- setdiff(sub_seq_id,c(543,544))
}

if(arrayid==113)
{
  sub_seq_id <- setdiff(sub_seq_id,c(575,576,577,578,579,580,582))
}

genes <- genes_info

results_coding <- c()

sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat_",arrayid,".Rdata")
cov.file.path <- paste0(file.dir,file.prefix,"_cov_",arrayid,".Rdata")
coding_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)
coding_cov_list <- sapply(cov.file.path, function(x) mget(load(x)), simplify = TRUE)

for(kk in sub_seq_id)
{
  print(kk)
  gene_name <- genes_info_chr[kk,1]
  coding_sumstat_gene_list <- lapply(sumstat.file.path, function(x) {
    coding_sumstat_list[[paste0(x,".coding_sumstat")]][[gene_name]]
  })
  coding_cov_gene_list <- lapply(cov.file.path, function(x) {
    coding_cov_list[[paste0(x,".coding_cov")]][[gene_name]]
  })
  results <- coding_MetaSTAARlite(chr=chr,gene_name=gene_name,genes=genes,
                                  sample.sizes=sample.sizes,coding_sumstat_gene_list=coding_sumstat_gene_list,
                                  coding_cov_gene_list=coding_cov_gene_list,
                                  cov_maf_cutoff=cov_maf_cutoff,
                                  rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                  check_qc_label=TRUE,variant_type=variant_type,
                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  results_coding <- append(results_coding,results)
}

save(results_coding,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

