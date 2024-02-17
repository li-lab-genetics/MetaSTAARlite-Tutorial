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
## aGDS directory
agds_dir <- get(load("/path_to_the_file/agds_dir.Rdata"))
## Null model
obj_nullmodel <- get(load("/path_to_the_file/obj_nullmodel.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "SNV"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/path_to_the_file/Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "JHS_noncoding"
## input array id from batch file
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
## gene number in job
gene_num_in_array <- 50
group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
sum(group.num.allchr)

chr <- which.max(arrayid <= cumsum(group.num.allchr))
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

## exclude large noncoding masks
jobid_exclude <- c(21,39,44,45,46,53,55,83,88,103,114,127,135,150,154,155,163,164,166,180,189,195,200,233,280,285,295,313,318,319,324,327,363,44,45,54)
sub_seq_id_exclude <- c(1009,1929,182,214,270,626,741,894,83,51,611,385,771,493,671,702,238,297,388,352,13,303,600,170,554,207,724,755,1048,319,324,44,411,195,236,677)

for(i in 1:length(jobid_exclude))
{
  if(arrayid==jobid_exclude[i])
  {
    sub_seq_id <- setdiff(sub_seq_id,sub_seq_id_exclude[i])
  }
}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

noncoding_sumstat <- list()
noncoding_cov <- list()
for(kk in sub_seq_id)
{
  print(kk)
  gene_name <- genes_info_chr[kk,1]
  results_temp <- noncoding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                 cov_maf_cutoff=0.05,signif.digits=NULL,
                                                 QC_label=QC_label,check_qc_label=TRUE,variant_type=variant_type,
                                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  noncoding_sumstat[[gene_name]] <- results_temp$summary_stat_list
  noncoding_cov[[gene_name]] <- results_temp$GTSinvG_rare_list
}

save(noncoding_sumstat,file=paste0(output_path,output_file_name,"_sumstat_",arrayid,".Rdata"),compress = "xz")
save(noncoding_cov,file=paste0(output_path,output_file_name,"_cov_",arrayid,".Rdata"),compress = "xz")

seqClose(genofile)

