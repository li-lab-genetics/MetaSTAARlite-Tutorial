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
output_file_name <- "JHS_ncRNA"

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

ncRNA_sumstat <- list()
ncRNA_cov <- list()
for(kk in 1:dim(region_spec)[1])
{
  arrayid <- region_spec$arrayid[kk]
  sub_seq_id <- region_spec$sub_seq_id[kk]
  
  chr <- which.max(arrayid <= cumsum(group.num.allchr))
  ncRNA_gene_chr <- ncRNA_gene[ncRNA_gene[,1]==chr,]
  
  ## aGDS file
  agds.path <- agds_dir[chr]
  genofile <- seqOpen(agds.path)
  
  gene_name <- ncRNA_gene_chr[sub_seq_id,2]
  results_temp <- ncRNA_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                             cov_maf_cutoff=0.05,signif.digits=NULL,
                                             QC_label=QC_label,check_qc_label=TRUE,variant_type=variant_type,
                                             Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                             Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  ncRNA_sumstat[[gene_name]] <- results_temp$summary_stat
  ncRNA_cov[[gene_name]] <- results_temp$GTSinvG_rare
  
  seqClose(genofile)
}

save(ncRNA_sumstat,file=paste0(output_path,output_file_name,"_sumstat_",223,".Rdata"),compress = "xz")
save(ncRNA_cov,file=paste0(output_path,output_file_name,"_cov_",223,".Rdata"),compress = "xz")

