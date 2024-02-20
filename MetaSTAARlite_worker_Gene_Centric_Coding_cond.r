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
output_file_name <- "JHS_LDLR_coding"

###############################
# LDLR coding
###############################
chr <- 19

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

genes <- genes_info

coding_sumstat <- list()
coding_cov <- list()
coding_cov_cond <- list()

known_loci_LDLR <- read.csv("/path_to_the_file/known_loci_LDLR.csv",colClasses=c("integer",
                                                                                 "integer",
                                                                                 "character",
                                                                                 "character"))
gene_name <- "LDLR"
results_temp <- coding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                            genes=genes,known_loci=known_loci_LDLR,
                                            cov_maf_cutoff=0.05,signif.digits=NULL,
                                            QC_label=QC_label,check_qc_label=TRUE,variant_type=variant_type,
                                            Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                            Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
coding_sumstat[[gene_name]] <- results_temp$summary_stat_list
coding_cov[[gene_name]] <- results_temp$GTSinvG_rare_list
coding_cov_cond[[gene_name]] <- results_temp$cov_cond_list

save(coding_sumstat,file=paste0(output_path,output_file_name,"_sumstat.Rdata"),compress = "xz")
save(coding_cov,file=paste0(output_path,output_file_name,"_cov.Rdata"),compress = "xz")
save(coding_cov_cond,file=paste0(output_path,output_file_name,"_cov_cond.Rdata"),compress = "xz")

seqClose(genofile)

