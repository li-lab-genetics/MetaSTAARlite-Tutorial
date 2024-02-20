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
variant_type <- "variant"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/path_to_the_file/Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- FALSE
## Annotation name
Annotation_name <- NULL

## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "JHS_LDLR_individual_analysis"

###############################
# LDLR individual analysis
###############################
chr <- 19

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

known_loci_LDLR <- read.csv("/path_to_the_file/known_loci_LDLR.csv",colClasses=c("integer",
                                                                                 "integer",
                                                                                 "character",
                                                                                 "character"))
start_loc <- 10500001
end_loc <- 11500000
results_temp <- individual_analysis_MetaSTAARlite_worker(chr=chr,start_loc=start_loc,end_loc=end_loc,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                         known_loci=known_loci_LDLR,subsegment.size=5e4,
                                                         QC_label=QC_label,check_qc_label=TRUE,variant_type=variant_type,
                                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
individual_analysis_sumstat <- results_temp$summary_stat
individual_analysis_cov_cond <- results_temp$cov_cond

save(individual_analysis_sumstat,file=paste0(output_path,output_file_name,"_sumstat.Rdata"),compress = "xz")
save(individual_analysis_cov_cond,file=paste0(output_path,output_file_name,"_cov_cond.Rdata"),compress = "xz")

seqClose(genofile)

