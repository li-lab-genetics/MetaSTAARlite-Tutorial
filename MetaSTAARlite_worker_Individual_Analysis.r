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
## Number of jobs for each chromosome
jobs_num <- read.csv("/path_to_the_file/jobs_num.csv")
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
output_file_name <- "JHS_individual_analysis"
## input array id from batch file
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
chr <- which.max(arrayid <= cumsum(jobs_num$individual_analysis_num))
group.num <- jobs_num$individual_analysis_num[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(jobs_num$individual_analysis_num)[chr-1]
}

start_loc <- (groupid-1)*10e6 + jobs_num$start_loc[chr]
end_loc <- start_loc + 10e6 - 1
end_loc <- min(end_loc,jobs_num$end_loc[chr])

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

individual_analysis_sumstat <- individual_analysis_MetaSTAARlite_worker(chr=chr,start_loc=start_loc,end_loc=end_loc,genofile=genofile,
                                                                        obj_nullmodel=obj_nullmodel,subsegment.size=5e4,
                                                                        QC_label=QC_label,check_qc_label=TRUE,variant_type=variant_type,
                                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

save(individual_analysis_sumstat,file=paste0(output_path,output_file_name,"_sumstat_",arrayid,".Rdata"),compress = "xz")

seqClose(genofile)

