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
## file directory of GDS file (genotype and annotation data) 
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
output_file_name <- "JHS_custom"

###########################################################
#           Main Function 
###########################################################
chr <- 19

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

position <- as.numeric(seqGetData(genofile, "position"))
REF <- as.character(seqGetData(genofile, "$ref"))
ALT <- as.character(seqGetData(genofile, "$alt"))
agds_variant_list <- data.frame(CHR=chr,POS=position,REF=REF,ALT=ALT)

## Mask file
custom_mask <- read.csv("/path_to_the_file/custom_mask_LDLR.csv",colClasses=c("integer",
                                                                              "integer",
                                                                              "character",
                                                                              "character",
                                                                              "character"))
mask_names <- unique(custom_mask$MaskName)

custom_sumstat <- list()
custom_cov <- list()
for(mask_name in mask_names)
{
  print(mask_name)
  variant_list <- custom_mask[custom_mask$mask == mask_name,]
  results_temp <- custom_MetaSTAARlite_worker(chr=chr,variant_list=variant_list,agds_variant_list=agds_variant_list,
                                              genofile=genofile,obj_nullmodel=obj_nullmodel,
                                              cov_maf_cutoff=0.05,signif.digits=NULL,
                                              QC_label=QC_label,check_qc_label=TRUE,variant_type=variant_type,
                                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  custom_sumstat[[mask_name]] <- results_temp$summary_stat
  custom_cov[[mask_name]] <- results_temp$GTSinvG_rare
}

save(custom_sumstat,file=paste0(output_path,output_file_name,"_sumstat.Rdata"),compress = "xz")
save(custom_cov,file=paste0(output_path,output_file_name,"_cov.Rdata"),compress = "xz")

seqClose(genofile)

