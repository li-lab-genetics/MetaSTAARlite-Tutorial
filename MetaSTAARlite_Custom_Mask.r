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
file.dir <- c("/path_to_JHS_custom/",
              "/path_to_MESA_custom/")
file.prefix <- c("JHS_custom","MESA_custom")
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
output_file_name <- "custom"

###########################################################
#           Main Function 
###########################################################
chr <- 19

## Mask file
custom_mask <- read.csv("/path_to_the_file/custom_mask_LDLR.csv",colClasses=c("integer",
                                                                              "integer",
                                                                              "character",
                                                                              "character",
                                                                              "character"))
mask_names <- unique(custom_mask$MaskName)

results_custom <- c()

sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat.Rdata")
cov.file.path <- paste0(file.dir,file.prefix,"_cov.Rdata")
custom_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)
custom_cov_list <- sapply(cov.file.path, function(x) mget(load(x)), simplify = TRUE)

for(mask_name in mask_names)
{
  print(mask_name)
  custom_sumstat_mask_list <- lapply(sumstat.file.path, function(x) {
    custom_sumstat_list[[paste0(x,".custom_sumstat")]][[mask_name]]
  })
  custom_cov_mask_list <- lapply(cov.file.path, function(x) {
    custom_cov_list[[paste0(x,".custom_cov")]][[mask_name]]
  })
  results <- custom_MetaSTAARlite(chr=chr,mask_name=mask_name,
                                  sample.sizes=sample.sizes,custom_sumstat_mask_list=custom_sumstat_mask_list,
                                  custom_cov_mask_list=custom_cov_mask_list,
                                  cov_maf_cutoff=cov_maf_cutoff,
                                  rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                  check_qc_label=TRUE,variant_type=variant_type,
                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  
  results_custom <- rbind(results_custom,results)
}

save(results_custom,file=paste0(output_path,output_file_name,".Rdata"))

