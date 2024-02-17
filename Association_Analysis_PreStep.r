rm(list=ls())
gc()

###########################################################
#           User Input
###########################################################
## file directory of aGDS file (genotype and annotation data) 
dir.geno <- "/path_to_the_aGDS_file/"
## file name of aGDS, seperate by chr number 
agds_file_name_1 <- "freeze.5.chr"
agds_file_name_2 <- ".pass_and_fail.gtonly.minDP0.gds"
## channel name of the QC label in the GDS/aGDS file
QC_label <- "annotation/filter"
## file directory for the output files
output_path <- "/path_to_the_output_file/" 
## annotation name. The first eight names are used to define masks in gene-centric analysis, do not change them!! 
## The others are the annotation you want to use in the STAAR procedure, and they are flexible to change.
name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category",
          "MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF",
          "aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
          "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
## channel name of the annotations. Make sure they are matched with the name, especially for the first eight one!! 
dir <- c("/rsid","/genecode_comprehensive_category","/genecode_comprehensive_info",
         "/genecode_comprehensive_exonic_category","/metasvm_pred",
         "/genehancer","/cage_tc","/rdhs","/cadd_phred","/linsight","/fathmm_xf",
         "/apc_epigenetics_active","/apc_epigenetics_repressed","/apc_epigenetics_transcription",
         "/apc_conservation","/apc_local_nucleotide_diversity","/apc_mappability",
         "/apc_transcription_factor","/apc_protein_function")

###########################################################
#           Main Function 
###########################################################
## aGDS directory
agds_dir <- paste0(dir.geno,agds_file_name_1,seq(1,22),agds_file_name_2) 
save(agds_dir,file=paste0(output_path,"agds_dir.Rdata",sep=""))

## Annotation name catalog (alternatively, can skip this part by providing Annotation_name_catalog.csv with the same information)
Annotation_name_catalog <- data.frame(name=name,dir=dir)
save(Annotation_name_catalog,file=paste0(output_path,"Annotation_name_catalog.Rdata",sep=""))

