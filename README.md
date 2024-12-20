# MetaSTAARlite-Tutorial
This is a tutorial for (1) automatically functionally annotating the variants of each participating whole-genome/whole-exome sequencing (WGS/WES) study and integrating the functional annotations with the genotype data using **FAVORannotator**, (2) generating study-specific variant summary statistics of each participating WGS/WES study using **MetaSTAARlite Worker**, and (3) performing association meta-analysis of WGS/WES studies using **MetaSTAARlite**. The software prerequisites, dependencies and installation can be found in the <a href="https://github.com/li-lab-genetics/MetaSTAARlite">**MetaSTAARlite**</a> package.

## Pre-step of association meta-analysis using MetaSTAARlite (same as <a href="https://github.com/xihaoli/STAARpipeline-Tutorial">STAARpipeline</a>)
### Generate study-specific Genomic Data Structure (GDS) file
R/Bioconductor package **SeqArray** provides functions to convert the genotype data (in VCF/BCF/PLINK BED/SNPRelate format) to SeqArray GDS format. For more details on usage, please see the R/Bioconductor package <a href="https://bioconductor.org/packages/release/bioc/html/SeqArray.html">**SeqArray**</a> [<a href="https://bioconductor.org/packages/release/bioc/manuals/SeqArray/man/SeqArray.pdf">manual</a>]. A wrapper for the `seqVCF2GDS`/`seqBCF2GDS` function in the SeqArray package can be found <a href="convertVCF2GDS.R">**here**</a> (**Credit: Michael R. Brown and Jennifer A. Brody**).

R package **gds2bgen** provides functions to convert the genotype data (in BGEN format) to SeqArray GDS format. For more details on usage, please see the R package <a href="https://github.com/zhengxwen/gds2bgen">**gds2bgen**</a>. An example for the `seqBGEN2GDS` function in the gds2bgen package can be found <a href="https://github.com/zhengxwen/gds2bgen#examples">**here**</a> (**Credit: Xiuwen Zheng**).

Note 1: As a file integrity check, it is expected that variant in the GDS file can be **uniquely identified** based on its **CHR-POS-REF-ALT** combination. That is, there shouldn't be two variants in the GDS file with identical CHR-POS-REF-ALT records. It is also expected that the physical positions of variants in the GDS file (of each chromosome) should be sorted in **ascending order**.

Note 2: After the (study-specific) GDS file is generated, there is supposed to be a channel in the GDS file (default is `annotation/filter`) where all variants passing the quality control (QC) should be labeled as `"PASS"`. If there is no such channel for a given post-QC GDS file (where all variants in the GDS file are pass variants), one can create a new channel in the GDS file by setting the value of all variants as `"PASS"`. An example script can be found <a href="https://github.com/xihaoli/STAARpipeline-Tutorial/blob/main/Add_QC_label.R">**here**</a>. Then, in all scripts of MetaSTAARlite, `QC_label <- "annotation/filter"` should be updated to `QC_label <- "annotation/info/QC_label"`.

### Generate study-specific annotated GDS (aGDS) file using FAVORannotator
#### Prerequisites:
**FAVORannotator** (CSV version 1.0.0) depends on the **xsv software** and the **FAVOR database** in CSV format. Please install the <a href="https://github.com/BurntSushi/xsv">**xsv software**</a> and download the **FAVOR essential database CSV files** from <a href="http://favor.genohub.org">**FAVOR website**</a> (under the "FAVORannotator" tab's top panel, 31.2 GB for chr1 CSV) or <a href="https://doi.org/10.7910/DVN/1VGTJI">**Harvard Dataverse**</a> before using **FAVORannotator** (CSV version 1.0.0).
#### Step 0: Install xsv
The following steps are for the widely used operating system (Ubuntu) on a virtual machine.

1. Install Rust and Cargo:
 - ```$ curl https://sh.rustup.rs -sSf | sh```
2. Source the environment: 
 - ```$ source $HOME/.cargo/env``` 
3. Install xsv using Cargo:
 - ```$ cargo install xsv```
#### Step 1: Generate the variants list to be annotated
##### Script: <a href="https://github.com/xihaoli/STAARpipeline-Tutorial/tree/main/FAVORannotator_csv/Varinfo_gds.R">**Varinfo_gds.R**</a>
##### Input: GDS files of each chromosome and the FAVOR database information <a href="https://github.com/xihaoli/STAARpipeline-Tutorial/tree/main/FAVORannotator_csv/FAVORdatabase_chrsplit.csv">**FAVORdatabase_chrsplit.csv**</a>. For more details, please see the R script.
##### Output: CSV files of the variants list. For each chromosome, the number of CSV files is listed in <a href="https://github.com/xihaoli/STAARpipeline-Tutorial/tree/main/FAVORannotator_csv/FAVORdatabase_chrsplit.csv">**FAVORdatabase_chrsplit.csv**</a>.
Note: The physical positions of variants in the GDS file (of each chromosome) should be sorted in ascending order.

#### Step 2: Annotate the variants using the FAVOR database through xsv software
##### Script: <a href="https://github.com/xihaoli/STAARpipeline-Tutorial/tree/main/FAVORannotator_csv/Annotate.R">**Annotate.R**</a>
##### Input: CSV files of the variants list to be annotated, the FAVOR database information <a href="https://github.com/xihaoli/STAARpipeline-Tutorial/tree/main/FAVORannotator_csv/FAVORdatabase_chrsplit.csv">**FAVORdatabase_chrsplit.csv**</a>,
the FAVOR database, and the directory xsv software. For more details, please see the R script.
##### Output: CSV files of the annotated variants list. 
* `Anno_chrXX.csv`: a CSV file containing annotated variants list of chromosome XX. <br>
* `Anno_chrXX_STAARpipeline.csv`: a CSV file containing the variants list with annotations required for STAARpipeline of chromosome XX. 
The annotations in this file is a subset of `Anno_chrXX.csv`. <br>

#### Step 3: Generate the annotated GDS (aGDS) file
##### Script: <a href="https://github.com/xihaoli/STAARpipeline-Tutorial/tree/main/FAVORannotator_csv/gds2agds.R">**gds2agds.R**</a>
##### Input: GDS files and the CSV files of annotated variants list (`Anno_chrXX.csv` or `Anno_chrXX_STAARpipeline.csv`). For more details, please see the R script.
##### Output: aGDS files including both the genotype and annotation information.
Note: FAVORannotator also supports the database in SQL format. Please see the <a href="https://github.com/zhouhufeng/FAVORannotator">**FAVORannotator** tutorial</a> for detailed usage of **FAVORannotator** (SQL version).

### Generate study-specific sparse Genetic Relatedness Matrix (GRM)
R package **FastSparseGRM** provides functions and a pipeline to efficiently calculate genetic principal components (PCs) and the ancestry-adjusted sparse genetic relatedness matrix (GRM). It accounts for population heterogeneity using genetic PCs which are automatically calculated as part of the pipeline. The genetic PCs can be used as fixed effect covariates to account for the population stratification and the sparse GRM can be used to model the random effects to account for the sample relatedness in a mixed effects phenotype-genotype association testing model implemented in MetaSTAARlite. For more details on usage, please see the R package <a href="https://github.com/rounakdey/FastSparseGRM">**FastSparseGRM**</a> and <a href="https://doi.org/10.21203/rs.3.rs-5343361/v1">manuscript</a>.

## Generate study-specific variant summary statistics using MetaSTAARlite Worker
### Step 0: Preparation for association analysis of whole-genome/whole-exome sequencing studies
#### Script: <a href="Association_Analysis_PreStep.r">**Association_Analysis_PreStep.r**</a>
#### Input: aGDS files of all 22 chromosomes. For more details, please see the R script.
#### Output: `agds_dir.Rdata`, `Annotation_name_catalog.Rdata`.
* `agds_dir.Rdata`: a vector containing directory of GDS/aGDS files of all chromosomes. <br>
* `Annotation_name_catalog.Rdata`: a data frame containing the annotation name and the corresponding channel name in the aGDS file. Alternatively, one can skip this part in the R script by providing `Annotation_name_catalog.csv` with the same information. An example of `Annotation_name_catalog.csv` can be found <a href="https://github.com/xihaoli/STAARpipeline-Tutorial/blob/main/FAVORannotator_csv/Annotation_name_catalog.csv">here</a>. <br>

### Step 1: Fit STAAR null model
#### Script: <a href="https://github.com/xihaoli/STAARpipeline-Tutorial/blob/main/STAARpipeline_Null_Model.r">**STAARpipeline_Null_Model.r**</a> or <a href="https://github.com/xihaoli/STAARpipeline-Tutorial/blob/main/STAARpipeline_Null_Model_GENESIS.r">**STAARpipeline_Null_Model_GENESIS.r**</a>
* `STAARpipeline_Null_Model.r` fits the STAAR null model using the STAARpipeline package. <br>
* `STAARpipeline_Null_Model_GENESIS.r` fits the null model using the GENESIS package and convert it to the STAAR null model using the STAARpipeline package.
#### Input: Phenotype data and sparse genetic relatedness matrix. For more details, please see the R scripts.
#### Output: a Rdata file of the STAAR null model.

### Step 2.1: Generate variant summary statistics for individual (single-variant) analysis using MetaSTAARlite Worker
#### Script: <a href="MetaSTAARlite_worker_Individual_Analysis.r">**MetaSTAARlite_worker_Individual_Analysis.r**</a>
Generate and store variant summary statistics (*score statistics*) using MetaSTAARlite Worker.
#### Input: aGDS files and the STAAR null model. For more details, please see the R script.
#### Output: 294 Rdata files with the user-defined names.
Note: The number of output files is the summation of the column "individual_analysis_num" for the object in <a href="jobs_num.csv">`jobs_num.csv`</a> (**please use this file for all studies)**, which is 294.

### Step 2.2: Generate variant summary statistics for gene-centric coding analysis using MetaSTAARlite Worker
#### Script: <a href="MetaSTAARlite_worker_Gene_Centric_Coding.r">**MetaSTAARlite_worker_Gene_Centric_Coding.r**</a> and <a href="MetaSTAARlite_worker_Gene_Centric_Coding_Long_Masks.r">**MetaSTAARlite_worker_Gene_Centric_Coding_Long_Masks.r**</a>
Generate and store variant summary statistics (*score statistics*, *functional annotations*, and *sparse weighted LD matrices*) for coding rare variants using the MetaSTAARlite package. The gene-centric coding analysis provides five functional categories to aggregate coding rare variants of each protein-coding gene: (1) putative loss of function (stop gain, stop loss and splice) RVs, (2) missense RVs, (3) disruptive missense RVs, (4) putative loss of function and disruptive missense RVs, and (5) synonymous RVs. <br>
* `MetaSTAARlite_worker_Gene_Centric_Coding.r` generates and stores variant summary statistics for all protein-coding genes across the genome. There are 379 jobs using this script. <br>
* `MetaSTAARlite_worker_Gene_Centric_Coding_Long_Masks.r` generates and stores variant summary statistics for some specific long masks, and might require larger memory compared to `MetaSTAARlite_worker_Gene_Centric_Coding.r`. There are 2 jobs using this script.
#### Input: aGDS files and the STAAR null model. For more details, please see the R scripts.
#### Output: 381 * 2 = 762 Rdata files with the user-defined names for protein-coding genes.
Note: For rare variant meta-analysis (e.g. combined MAF < 1%), one can set `cov_maf_cutoff = 0.05` (by default) when generating sparse weighted LD matrices for each study.

### Step 2.3: Generate variant summary statistics for gene-centric noncoding analysis using MetaSTAARlite Worker
#### Script: <a href="MetaSTAARlite_worker_Gene_Centric_Noncoding.r">**MetaSTAARlite_worker_Gene_Centric_Noncoding.r**</a>, <a href="MetaSTAARlite_worker_Gene_Centric_Noncoding_Long_Masks.r">**MetaSTAARlite_worker_Gene_Centric_Noncoding_Long_Masks.r**</a>, <a href="MetaSTAARlite_worker_Gene_Centric_ncRNA.r">**MetaSTAARlite_worker_Gene_Centric_ncRNA.r**</a> and <a href="MetaSTAARlite_worker_Gene_Centric_ncRNA_Long_Masks.r">**MetaSTAARlite_worker_Gene_Centric_ncRNA_Long_Masks.r**</a>
Generate and store variant summary statistics (*score statistics*, *functional annotations*, and *sparse weighted LD matrices*) for noncoding rare variants using the MetaSTAARlite package. The gene-centric noncoding meta-analysis provides eight functional categories of regulatory regions to aggregate noncoding rare variants: (1) promoter RVs overlaid with CAGE sites, (2) promoter RVs overlaid with DHS sites, (3) enhancer RVs overlaid with CAGE sites, (4) enhancer RVs overlaid with DHS sites, (5) untranslated region (UTR) RVs, (6) upstream region RVs, (7) downstream region RVs, and (8) noncoding RNA (ncRNA) RVs. <br>
* `MetaSTAARlite_worker_Gene_Centric_Noncoding.r` generates and stores variant summary statistics for all protein-coding genes across the genome. There are 379 jobs using this script. <br>
* `MetaSTAARlite_worker_Gene_Centric_Noncoding_Long_Masks.r` generates and stores variant summary statistics for some specific long masks, and might require larger memory compared to `MetaSTAARlite_worker_Gene_Centric_Noncoding.r`. There are 8 jobs using this script. <br>
* `MetaSTAARlite_worker_Gene_Centric_ncRNA.r` generates and stores variant summary statistics for ncRNA genes across the genome. There are 222 jobs using this script. <br> 
* `MetaSTAARlite_worker_Gene_Centric_ncRNA_Long_Masks.r` generates and stores variant summary statistics  for some specific long masks, and might require larger memory compared to `MetaSTAARlite_worker_Gene_Centric_ncRNA`. There is 1 job using this script. 
#### Input: aGDS files and the STAAR null model. For more details, please see the R scripts.
#### Output: 387 * 2 = 774 Rdata files with the user-defined names for protein-coding genes and 223 * 2 = 446 Rdata files with the user-defined names for ncRNA genes.
Note: For rare variant meta-analysis (e.g. combined MAF < 1%), one can set `cov_maf_cutoff = 0.05` (by default) when generating sparse weighted LD matrices for each study.

### Step 2.C: Generate variant summary statistics for custom mask analysis using MetaSTAARlite Worker
#### Script: <a href="MetaSTAARlite_worker_Custom_Mask.r">**MetaSTAARlite_worker_Custom_Mask.r**</a>
Generate and store variant summary statistics (*score statistics*, *functional annotations*, and *sparse weighted LD matrices*) for user-defined custom masks of rare variants using the MetaSTAARlite package. <br>
#### Input: aGDS files, the STAAR null model, and the custom mask definition file. For more details, please see the R script.
#### Output: Rdata files with the user-defined names for custom masks.
An example file for a list of custom masks (5-column "CHR-POS-REF-ALT-MaskName" format) in the *LDLR* locus is given in <a href="custom_mask_LDLR.csv">`custom_mask_LDLR.csv`</a>.

Note: For rare variant meta-analysis (e.g. combined MAF < 1%), one can set `cov_maf_cutoff = 0.05` (by default) when generating sparse weighted LD matrices for each study. Users may want to perform parallel computation by splitting the analysis of custom masks into multiple jobs.

## Association meta-analysis using MetaSTAARlite
### Step 3: Individual (single-variant) meta-analysis
#### Script: <a href="MetaSTAARlite_Individual_Analysis.r">**MetaSTAARlite_Individual_Analysis.r**</a>
Perform single-variant meta-analysis for common and low-frequency variants across the genome using the MetaSTAARlite package. 
#### Input: Variant summary statistics files from Step 2.1 for each participating study. For more details, please see the R script.
#### Output: 294 Rdata files with the user-defined names.
Note: The number of output files is the summation of the column "individual_analysis_num" for the object in <a href="jobs_num.csv">`jobs_num.csv`</a>, which is 294.

### Step 4.1: Gene-centric coding meta-analysis
#### Script: <a href="MetaSTAARlite_Gene_Centric_Coding.r">**MetaSTAARlite_Gene_Centric_Coding.r**</a> and <a href="MetaSTAARlite_Gene_Centric_Coding_Long_Masks.r">**MetaSTAARlite_Gene_Centric_Coding_Long_Masks.r**</a>
Perform gene-centric meta-analysis for coding rare variants using the MetaSTAARlite package. The gene-centric coding meta-analysis provides five functional categories to aggregate coding rare variants of each protein-coding gene: (1) putative loss of function (stop gain, stop loss and splice) RVs, (2) missense RVs, (3) disruptive missense RVs, (4) putative loss of function and disruptive missense RVs, and (5) synonymous RVs. <br>
* `MetaSTAARlite_Gene_Centric_Coding.r` performs gene-centric coding meta-analysis for all protein-coding genes across the genome. There are 379 jobs using this script. <br>
* `MetaSTAARlite_Gene_Centric_Coding_Long_Masks.r` performs gene-centric coding meta-analysis for some specific long masks, and might require larger memory compared to `MetaSTAARlite_Gene_Centric_Coding.R`. There are 2 jobs using this script.
#### Input: Variant summary statistics files from Step 2.2 for each participating study. For more details, please see the R scripts.
#### Output: 381 Rdata files with the user-defined names.

### Step 4.2: Gene-centric noncoding meta-analysis
#### Script: <a href="MetaSTAARlite_Gene_Centric_Noncoding.r">**MetaSTAARlite_Gene_Centric_Noncoding.r**</a>, <a href="MetaSTAARlite_Gene_Centric_Noncoding_Long_Masks.r">**MetaSTAARlite_Gene_Centric_Noncoding_Long_Masks.r**</a>, <a href="MetaSTAARlite_Gene_Centric_ncRNA.r">**MetaSTAARlite_Gene_Centric_ncRNA.r**</a> and <a href="MetaSTAARlite_Gene_Centric_ncRNA_Long_Masks.r">**MetaSTAARlite_Gene_Centric_ncRNA_Long_Masks.r**</a>
Perform gene-centric meta-analysis for noncoding rare variants using the MetaSTAARlite package. The gene-centric noncoding meta-analysis provides eight functional categories of regulatory regions to aggregate noncoding rare variants: (1) promoter RVs overlaid with CAGE sites, (2) promoter RVs overlaid with DHS sites, (3) enhancer RVs overlaid with CAGE sites, (4) enhancer RVs overlaid with DHS sites, (5) untranslated region (UTR) RVs, (6) upstream region RVs, (7) downstream region RVs, and (8) noncoding RNA (ncRNA) RVs. <br>
* `MetaSTAARlite_Gene_Centric_Coding.r` performs gene-centric noncoding meta-analysis for all protein-coding genes across the genome. There are 379 jobs using this script. <br>
* `MetaSTAARlite_Gene_Centric_Coding_Long_Masks.r` performs gene-centric noncoding meta-analysis for some specific long masks, and might require larger memory compared to `MetaSTAARlite_Gene_Centric_Coding.r`. There are 8 jobs using this script. <br>
* `MetaSTAARlite_Gene_Centric_ncRNA.r` performs gene-centric noncoding meta-analysis for ncRNA genes across the genome. There are 222 jobs using this script. <br> 
* `MetaSTAARlite_Gene_Centric_ncRNA_Long_Masks.r` performs gene-centric noncoding meta-analysis for some specific long masks, and might require larger memory compared to `MetaSTAARlite_Gene_Centric_ncRNA`. There is 1 job using this script. 
#### Input: Variant summary statistics files from Step 2.3 for each participating study. For more details, please see the R scripts.
#### Output: 387 Rdata files with the user-defined names for protein-coding genes and 223 Rdata files with the user-defined names for ncRNA genes.

### Step C: Custom mask meta-analysis
#### Script: <a href="MetaSTAARlite_Custom_Mask.r">**MetaSTAARlite_Custom_Mask.r**</a>
Perform meta-analysis for user-defined custom masks of rare variants using the MetaSTAARlite package. <br>
#### Input: Variant summary statistics files from Step 2.C for each participating study. For more details, please see the R script.
#### Output: Rdata files with the user-defined names for custom masks.
Note: Users may want to perform parallel computation by splitting the analysis of custom masks into multiple jobs.

## Summarization and visualization of association analysis results using MetaSTAARliteSummary
### Step 5.0 (Optional): Select independent variants from a known variants list to be used in conditional meta-analysis
An example file for a list of known total cholesterol associated variants (4-column "CHR-POS-REF-ALT" format) in the *LDLR* locus is given in <a href="known_loci_LDLR.csv">`known_loci_LDLR.csv`</a>.

Note: It is typically assumed that variants in `known_loci` for conditional analysis are present in each participating study of the meta-analysis, which implies that the MAFs for these variants would not be extremely rare.

### Step 5.1.1: Summarize individual (single-variant) meta-analysis results
#### Script: <a href="MetaSTAARliteSummary_Individual_Analysis.r">**MetaSTAARliteSummary_Individual_Analysis.r**</a>
Summarize single-variant meta-analysis results.
#### Input: Individual analysis results generated by MetaSTAARlite. For more details, please see the R script.
#### Output: The summary includes the Manhattan plot and Q-Q plot.

### Step 5.1.2: Generate variant summary statistics for individual (single-variant) conditional analysis using MetaSTAARlite Worker
#### Script: <a href="MetaSTAARlite_worker_Individual_Analysis_cond.r">**MetaSTAARlite_worker_Individual_Analysis_cond.r**</a>
Generate and store variant summary statistics (*score statistics* and *conditional LD matrices*) using MetaSTAARlite Worker.
#### Input: aGDS files, the STAAR null model and a list of known variants. For more details, please see the R script.
#### Output: Rdata files with the user-defined names.

### Step 5.1.3: Conditional individual (single-variant) meta-analysis
#### Script: <a href="MetaSTAARlite_Individual_Analysis_cond.r">**MetaSTAARlite_Individual_Analysis_cond.r**</a>
Perform conditional meta-analysis of unconditionally significant variants by adjusting a list of known variants.
#### Input: Variant summary statistics files from Step 5.1.2 for each participating study and (significant) individual analysis results from Step 5.1.1. For more details, please see the R script.
#### Output: Conditional p-values of unconditionally significant variants.

### Step 5.2.1: Summarize gene-centric coding meta-analysis results
#### Script: <a href="MetaSTAARliteSummary_Gene_Centric_Coding.r">**MetaSTAARliteSummary_Gene_Centric_Coding.r**</a>
Summarize gene-centric coding meta-analysis results.
#### Input: Gene-centric coding analysis results generated by MetaSTAARlite. For more details, please see the R script.
#### Output: The summary includes the Manhattan plot, Q-Q plot.

### Step 5.2.2: Generate variant summary statistics for gene-centric coding conditional analysis using MetaSTAARlite Worker
#### Script: <a href="MetaSTAARlite_worker_Gene_Centric_Coding_cond.r">**MetaSTAARlite_worker_Gene_Centric_Coding_cond.r**</a>
Generate and store variant summary statistics (*score statistics*, *functional annotations*, *sparse weighted LD matrices*, and *conditional LD matrices*) for coding rare variants using MetaSTAARlite Worker.
#### Input: aGDS files, the STAAR null model, and a list of known variants. For more details, please see the R script.
#### Output: Rdata files with the user-defined names.

### Step 5.2.3: Conditional gene-centric coding meta-analysis
#### Script: <a href="MetaSTAARlite_Gene_Centric_Coding_cond.r">**MetaSTAARlite_Gene_Centric_Coding_cond.r**</a>
Perform conditional meta-analysis of unconditionally significant coding masks by adjusting a list of known variants.
#### Input: Variant summary statistics files from Step 5.2.2 for each participating study. For more details, please see the R script.
#### Output: Conditional p-values of unconditionally significant coding masks.

### Step 5.3.1: Summarize gene-centric noncoding meta-analysis results
#### Script: <a href="MetaSTAARliteSummary_Gene_Centric_Noncoding.r">**MetaSTAARliteSummary_Gene_Centric_Noncoding.r**</a>
Summarize gene-centric noncoding meta-analysis results.
#### Input: Gene-centric noncoding analysis results generated by MetaSTAARlite. For more details, please see the R script.
#### Output: The summary includes the Manhattan plots and Q-Q plots.

### Step 5.3.2: Generate variant summary statistics for gene-centric noncoding conditional analysis using MetaSTAARlite Worker
#### Script: <a href="MetaSTAARlite_worker_Gene_Centric_Noncoding_cond.r">**MetaSTAARlite_worker_Gene_Centric_Noncoding_cond.r**</a> and <a href="MetaSTAARlite_worker_Gene_Centric_ncRNA_cond.r">**MetaSTAARlite_worker_Gene_Centric_ncRNA_cond.r**</a>
Generate and store variant summary statistics (*score statistics*, *functional annotations*, *sparse weighted LD matrices*, and *conditional LD matrices*) for noncoding rare variants using MetaSTAARlite Worker.
#### Input: aGDS files, the STAAR null model, and a list of known variants. For more details, please see the R scripts.
#### Output: Rdata files with the user-defined names.

### Step 5.3.3: Conditional gene-centric noncoding meta-analysis
#### Script: <a href="MetaSTAARlite_Gene_Centric_Noncoding_cond.r">**MetaSTAARlite_Gene_Centric_Noncoding_cond.r**</a> and <a href="MetaSTAARlite_Gene_Centric_ncRNA_cond.r">**MetaSTAARlite_Gene_Centric_ncRNA_cond.r**</a>
Perform conditional meta-analysis of unconditionally significant noncoding masks by adjusting a list of known variants.
#### Input: Variant summary statistics files from Step 5.3.2 for each participating study. For more details, please see the R scripts.
#### Output: Conditional p-values of unconditionally significant noncoding masks.

### Step 5.C.1: Summarize custom mask meta-analysis results
Note: Users may develop custom scripts to summarize the analysis results by following Step 5.2.1/5.3.1, as the number of output files from Step C may vary.

### Step 5.C.2: Generate variant summary statistics for custom mask conditional analysis using MetaSTAARlite Worker
#### Script: <a href="MetaSTAARlite_worker_Custom_Mask_cond.r">**MetaSTAARlite_worker_Custom_Mask_cond.r**</a>
Generate and store variant summary statistics (*score statistics*, *functional annotations*, *sparse weighted LD matrices*, and *conditional LD matrices*) for user-defined custom masks of rare variants using MetaSTAARlite Worker.
#### Input: aGDS files, the STAAR null model, the custom mask definition file, and a list of known variants. For more details, please see the R script.
#### Output: Rdata files with the user-defined names.

### Step 5.C.3: Conditional custom mask meta-analysis
#### Script: <a href="MetaSTAARlite_Custom_Mask_cond.r">**MetaSTAARlite_Custom_Mask_cond.r**</a>
Perform conditional meta-analysis of unconditionally significant user-defined custom masks by adjusting a list of known variants.
#### Input: Variant summary statistics files from Step 5.C.2 for each participating study. For more details, please see the R script.
#### Output: Conditional p-values of unconditionally significant user-defined custom masks.

### An example of batch job submission scripts for these analyses can be found <a href="/batch jobs">**here**</a>.

