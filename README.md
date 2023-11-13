# PRS_tutorial
## INTRODUCTION

The purpose of this tutorial is to furnish a step-by-step guide on the execution of a Polygenic Risk Score (PRS) analysis, accompanied by an introduction to fundamental genetic principles germane to this methodology.

The tutorial is divided into two sections, First, a series of basic concepts about genetics and PRS are introduced at a theoretical level. Secondly, the tutorial itself consists of 6 steps which are detailed below.

**Step-by-step PRS analysis:**

0. [STEP 0: Data selection + QC check](#step0)
1. [STEP 1: SNPs Extraction](#step1)
2. [STEP 2: SNPs Alignament check](#step2)
3. [STEP 3: Quality control (QC)](#step3)
4. [STEP 4: PRS calculation](#step4)
5. [STEP 5: PRS resultsn](#step5)



## Software requirements


- **Linux terminal** access with the following tools installed:
    - **PLINK 2** - https://www.cog-genomics.org/plink/2.0/ 
    - **bcftools** - http://www.htslib.org/download/ 
<br>    
- **R** 

### - R Library

R packages required:

- library(readxl)
- library(data.table)
- library(dplyr)
- library(ggplot2)


## Required Data

To carry out this PRS analysis it is necessary two types of data. Here is the link to download the data used in the tutorial. If you want to use your own data, see what is necessary in the section [STEP 0 - Data selection + QC check](#step0).

  - **SNPs DATA**. 
    List of SNPs associated with a characteristic of interest. This tutorial will use a list of 205 SNPs associated with colon cancer (CRC) from the article: [Deciphering colorectal cancer genetics through multi-omic analysis of 100,204 cases and 154,587 controls of European and east Asian ancestries](https://www.nature.com/articles/s41588-022-01222-9)
    
   You can download it directly here: [**Download - 1_CRC_SNPs_list_205_Nature_2023**](./resources)
  
    
<br>

  - **TARGET DATA**. 
  Imputed genotyping data from a population of interest. Simulated data based on GenRisk data will be used here. For reasons of space, only the data of the patients with the SNPs associated with CRC can be downloaded (document obtained in STEP 1) [**Download - GenRisk_all_CRC_SNPs_perm.vcf **](./resources).
  
  To be used from the section [Step 1 - SNP extraction. 2 VCF files concatenation](#vcf_conca)
  
  - **METADATA**.
  Phenotypic data of the target data. In this case which of the samples is controlled and which CRC. [Metadata_GenRIsk - ](./resources)


<div class="alert alert-danger">
**Attention**: Whether you wish to utilize your own data or existing datasets, it's important to note that the data used should have already undergone a minimum quality control (QC) process, with hg37 as the reference genome. For the target data, it has been imputed using the 1000 Genomes reference.

For more information regarding the format of the data to be used or how to adapt see [here](#step0)

</div>
