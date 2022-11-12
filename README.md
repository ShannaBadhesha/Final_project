# Final Project Outline

## Title
Differential Gene Expression in TCGA within stage 1-4 Ovarian Serous Cystadenocarcinoma comparing radiation treatment vs. non-radiation treatment using DeSEQ2

## Author
Shanna Badhesha

## Overview of the Project
I will identify differentially expressed genes between Ovarian Serous Cystadenocarcinoma radiation treatment vs. non-radiation treatment. This analysis will utilize the package DeSEQ2 and the http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html vignette. For this analysis, I will be using the TCGA-OV cohort. I have identified a total of 585 STAR-counts files with 450 radiation treatment and 135 non-radiation treatment files. 

## Data
I will use the data from the following website: https://portal.gdc.cancer.gov/repository. There are 585 STAR-count files with 450 radiation treatment and 135 non-radiation treatment files. The specific files are available [here](https://portal.gdc.cancer.gov/repository?facetTab=files&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22ovary%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.analysis.workflow_type%22%2C%22value%22%3A%5B%22STAR%20-%20Counts%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%5D%7D).

## Milestone 1
I will fully load my data into a dataframe using the the HT-SEQ steps in the vignette listed above. 

## Milestone 2
I will do an initial run of the entire vignette listed above for analysis. This will include having data loaded in (completed during milestone 1) and run through the entire vignette. I will be seeking feedback during this step. I will begin the process of uploading documentation and defining my results.

## Deliverable

Upon completion, I will be uploading documentation and descriptions of my analysis and results. 
