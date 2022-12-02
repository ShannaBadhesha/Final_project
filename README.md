# Final Project Outline

## Title
Differential Gene Expression in TCGA within stage 1-4 Ovarian Serous Cystadenocarcinoma comparing radiation treatment using carboplatin vs. non-radiation treatment without carboplatin using DeSEQ2

## Author
Shanna Badhesha

## Overview of the Project
I will identify differentially expressed genes between Ovarian Serous Cystadenocarcinoma radiation treatment using carboplatin vs. non-radiation treatment without carboplatin. This analysis will utilize the package DeSEQ2 and the http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html vignette. For this analysis, I will be using the TCGA-OV cohort. I have identified a total of 40 STAR-counts files with 20 radiation treatment and 20 non-radiation treatment files. 

## Data
I will use the data from the following website: https://portal.gdc.cancer.gov/repository. There are 585 STAR-count files with 450 radiation treatment and 135 non-radiation treatment files in total. I will be using 40 STAR-counts files with 20 samples per group. The specific files are available [here](https://portal.gdc.cancer.gov/repository?facetTab=files&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22ovary%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.analysis.workflow_type%22%2C%22value%22%3A%5B%22STAR%20-%20Counts%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%5D%7D).

## Milestone 1
I have downloaded the data and matched the file names to the corresponding TCGA IDs. I have isolated the necessary unstranded column from each file and merged them together for the analysis. I have created another spreadsheet identifying the TCGA IDs for each of the two treatment groups.

## Known Issues
The original dataset included over 500 samples. I narrowed it down to 300 and then further narrowed it down to 40 samples for the analysis. Since there were no filters avaliable, I chose 20 samples for each group at random to speed up the process of the analysis. 

## Milestone 2
I will do an initial run of the entire vignette listed above for analysis. This will include having data loaded in (completed during milestone 1) and run through the entire vignette. I will be seeking feedback during this step. I will begin the process of uploading documentation and defining my results.

## Deliverable
Below is the completed vignette with my 40 samples. 

## Organzing Sample Files

The STAR-Count files provide a lot of information that I will not be needing for the analysis. For the purpose of this vignette, we will extract the unstranded gene count column from each file and merge them together.

Create a ```gene.id``` file that includes all the gene id's from one of the sample files: 
	```awk '{print $1}’ TCGA-04-1331.tsv > TCGA-04-1331.txt```

Extract the unstranded column (column 4) from each tsv file:
	```awk '{print $4}’ TCGA-04-1331.tsv > TCGA-04-1331.txt```

Remove the first 6 lines from each tsv file since they will not be needed for our analysis:
	```sed -i .bak '1,6d' 'TCGA-04-1332.txt'```

Add a header to each of the file using ```vi```. The header should be the sample ID. 

Merge the files together using the ```paste``` command: 
	```paste gene.id *.txt > merged_files.txt```

### Sample Information Table 

Create a txt file using Microsoft Excel with the Sample IDs and their condition. The Sample IDs should be in the same order as the merged_file.txt otherwise DESeq2 will not accept the file. 

## DESeq2 Vignette




