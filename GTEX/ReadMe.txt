# GTEx Sex Chromosome Complement Analysis

## Overview
This project analyzes gene expression data from the Genotype-Tissue Expression (GTEx) project, focusing on gene expression related to the X and Y chromosomes. The analysis aims to explore sex chromosome complements and their relationship with gene expression.

## Contents
1. **RMarkdown Files**
   - `GTEx_Sex_Chromosome_Analysis.Rmd`: Generates visualizations that illustrate the relationship between gene expression levels and inferred sex chromosome complements.

2. **Data**
   - `gtex_expression_data.csv`: Gene expression dataset containing expression levels for various genes, including XIST and Y chromosome genes across different tissues and samples.
   - `inferred_sex_chromosome_complement.csv`: Output from the inferred sex chromosome complements determination based on the analysis of gene expression data.

3. **Output**
   - Generated plots in PDF format, saved in the output directory.

## Prerequisites
Ensure you have the following R packages installed:
- `ggplot2`
- `dplyr`
- `readr`
- `tidyverse`

## How to Run the Analysis
1. **Set the Working Directory**
   Update the paths in the RMarkdown file to point to the correct directories on your system.

2. **Load the Data**
   The data is loaded from CSV files. Ensure that the specified paths in the `data_directory` variable are correct.

3. **Execute the RMarkdown File**
   You can run the analysis in RStudio by opening the RMarkdown file and clicking on the "Knit" button to generate an HTML report.

## Analysis Steps
### Step 1: Data Loading and Cleaning
- Load the gene expression data and inferred sex chromosome complements.
- Merge and clean the datasets, retaining only the relevant columns for analysis.

### Step 2: Generate Visualizations
- Create visualizations to illustrate the relationship between gene expression levels and inferred sex chromosome complements.
- Save the generated plots in the output directory.

## Output
- Visualizations showing the relationship between gene expression and sex chromosome complements.

## License


## Acknowledgments
Thanks to the GTEx project for providing the data used in this analysis and to the researchers who contributed to the understanding of sex chromosomes and gene expression.
