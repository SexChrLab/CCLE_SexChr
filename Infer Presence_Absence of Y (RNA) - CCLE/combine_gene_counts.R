
# Load required packages
library(dplyr)

# Define the directory where your CSV files are located
data_dir <- "D:/CCLE_Project/Infer_Y_output/"

# Output directory
output_dir <- "D:/CCLE_Project/data/"

# List of genes in order of their location
genes <- c("XIST","SRY", "RPS4Y1", "ZFY", "TGIF2LY", "PCDH11Y", "TBL1Y", "USP9Y",
           "DDX3Y", "UTY", "TMSB4Y", "NLGN4Y", "KDM5D", "EIF1AY", "RPS4Y2")

# Initialize an empty list to store data frames
gene_data_list <- list()

# Loop through each gene and read its corresponding CSV file
for (gene in genes) {
  gene_file <- paste0(data_dir, "predicted_sex_", gene, ".csv")
  gene_data <- read.csv(gene_file, header = TRUE)
  
  # Rename columns to include the gene name
  colnames(gene_data) <- c("cell_line", 
                           "reported_sex", 
                           paste0(gene, "_expression"), 
                           paste0(gene, "_log_expression"), 
                           paste0(gene, "_expression_level"))
  
  # Select and reorder columns
  gene_data <- gene_data %>%
    select(cell_line, reported_sex, ends_with("expression_level"), ends_with("expression"), ends_with("log_expression"))
  
  # Add the gene data frame to the list
  gene_data_list[[gene]] <- gene_data
}

# Merge all gene data frames into one by "cell_line"
combined_gene_data <- Reduce(function(x, y) left_join(x, y, by = c("cell_line", "reported_sex")), gene_data_list)

# Export the combined data frame to a CSV file
write.csv(combined_gene_data, paste0(output_dir, "CCLE_all_combined_gene_expression.csv"), row.names = FALSE)
