---
title: "R Notebook"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 

```{r}
# 
# title: "Using Combined TCGA data for gene loop"
# author: "Ilsa Rodriguez"
# date: "`r Sys.Date()`"
TCGA_RNA_data_all <- read_csv("TCGA_RNA_data_all.csv")


CancerType <- TCGA_RNA_data_all[TCGA_RNA_data_all$project_id == "LIHC", ]

#Calling XY and XX genes
# Define female genes
female_genes <- c("XIST")

# Define male genes
male_genes <- c("USP9Y","UTY" ,"ZFY" , "SRY" ,"AMELY" , "DDX3Y",  "EIF1AY" ,  "KDM5D", "NLGN4Y" ,  "PRKY"  ,"TMSB4Y")


# Initialize a dataframe to store results for all genes
all_genes_data <- data.frame(RowID = CancerType$rowname, 
                             annotated_sex = CancerType$annotated_sex, 
                             XIST_status = CancerType$status_XIST, 
                             Y_status = CancerType$status_Y)

for (gene in c(female_genes, male_genes)) {
  if (gene %in% female_genes) {
    # Assuming you want to use log1p for female genes
    gene_expression <- log1p(CancerType[[gene]]) 
  } else {
    # Assuming you want to use  log1p for male genes
    gene_expression <- log1p(CancerType[[gene]]) 
  }
  gene_expression[is.infinite(gene_expression)] <- NA
  
  # Set thresholds based on the sex
  
  # low thresholds in XY means that it is an XX gene
  if (gene %in% female_genes) {
    high_threshold <- mean(CancerType[[gene]][CancerType$annotated_sex == 'male'], na.rm = TRUE) + (4.5 * sd(CancerType[[gene]][CancerType$annotated_sex == 'male'], na.rm = TRUE))
    low_threshold <- mean(CancerType[[gene]][CancerType$annotated_sex == 'male'], na.rm = TRUE) + (3 * sd(CancerType[[gene]][CancerType$annotated_sex == 'male'], na.rm = TRUE))
    
    
  } else {
    high_threshold <- mean(CancerType[[gene]][CancerType$annotated_sex == 'female'], na.rm = TRUE) + (4.5 * sd(CancerType[[gene]][CancerType$annotated_sex == 'female'], na.rm = TRUE))
    low_threshold <- mean(CancerType[[gene]][CancerType$annotated_sex == 'female'], na.rm = TRUE) + (3* sd(CancerType[[gene]][CancerType$annotated_sex == 'female'], na.rm = TRUE))
    
  }
  
  
  
  # Predict sex based on gene expression thresholds
  if (gene %in% male_genes) {
    # For Y chromosome associated genes, higher expression suggests "male"
    predicted_sex <- ifelse(gene_expression >= high_threshold, "male",
                            ifelse(gene_expression <= low_threshold, "female", "cannot_predict"))
  } else {
    # For X chromosome gene expression, higher expression suggests "female"
    predicted_sex <- ifelse(gene_expression >= high_threshold, "female",
                            ifelse(gene_expression <= low_threshold, "male", "cannot_predict"))
  }
  # Prepare data for plotting
  gene_data_for_plot <- data.frame(annotated_sex = CancerType$annotated_sex, expression = gene_expression)
  
  p <- ggplot(gene_data_for_plot, aes(x = annotated_sex, y = expression, fill = annotated_sex)) +
    geom_violin(trim = FALSE) +
    scale_y_continuous(trans = "log10") +
    scale_fill_manual(values = c("grey", "orange")) +
    geom_jitter(size = 0.75) +
    ylab(paste0(gene, " Expression (", if(gene %in% female_genes) "log10" else "log", "-transformed)")) +
    geom_hline(yintercept = high_threshold, linetype="dashed", color = "maroon") +
    geom_hline(yintercept = low_threshold, linetype="dashed", color = "blue") +
    theme_light() 
  # Change scale if you want log10
  print(p)
  
  
  
    comparison_reported_predicted <- table(CancerType$annotated_sex, predicted_sex)

    # Visualize the comparison
   x <- barplot(comparison_reported_predicted, 
            beside = TRUE, 
            legend.text = TRUE,
            main = paste(gene, "gene for predicting sex in TCGA-LIHC"),
            xlab = "Reported Sex", 
            ylab = "Number of Samples",
            col = c("grey", "orange"),
            args.legend = list(x = "topleft", bty = "n", inset = c(0.05, 0))) 
    y <- as.matrix(comparison_reported_predicted)

text(x, y + 20, labels = as.character(y))
  
  
  # Add expression data to the dataframe
  all_genes_data[[paste0(gene, "_expression")]] <- CancerType[[gene]]
  all_genes_data[[paste0(gene, "_predicted_sex")]] <- predicted_sex
  all_genes_data[[paste0(gene, "_predicted_expression")]] <- ifelse(CancerType[[gene]] >= high_threshold, "high_expression", ifelse(CancerType[[gene]] <= low_threshold, "low_expression", "intermediate_expression"))
    # Comparison of reported_sex to predicted_sex
  
  
violin6_fname <- paste("Log_Violin", gene, ".png", sep = "_")
ggsave(violin6_fname, p, width = 10, height = 8, units = "in", dpi = 300)

}

# Save the comprehensive results to a CSV file
write.csv(all_genes_data, "expression_data_all_genes_with_predicted_sex.csv", row.names = FALSE)
```
```{r}
# Load your data here
counts_plus <- read.table(file = "C:/Users/livin/ASU FILES/Reseach Files/Genomic Sex Chromosome/TCGA_RNAseq_counts/outcomes TCGA/TCGA_RNA_data_all_modified.csv", sep = '\t', header = TRUE)
# Initialize a dataframe to store results for all genes

library(ggplot2)

#Calling XY and XX genes
# Define female genes
female_genes <- c("XIST")

# Define male genes
male_genes <- c("USP9Y","UTY"  ,   "ZFY" ,   "SRY" ,  "AMELY" ,   "DDX3Y",  "EIF1AY" ,  "KDM5D", "NLGN4Y" ,  "PRKY"  ,"TMSB4Y" )


# Initialize a dataframe to store results for all genes
all_genes_data <- data.frame(RowID  = rownames(counts_plus), 
                             annotated_sex = counts_plus$annotated_sex, 
                             XIST_status = counts_plus$status_XIST, 
                             Y_status = counts_plus$status_Y)

for (gene in c(female_genes, male_genes)) {
  if (gene %in% female_genes) {
    # Assuming you want to use log1p for female genes
    gene_expression <- log1p(counts_plus[[gene]]) 
  } else {
    # Assuming you want to use  log1p for male genes
    gene_expression <- log1p(counts_plus[[gene]]) 
  }
  gene_expression[is.infinite(gene_expression)] <- NA

  # Set thresholds based on the sex
  
  # low thresholds in XY means that it is an XX gene

  if (gene %in% female_genes) {
    high_threshold <- mean(counts_plus[[gene]][counts_plus$annotated_sex == 'male'], na.rm = TRUE) + (1 * sd(counts_plus[[gene]][counts_plus$annotated_sex == 'male'], na.rm = TRUE))
    low_threshold <- mean(counts_plus[[gene]][counts_plus$annotated_sex == 'male'], na.rm = TRUE) + (3 * sd(counts_plus[[gene]][counts_plus$annotated_sex == 'male'], na.rm = TRUE))

  
  } else {
    high_threshold <- mean(counts_plus[[gene]][counts_plus$annotated_sex == 'female'], na.rm = TRUE) + (4.5 * sd(counts_plus[[gene]][counts_plus$annotated_sex == 'female'], na.rm = TRUE))
  low_threshold <- mean(counts_plus[[gene]][counts_plus$annotated_sex == 'female'], na.rm = TRUE) + (3* sd(counts_plus[[gene]][counts_plus$annotated_sex == 'female'], na.rm = TRUE))
  
  }
 
 

  # Predict sex based on gene expression thresholds
 if (gene %in% male_genes) {
    # For Y chromosome associated genes, higher expression suggests "male"
    predicted_sex <- ifelse(gene_expression >= high_threshold, "male",
                            ifelse(gene_expression <= low_threshold, "female", "cannot_predict"))
  } else {
    # For X chromosome gene expression, higher expression suggests "female"
    predicted_sex <- ifelse(gene_expression >= high_threshold, "female",
                            ifelse(gene_expression <= low_threshold, "male", "cannot_predict"))
  }
  # Prepare data for plotting
  gene_data_for_plot <- data.frame(annotated_sex = counts_plus$annotated_sex, expression = gene_expression)

  p <- ggplot(gene_data_for_plot, aes(x = annotated_sex, y = expression, fill = annotated_sex)) +
    geom_violin(trim = FALSE) +
    scale_y_continuous(trans = "log10") +
    scale_fill_manual(values = c("grey", "orange")) +
    geom_jitter(size = 0.75) +
    ylab(paste0(gene, " Expression (", if(gene %in% female_genes) "log10" else "log", "-transformed)")) +
    geom_hline(yintercept = high_threshold, linetype="dashed", color = "maroon") +
   geom_hline(yintercept = low_threshold, linetype="dashed", color = "blue") +
    theme_light() 
     # Change scale if you want log10
  print(p)


  
  # Add expression data to the dataframe
  all_genes_data[[paste0(gene, "_log_expression")]] <- counts_plus[[gene]]
  all_genes_data[[paste0(gene, "_predicted_sex")]] <- predicted_sex
  all_genes_data[[paste0(gene, "_predicted_expression")]] <- ifelse(counts_plus[[gene]] >= high_threshold, "high_expression", ifelse(counts_plus[[gene]] <= low_threshold, "low_expression", "intermediate_expression"))
  
comparison_reported_predicted <- table(counts_plus$annotated_sex, predicted_sex)

    # Visualize the comparison
   x <- barplot(comparison_reported_predicted, 
            beside = TRUE, 
            legend.text = TRUE,
            main = paste(gene, "gene for predicting sex in TCGA-LIHC"),
            xlab = "Reported Sex", 
            ylab = "Number of Samples",
            col = c("grey", "orange"),
            args.legend = list(x = "topleft", bty = "n", inset = c(0.05, 0))) 
    y <- as.matrix(comparison_reported_predicted)

text(x, y + 20, labels = as.character(y))
  
  
  # Add expression data to the dataframe
  all_genes_data[[paste0(gene, "_expression")]] <- counts_plus[[gene]]
  all_genes_data[[paste0(gene, "_predicted_sex")]] <- predicted_sex
  all_genes_data[[paste0(gene, "_predicted_expression")]] <- ifelse(counts_plus[[gene]] >= high_threshold, "high_expression", ifelse(counts_plus[[gene]] <= low_threshold, "low_expression", "intermediate_expression"))
    # Comparison of reported_sex to predicted_sex
  
  
violin6_fname <- paste("Log_Violin", gene, ".png", sep = "_")
ggsave(violin6_fname, p, width = 10, height = 8, units = "in", dpi = 300)
} 


# Save the comprehensive results to a CSV file
write.csv(all_genes_data, "expression_data_all_genes_with_predicted_sex.csv", row.names = FALSE)


```
```{r}
##CCLE code being applied to TCGA, for predicted sex using calucated thresholds

#genes_of_interest <- c("XIST", "AR")
# Define the list of genes you're interested in
#genes_of_interest <- c( "XIST", "AR","AMELY", "DDX3Y", "EIF1AY", "KDM5D", "NLGN4Y", "PRKY", "TMSB4Y", "USP9Y", "UTY", "ZFY", "SRY")
genes_of_interest <- c("XIST", "DDX3Y" ,"USP9Y", "UTY" ,"ZFY")

# Initialize a dataframe to store results for all genes

library(ggplot2)

# Load your data here


female_genes <- c("XIST")
male_genes <- c( "DDX3Y", "EIF1AY", "KDM5D", "NLGN4Y", "PRKY", "TMSB4Y", "USP9Y", "UTY", "ZFY")

all_genes_data <- data.frame(RowID = rownames(counts_plus), 
                             annotated_sex = counts_plus$annotated_sex, 
                             XIST_status = counts_plus$status_XIST, 
                             Y_status = counts_plus$status_Y)

for (gene in c(female_genes, male_genes)) {
  # Apply appropriate log transformation based on the type of gene
  if (gene %in% female_genes) {
    gene_expression <- log10(counts_plus[[gene]]) # for natural log use log1p
  } else {
    gene_expression <- log10(counts_plus[[gene]]) # for natural log use log1p, adjust if you want log10
  }
  gene_expression[is.infinite(gene_expression)] <- NA 
  # Set thresholds based on the sex
  if (gene %in% female_genes) {
    # Calculate thresholds for female-associated genes
    high_threshold <- mean(gene_expression[counts_plus$annotated_sex == 'male'], na.rm = TRUE) +
                      (4.5 * sd(gene_expression[counts_plus$annotated_sex == 'male'], na.rm = TRUE))
    low_threshold <- mean(gene_expression[counts_plus$annotated_sex == 'male'], na.rm = TRUE) +
                     (3 * sd(gene_expression[counts_plus$annotated_sex == 'male'], na.rm = TRUE))
  } else {
    # Calculate thresholds for male-associated genes
    high_threshold <- mean(gene_expression[counts_plus$annotated_sex == 'female'], na.rm = TRUE) +
                      (4.5* sd(gene_expression[counts_plus$annotated_sex == 'female'], na.rm = TRUE))
    low_threshold <- mean(gene_expression[counts_plus$annotated_sex == 'female'], na.rm = TRUE)+
                     (3 * sd(gene_expression[counts_plus$annotated_sex == 'female'], na.rm = TRUE))
  }
  
  # Predict sex based on gene expression thresholds
 if (gene %in% male_genes) {
    # For Y chromosome associated genes, higher expression suggests "male"
    predicted_sex <- ifelse(gene_expression >= high_threshold, "male",
                            ifelse(gene_expression <= low_threshold, "female", "cannot_predict"))
  } else {
    # For other genes (assuming X chromosome), higher expression suggests "female"
    predicted_sex <- ifelse(gene_expression >= high_threshold, "female",
                            ifelse(gene_expression <= low_threshold, "male", "cannot_predict"))
  }
  # Prepare data for plotting
  gene_data_for_plot <- data.frame(annotated_sex = counts_plus$annotated_sex, expression = gene_expression)
  
  p <-ggplot(gene_data_for_plot, aes(x = annotated_sex, y = expression, fill = annotated_sex)) +
      geom_violin(trim = FALSE) + 
      scale_fill_manual(values=c("grey", "orange"))+ 
      geom_jitter(size = 0.75) + 
      ylab(paste0(gene," Expression (log counts)")) + 
      geom_hline(yintercept = high_threshold, linetype="dashed", color = "maroon") + 
      geom_hline(yintercept = low_threshold, linetype="dashed", color = "blue") + 
      theme_light() # Change scale if you want log10
  print(p)
  
   comparison_reported_predicted <- table(counts_plus$annotated_sex, predicted_sex)

    # Visualize the comparison
   x <- barplot(comparison_reported_predicted, 
            beside = TRUE, 
            legend.text = TRUE,
            main = paste(gene, "gene for predicting sex in TCGA-LIHC"),
            xlab = "Reported Sex", 
            ylab = "Number of Samples",
            col = c("grey", "orange"),
            args.legend = list(x = "topleft", bty = "n", inset = c(0.05, 0))) 
    y <- as.matrix(comparison_reported_predicted)

text(x, y + 20, labels = as.character(y))
  
  
  # Add expression data to the dataframe
  all_genes_data[[paste0(gene, "_expression")]] <- counts_plus[[gene]]
  all_genes_data[[paste0(gene, "_predicted_sex")]] <- predicted_sex
  all_genes_data[[paste0(gene, "_predicted_expression")]] <- ifelse(counts_plus[[gene]] >= high_threshold, "high_expression", ifelse(counts_plus[[gene]] <= low_threshold, "low_expression", "intermediate_expression"))
    # Comparison of reported_sex to predicted_sex
  
#   
# violin6_fname <- paste("Log_Violin", gene, ".png", sep = "_")
# ggsave(violin6_fname, p, width = 10, height = 8, units = "in", dpi = 300)
}

# Save the comprehensive results to a CSV file
write.csv(all_genes_data, "expression_data_all_genes_with_predicted_sex.csv", row.names = FALSE)
```
```{r}
library(ggplot2)

# Load your data
counts_plus <- read.table(file = "C:/Users/livin/ASU FILES/Reseach Files/Genomic Sex Chromosome/TCGA_RNAseq_counts/TCGA_RNA_data_all_modified.csv", sep = '\t', header = TRUE)

# Define female and male genes
female_genes <- c("XIST")
male_genes <- c("USP9Y", "UTY", "ZFY", "SRY", "AMELY", "DDX3Y", "EIF1AY", "KDM5D", "NLGN4Y", "PRKY", "TMSB4Y")

# Initialize a dataframe to store results for all genes
all_genes_data <- data.frame(RowID  = rownames(counts_plus), 
                             annotated_sex = counts_plus$annotated_sex, 
                             XIST_status = counts_plus$status_XIST, 
                             Y_status = counts_plus$status_Y)

for (gene in c(female_genes, male_genes)) {
    gene_expression <- log1p(counts_plus[[gene]]) 
    gene_expression[is.infinite(gene_expression)] <- NA

    # Determine thresholds based on the sex
    if (gene %in% female_genes) {
        mean_expr <- mean(counts_plus[[gene]][counts_plus$annotated_sex == 'male'], na.rm = TRUE)
        sd_expr <- sd(counts_plus[[gene]][counts_plus$annotated_sex == 'male'], na.rm = TRUE)
    } else {
        mean_expr <- mean(counts_plus[[gene]][counts_plus$annotated_sex == 'female'], na.rm = TRUE)
        sd_expr <- sd(counts_plus[[gene]][counts_plus$annotated_sex == 'female'], na.rm = TRUE)
    }
    
    high_threshold <- mean_expr + 4.5 * sd_expr
    low_threshold <- mean_expr + 3 * sd_expr

    # Prepare data for plotting
    gene_data_for_plot <- data.frame(annotated_sex = counts_plus$annotated_sex, expression = gene_expression)

    # Plotting with specified colors and thresholds
    p <-ggplot(gene_data_for_plot, aes(x = annotated_sex, y = expression, fill = annotated_sex)) +
      geom_violin(trim = FALSE) + 
      scale_fill_manual(values=c("grey", "orange"))+ 
      geom_jitter(size = 0.75) + 
      ylab(paste0(gene," Expression (log counts)")) + 
      scale_y_continuous(trans = "log10") +
      # geom_hline(yintercept = high_threshold, linetype="dashed", color = "maroon") +
      # geom_hline(yintercept = low_threshold, linetype="dashed", color = "blue") +
      geom_hline(yintercept = mean_expr, color = "purple", linetype = "solid") +
      geom_hline(yintercept = sd_expr, color = "darkgreen", linetype = "dashed") +
      #geom_hline(yintercept = mean_expr - sd_expr, color = "darkgreen", linetype = "dashed") +
      labs(y = paste(gene, " Expression (log-transformed)"), title = paste(gene, " Expression by Annotated Sex")) +
    theme_minimal()

# Printing and saving the plot
print(p)
# ggsave(paste0("violin_plot_", gene, ".png"), plot = p, width = 10, height = 8, dpi = 300)

   # Add expression data to the dataframe
  all_genes_data[[paste0(gene, "_expression")]] <- counts_plus[[gene]]
  all_genes_data[[paste0(gene, "_predicted_sex")]] <- predicted_sex
  all_genes_data[[paste0(gene, "_predicted_expression")]] <- ifelse(counts_plus[[gene]] >= high_threshold, "high_expression", ifelse(counts_plus[[gene]] <= low_threshold, "low_expression", "intermediate_expression"))
    # Comparison of reported_sex to predicted_sex
  
#   
violin6_fname <- paste("violin_plot_", gene, ".png", sep = "_")
ggsave(violin6_fname, p, width = 10, height = 8, units = "in", dpi = 300)
}

# Save the comprehensive results
write.csv(all_genes_data, "expression_data_all_genes_with_predicted_sex.csv", row.names = FALSE)



```

```{r}
# Load your data here
library(ggplot2)
counts_plus <- read_csv("TCGA_RNA_data_all_modified.csv")
# Load required library


# Define female and male genes
female_genes <- c("XIST")
male_genes <- c("USP9Y", "UTY", "ZFY", "SRY", "AMELY", "DDX3Y", "EIF1AY", "KDM5D", "NLGN4Y", "PRKY", "TMSB4Y")

# Initialize a dataframe to store results for all genes
all_genes_data <- data.frame(RowID = rownames(counts_plus), 
                             annotated_sex = counts_plus$sex, 
                             XIST_status = counts_plus$status_XIST, 
                             Y_status = counts_plus$status_Y)

for (gene in c(female_genes, male_genes)) {
  # Compute gene expression with log1p transformation
  gene_expression <- log1p(counts_plus[[gene]])
  gene_expression[is.infinite(gene_expression)] <- NA

  # Set thresholds based on the sex
  if (gene %in% female_genes) {
    high_threshold <- mean(counts_plus[[gene]][counts_plus$annotated_sex == 'male'], na.rm = TRUE) + 
                      (2.5 * sd(counts_plus[[gene]][counts_plus$annotated_sex == 'male'], na.rm = TRUE))
    low_threshold <- mean(counts_plus[[gene]][counts_plus$annotated_sex == 'male'], na.rm = TRUE) + 
                     (3 * sd(counts_plus[[gene]][counts_plus$annotated_sex == 'male'], na.rm = TRUE))
  } else {
    high_threshold <- mean(counts_plus[[gene]][counts_plus$annotated_sex == 'female'], na.rm = TRUE) + 
                      (4.5 * sd(counts_plus[[gene]][counts_plus$annotated_sex == 'female'], na.rm = TRUE))
    low_threshold <- mean(counts_plus[[gene]][counts_plus$annotated_sex == 'female'], na.rm = TRUE) + 
                     (3 * sd(counts_plus[[gene]][counts_plus$annotated_sex == 'female'], na.rm = TRUE))
  }
  
  # Predict sex based on gene expression thresholds
  if (gene %in% male_genes) {
    predicted_sex <- ifelse(gene_expression >= high_threshold, "male",
                            ifelse(gene_expression <= low_threshold, "female", "cannot_predict"))
  } else {
    predicted_sex <- ifelse(gene_expression >= high_threshold, "female",
                            ifelse(gene_expression <= low_threshold, "male", "cannot_predict"))
  }
  
  # Prepare data for plotting
  gene_data_for_plot <- data.frame(annotated_sex = counts_plus$annotated_sex, expression = gene_expression)

  p <- ggplot(gene_data_for_plot, aes(x = annotated_sex, y = expression, fill = annotated_sex)) +
    geom_violin(trim = FALSE) +
    scale_y_continuous(trans = "log10") +
    scale_fill_manual(values = c("grey", "orange")) +
    geom_jitter(size = 0.75) +
    ylab(paste0(gene, " Expression (log1p-transformed)")) +
    geom_hline(yintercept = high_threshold, linetype="dashed", color = "maroon") +
    geom_hline(yintercept = low_threshold, linetype="dashed", color = "blue") +
    theme_light()
  
  print(p)
  
  # Save the plot
  violin6_fname <- paste("Log_Violin", gene, ".png", sep = "_")
  ggsave(violin6_fname, p, width = 10, height = 8, units = "in", dpi = 300)
  
  # Add expression data to the dataframe
  all_genes_data[[paste0(gene, "_log_expression")]] <- gene_expression
  all_genes_data[[paste0(gene, "_predicted_sex")]] <- predicted_sex
  all_genes_data[[paste0(gene, "_predicted_expression")]] <- ifelse(gene_expression >= high_threshold, "high_expression", 
                                                                     ifelse(gene_expression <= low_threshold, "low_expression", "intermediate_expression"))
  
  # Compare reported sex to predicted sex
  comparison_reported_predicted <- table(counts_plus$annotated_sex, predicted_sex)
  
  # Visualize the comparison
  x <- barplot(comparison_reported_predicted, 
               beside = TRUE, 
               legend.text = TRUE,
               main = paste(gene, "gene for predicting sex in TCGA-LIHC"),
               xlab = "Reported Sex", 
               ylab = "Number of Samples",
               col = c("grey", "orange"),
               args.legend = list(x = "topleft", bty = "n", inset = c(0.05, 0)))
  
  y <- as.matrix(comparison_reported_predicted)
  text(x, y + 20, labels = as.character(y))
}

# Save the comprehensive results to a CSV file
write.csv(all_genes_data, "expression_data_all_genes_with_predicted_sex.csv", row.names = FALSE)

```

