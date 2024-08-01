library(tidyverse)
library(data.table)
library(ggplot2)
library(plotly)

# List of TCGA studies to analyze
cancer_types <- c("LAML", "ACC", "BLCA", "LGG", "BRCA", "CESC", "CHOL", "COAD", "ESCA", 
                  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "DLBC", 
                  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
                  "TGCT", "THCA", "THYM", "UCS", "UCEC", "UVM")

# Set directories
working_dir <- 'C:/Users/livin/ASU FILES/Reseach Files/Genomic Sex Chromosome/TCGA_RNAseq_counts'
setwd(working_dir)
filepath <- "C:/Users/livin/ASU FILES/Reseach Files/Genomic Sex Chromosome/TCGA_RNAseq_counts/"

# Genes under consideration for inference:
ychr_genes <- c("ENSG00000099721", "ENSG00000067048", "ENSG00000198692",
                "ENSG00000012817", "ENSG00000165246", "ENSG00000099725",
                "ENSG00000154620", "ENSG00000114374", "ENSG00000183878",
                "ENSG00000067646", "ENSG00000184895","ENSG00000129824")

ychr_gnames <- c("AMELY", "DDX3Y", "EIF1AY", "KDM5D", "NLGN4Y", "PRKY",
                 "TMSB4Y", "USP9Y", "UTY", "ZFY", "SRY","RPS4Y1")

xchr_genes  <- c("ENSG00000229807", "ENSG00000169083")
xchr_gnames <- c("XIST", "AR")

for (study in cancer_types) {
  # Define file paths
  counts_fname <- paste("TCGA", study, "TPM.tsv", sep = "-")
  meta_fname <- paste("TCGA", study, "META.tsv", sep = "-")
  
  counts_path <- paste(filepath, counts_fname, sep = "")
  meta_path <- paste(filepath, meta_fname, sep = "")
  
  # Check if the files exist
  if (!file.exists(counts_path) | !file.exists(meta_path)) {
    cat("File not found for study:", study, "\n")
    next
  }
  
  # Read in data
  counts <- read.delim(counts_path, row.names = 1)
  metadf <- read.delim(meta_path, row.names = 1)
  
  # Modify row id's in metadf to correspond to similar format in counts
  meta_ids <- rownames(metadf)
  meta_ids <- gsub("[-]", ".", meta_ids)
  rownames(metadf) <- meta_ids
  
  ychr_counts <- counts[ychr_genes, ]
  xchr_counts <- counts[xchr_genes, ]
  
  new_counts1 <- data.frame(t(subset(xchr_counts, select = -c(gene_name, gene_type)))) #xchr genes 
  colnames(new_counts1) <- xchr_counts[, 1]
  
  new_counts2 <- data.frame(t(subset(ychr_counts, select = -c(gene_name, gene_type)))) #ychr genes
  colnames(new_counts2) <- ychr_counts[, 1]
  
  new_counts <- data.frame(c(new_counts1, new_counts2)) #combine new_count1(xchr genes) and new_count2 (ychr)
  row.names(new_counts) <- row.names(new_counts1)
  
  ################################################
  ###---CHECK FOR REPEATS AND LOOK AT VALUES---###
  ################################################
  
  # Return list of row labels (ids):
  id_list <- rownames(new_counts)
  
  # Extract file and case uuid's from combo:
  splt_id_list <- strsplit(id_list, "[_]")
  c_uuid <- vector(mode = "character", length = length(id_list))
  f_uuid <- vector(mode = "character", length = length(id_list))
  
  for (j in seq_len(length(id_list))) {
    f_uuid[j] <- splt_id_list[[j]][1]
    c_uuid[j] <- splt_id_list[[j]][2]
  }
  
  # Organize ids by case uuid (for viewing)
  id_sets <- split(id_list, c_uuid)
  
  repeats_df <- data.frame()
  rep_cases_list <- list()
  
  c_uuid <- unique(c_uuid)
  
  for (x in c_uuid) {
    if (length(id_sets[[x]]) > 1) {
      
      # extract labels needed for the next step
      rows_keep <- id_sets[[x]]
      
      # PUTTING THESE INTO A SEPARATE DF TO *LOOK* AT:
      repeats_df <- bind_rows(repeats_df, new_counts[rows_keep[1], ])
      repeats_df <- bind_rows(repeats_df, new_counts[rows_keep[2], ])
      
      if (length(id_sets[[x]]) > 2) {
        repeats_df <- bind_rows(repeats_df, new_counts[rows_keep[3], ])
      }
      
      rep_cases_list <- c(rep_cases_list, x)
    }
  }
  
  # repeats_df dataframe to tsv
  repeats_fname <- paste("TCGA", study, "repeat_samples_TPM_counts.tsv", sep = "_")
  write_tsv(
    repeats_df %>% rownames_to_column(),
    repeats_fname,
    na = "NA",
    append = FALSE,
    col_names = TRUE,
    quote = "none",
    eol = "\n",
    num_threads = readr_threads(),
    progress = show_progress()
  )
  
  ################################################
  ###---ANALYZE DATA AND CREATE FINAL DATAFRAMES---###
  ################################################
  
  # Create the sex_check dataframe
  sex_check <- data.frame(matrix("no", length(c_uuid), 12), row.names = c_uuid)
  colnames(sex_check) <- c("status_XIST", "status_Y",
                           "annotated_sex", "survival_status", "time_to_status", "age",
                           "Solid_Tissue_Normal", "Blood_Derived_Normal", "Primary_Tumor", "Recurrent_Tumor", "Primary Blood Derived Cancer - Peripheral Blood", "sample_types")
  
  # Create the counts_plus dataframe
  dummy_df <- data.frame(matrix("no", length(id_list), 12))
  counts_plus <- data.frame(c(new_counts, dummy_df))
  row.names(counts_plus) <- row.names(new_counts)
  colnames(counts_plus) <- c(colnames(new_counts), "status_XIST", "status_Y",
                             "annotated_sex", "survival", "followed", "age",
                             "Solid_Tissue_Normal", "Blood_Derived_Normal", "Primary_Tumor", "Recurrent_Tumor", "Primary Blood Derived Cancer - Peripheral Blood", "sample_types")
  
  # Process each case
  for (i in c_uuid) {
    # Set annotated_sex based on metadata
    sex_check[i, "annotated_sex"] <- metadf[i, "mf_list"]
    
    # Set survival_status and time_to_status based on status_list and follow_list/surv_times
    if (!is.na(metadf[i, "status_list"]) && metadf[i, "status_list"] == "Alive") {
      sex_check[i, "survival_status"] <- 0
      sex_check[i, "time_to_status"] <- metadf[i, "follow_list"]
    } else {
      sex_check[i, "survival_status"] <- 1
      sex_check[i, "time_to_status"] <- metadf[i, "surv_times"]
    }
    
    # Extract sample types for the current case
    sample_types <- strsplit(metadf[i, "sample_types"], ",")[[1]]
    
    # Set age based on metadata
    sex_check[i, "age"] <- metadf[i, "age_index"]
    
    # Handle replicates
    if (length(id_sets[[i]]) > 1) {
      iid <- id_sets[[i]]  # Extract all ids for the current case
      counts_plus[iid, "annotated_sex"] <- metadf[i, "mf_list"]
      counts_plus[iid, c("survival", "followed")] <- metadf[i, c("surv_times", "follow_list")]
      counts_plus[iid, "age"] <- metadf[i, "age_index"]
      counts_plus[iid, "sample_types"] <- paste(sample_types, collapse=",")
      sex_check[i, "sample_types"] <- paste(sample_types, collapse=",")
      
      # Initialize sample type columns
      counts_plus[iid, c("Solid_Tissue_Normal", "Blood_Derived_Normal", "Primary_Tumor", "Recurrent_Tumor", "Primary Blood Derived Cancer - Peripheral Blood")] <- "no"
      sex_check[i, c("Solid_Tissue_Normal", "Blood_Derived_Normal", "Primary_Tumor", "Recurrent_Tumor", "Primary Blood Derived Cancer - Peripheral Blood")] <- "no"
      
      for (j in seq_along(iid)) {
        for (sample_type in sample_types) {
          # Update sample types
          if (sample_type == "Solid Tissue Normal") {
            counts_plus[iid[j], "Solid_Tissue_Normal"] <- "yes"
            sex_check[i, "Solid_Tissue_Normal"] <- "yes"
          }
          if (sample_type == "Blood Derived Normal") {
            counts_plus[iid[j], "Blood_Derived_Normal"] <- "yes"
            sex_check[i, "Blood_Derived_Normal"] <- "yes"
          }
          if (sample_type == "Primary Tumor") {
            counts_plus[iid[j], "Primary_Tumor"] <- "yes"
            sex_check[i, "Primary_Tumor"] <- "yes"
          }
          if (sample_type == "Recurrent Tumor") {
            counts_plus[iid[j], "Recurrent_Tumor"] <- "yes"
            sex_check[i, "Recurrent_Tumor"] <- "yes"
          }
          if (sample_type == "Primary Blood Derived Cancer - Peripheral Blood") {
            counts_plus[iid[j], "Primary Blood Derived Cancer - Peripheral Blood"] <- "yes"
            sex_check[i, "Primary Blood Derived Cancer - Peripheral Blood"] <- "yes"
          }
        }
        
        # -- Y chromosome:
        if (all(new_counts[iid[j], c("DDX3Y", "USP9Y", "UTY", "ZFY")] < 1.0)) {
          counts_plus[iid[j], "status_Y"] <- "no"
        } else {
          counts_plus[iid[j], "status_Y"] <- "yes"
        }
        
        # -- XIST:
        counts_plus[iid[j], "status_XIST"] <- ifelse(new_counts[iid[j], "XIST"] > 1.0, "yes", "no")
      }
      
      # Evaluate the *pairs* of replicates (only written for up to 3 reps!)
      if (length(iid) == 2) {
        for (k in 1:length(iid)) {
          sex_check[i, "status_Y"] <- counts_plus[iid[k], "status_Y"]
          sex_check[i, "status_XIST"] <- counts_plus[iid[k], "status_XIST"]
          
          # Check sample types for two replicates
          for (sample_type in c("Solid_Tissue_Normal", "Blood_Derived_Normal", "Primary_Tumor", "Recurrent_Tumor", "Primary Blood Derived Cancer - Peripheral Blood")) {
            if (counts_plus[iid[k], sample_type] == "yes") {
              sex_check[i, sample_type] <- "yes"
            }
          }
        }
      } else if (length(iid) == 3) {
        Y_status <- counts_plus[iid, "status_Y"]
        XIST_status <- counts_plus[iid, "status_XIST"]
        sex_check[i, "status_Y"] <- if (all(Y_status == "yes") || all(Y_status == "no")) Y_status[1] else "ambiguous"
        sex_check[i, "status_XIST"] <- if (all(XIST_status == "yes") || all(XIST_status == "no")) XIST_status[1] else "ambiguous"
        
        # Check sample types for three replicates
        for (sample_type in c("Solid_Tissue_Normal", "Blood_Derived_Normal", "Primary_Tumor", "Recurrent_Tumor", "Primary Blood Derived Cancer - Peripheral Blood")) {
          sample_type_status <- counts_plus[iid, sample_type]
          sex_check[i, sample_type] <- if (all(sample_type_status == "yes") || all(sample_type_status == "no")) sample_type_status[1] else "ambiguous"
        }
      }
    } else {
      # Handle cases with only a single sample
      iid <- id_sets[[i]]
      counts_plus[iid, "annotated_sex"] <- metadf[i, "mf_list"]
      counts_plus[iid, c("survival", "followed")] <- metadf[i, c("surv_times", "follow_list")]
      counts_plus[iid, "age"] <- metadf[i, "age_index"]
      counts_plus[iid, "sample_types"] <- paste(sample_types, collapse=",")
      sex_check[i, "sample_types"] <- paste(sample_types, collapse=",")
      # Initialize sample type columns to "no"
      counts_plus[iid, c("Solid_Tissue_Normal", "Blood_Derived_Normal", "Primary_Tumor", "Recurrent_Tumor", "Primary Blood Derived Cancer - Peripheral Blood")] <- "no"
      sex_check[i, c("Solid_Tissue_Normal", "Blood_Derived_Normal", "Primary_Tumor", "Recurrent_Tumor", "Primary Blood Derived Cancer - Peripheral Blood")] <- "no"
      
      # Split sample_types into a vector and process each one
      sample_types_vector <- unlist(strsplit(sample_types, ",\\s*"))
      
      for (sample_type in sample_types_vector) {
        sample_type <- trimws(sample_type)  # Standardize and trim whitespace
        # print(paste("Processing sample type:", sample_type))  # Debugging output
        
        if (sample_type == "Solid Tissue Normal") {
          counts_plus[iid, "Solid_Tissue_Normal"] <- "yes"
          sex_check[i, "Solid_Tissue_Normal"] <- "yes"
        } else if (sample_type == "Blood Derived Normal") {
          counts_plus[iid, "Blood_Derived_Normal"] <- "yes"
          sex_check[i, "Blood_Derived_Normal"] <- "yes"
        } else if (sample_type == "Primary Tumor") {
          counts_plus[iid, "Primary_Tumor"] <- "yes"
          sex_check[i, "Primary_Tumor"] <- "yes"
        } else if (sample_type == "Recurrent Tumor") {
          counts_plus[iid, "Recurrent_Tumor"] <- "yes"
          sex_check[i, "Recurrent_Tumor"] <- "yes"
        } else if (sample_type == "Primary Blood Derived Cancer - Peripheral Blood") {
          counts_plus[iid, "Primary Blood Derived Cancer - Peripheral Blood"] <- "yes"
          sex_check[i, "Primary Blood Derived Cancer - Peripheral Blood"] <- "yes"
        }
      }
      
      # -- Y chromosome:
      if (all(new_counts[iid, c("DDX3Y", "USP9Y", "UTY", "ZFY")] < 1.0)) {
        counts_plus[iid, "status_Y"] <- "no"
        sex_check[i, "status_Y"] <- "no"
      } else {
        counts_plus[iid, "status_Y"] <- "yes"
        sex_check[i, "status_Y"] <- "yes"
      }
      # -- XIST:
      counts_plus[iid, "status_XIST"] <- ifelse(new_counts[iid, "XIST"] > 1.0, "yes", "no")
      sex_check[i, "status_XIST"] <- counts_plus[iid, "status_XIST"]
    }
  }
  
  
  # Write sex_check dataframe to tsv
  sex_check_fname <- paste("TCGA", study, "case_XIST-Y_outcomes.tsv", sep = "_")
  write_tsv(sex_check %>% rownames_to_column(), sex_check_fname, na = "NA", append = FALSE, col_names = TRUE, quote = "none", eol = "\n", num_threads = readr_threads(), progress = show_progress())
  
  # Write counts_plus dataframe to tsv
  counts_plus_fname <- paste("TCGA", study, "sample_XIST-Y_outcomes.tsv", sep = "_")
  write_tsv(counts_plus %>% rownames_to_column(), counts_plus_fname, na = "NA", append = FALSE, col_names = TRUE, quote = "none", eol = "\n", num_threads = readr_threads(), progress = show_progress())
  
  cat("Processed study:", study, "\n")
}

cat("All studies processed.\n")
