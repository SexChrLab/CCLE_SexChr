
title: "Sex-check-analysis-loop"
author: "Sussan Massey and edited by Ilsa Rodriguez "
date: "`r Sys.Date()`"


# Import libraries to use
library(tidyverse)
library(ggplot2)
library(plotly)

# # List of TCGA studies to analyze

#   "LAML", "ACC", "BLCA", "LGG", "BRCA", "CESC", "CHOL", "COAD", "ESCA", 
#   "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "DLBC", 
#   "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
#   "TGCT", "THCA", "THYM", "UCS", "UCEC", "UVM"

study <- 'PAAD' #call which cancer study to want to work with 

# Set directories( set to your working directories!! )
setwd('C:/Users/livin/ASU FILES/Reseach Files/Genomic Sex Chromosome/TCGA_RNAseq_counts')
filepath <- "C:/Users/livin/ASU FILES/Reseach Files/Genomic Sex Chromosome/TCGA_RNAseq_counts/"


#created count files and calls them into your R enviroment

counts_fname <- paste("TCGA", study, "TPM.tsv", sep = "-")
meta_fname <- paste("TCGA", study, "META.tsv", sep = "-")

counts_path <- paste(filepath, counts_fname, sep = "")
meta_path <- paste(filepath, meta_fname, sep = "")

# Read in data
counts <- read.delim(counts_path, row.names = 1)
metadf <- read.delim(meta_path, row.names = 1)

# Modify row id's in metadf to correspond to similar format in counts
meta_ids <- rownames(metadf)
meta_ids <- gsub("[-]", ".", meta_ids)
rownames(metadf) <- meta_ids




# Genes under consideration for inference:
#   chromosome Y 
#     AMELY  ENSG00000099721
#     DDX3Y  ENSG00000067048
#     EIF1AY ENSG00000198692
#     KDM5D  ENSG00000012817
#     NLGN4Y ENSG00000165246
#     PRKY   ENSG00000099725
#     TMSB4Y ENSG00000154620
#     USP9Y  ENSG00000114374
#     UTY    ENSG00000183878
#     ZFY    ENSG00000067646
#     SRY    ENSG00000184895
#     TSPY   (may be challenging due to multicopy, there are 10+, left out)
#
#   chromosome X
#     XIST   ENSG00000229807
#     AR     ENSG00000169083

ychr_genes <- c("ENSG00000099721", "ENSG00000067048", "ENSG00000198692",
                "ENSG00000012817", "ENSG00000165246", "ENSG00000099725",
                "ENSG00000154620", "ENSG00000114374", "ENSG00000183878",
                "ENSG00000067646", "ENSG00000184895")

ychr_gnames <- c("AMELY", "DDX3Y", "EIF1AY", "KDM5D", "NLGN4Y", "PRKY",
                 "TMSB4Y", "USP9Y", "UTY", "ZFY", "SRY")

xchr_genes  <- c("ENSG00000229807", "ENSG00000169083")
xchr_gnames <- c("XIST", "AR")

ychr_counts <- counts[ychr_genes, ]
xchr_counts <- counts[xchr_genes, ]

## may want a version of new_counts for x and y?:
new_counts1 <- data.frame(t(subset(xchr_counts,
                                   select = -c(gene_name, gene_type)))) #xchr genes 

colnames(new_counts1) <- xchr_counts[, 1]

new_counts2 <- data.frame(t(subset(ychr_counts,
                                   select = -c(gene_name, gene_type)))) #ychr genes
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



sex_check <- data.frame(matrix(0, length(c_uuid), 5),
                        row.names = c_uuid)
colnames(sex_check) <- c("status_XIST", "status_Y", "annotated_sex",
                         "survival_status", "time_to_status")

dummy_df <- data.frame(matrix(0, length(id_list), 5))

counts_plus <- data.frame(c(new_counts, dummy_df))
row.names(counts_plus) <- row.names(new_counts)
colnames(counts_plus) <- c(colnames(new_counts), "status_XIST", "status_Y",
                           "annotated_sex", "survival", "followed")



for (i in c_uuid) {
  # Set annotated_sex based on metadata
  sex_check[i, "annotated_sex"] <- metadf[i, "mf_list"]
  # Now when you do further operations with counts_plus, it won't include rows where annotated_sex is NA
  # Set survival_status and time_to_status based on status_list and follow_list/surv_times
  if (!is.na(metadf[i, "status_list"]) && metadf[i, "status_list"] == "Alive") {
    sex_check[i, "survival_status"] <- 0
    sex_check[i, "time_to_status"] <- metadf[i, "follow_list"]
  } else {
    sex_check[i, "survival_status"] <- 1
    sex_check[i, "time_to_status"] <- metadf[i, "surv_times"]
  }
  
  # Check if there are replicates and if so, handle those first:
  if (length(id_sets[[i]]) > 1) {
    iid <- id_sets[[i]]  # Extract all ids for the current case
    counts_plus[iid, "annotated_sex"] <- metadf[i, "mf_list"]
    counts_plus[iid, c("survival", "followed")] <- metadf[i, c("surv_times", "follow_list")]
    
    # Loop through each sample for the current case
    for (j in seq_along(iid)) {
      # -- Y chromosome:
      if (all(new_counts[iid[j], c("DDX3Y", "USP9Y", "UTY", "ZFY")] < 1.0)) {
        counts_plus[iid[j], "status_Y"] <- "no"
      } else {
        counts_plus[iid[j], "status_Y"] <- "yes"
      }
      
      # -- XIST:
      counts_plus[iid[j], "status_XIST"] <- ifelse(new_counts[iid[j], "XIST"] > 1.0, "yes", "no")
    }
    
    # Evaluate the *pairs* of replicates, should handle all pairs (only written for up to 3 replicates!) - Ilsa
    if (length(iid) == 2) {
      # For two samples, evaluate the pairs directly
      for (k in 1:length(iid)) {
        sex_check[i, "status_Y"] <- counts_plus[iid[k], "status_Y"]
        sex_check[i, "status_XIST"] <- counts_plus[iid[k], "status_XIST"]
      }
    } else if (length(iid) == 3) {
      # For three samples, evaluate based on majority vote
      Y_status <- counts_plus[iid, "status_Y"]
      XIST_status <- counts_plus[iid, "status_XIST"]
      if (all(Y_status == "yes") || all(Y_status == "no")) {
        sex_check[i, "status_Y"] <- Y_status[1]  # Take the first one as majority
      } else {
        sex_check[i, "status_Y"] <- "ambiguous"
      }
      if (all(XIST_status == "yes") || all(XIST_status == "no")) {
        sex_check[i, "status_XIST"] <- XIST_status[1]  
      } else {
        sex_check[i, "status_XIST"] <- "ambiguous"
      }
    }
    
  } else {
    # Handle cases with only a single sample
    iid <- id_sets[[i]]
    counts_plus[iid, "annotated_sex"] <- metadf[i, "mf_list"]
    counts_plus[iid, c("survival", "followed")] <- metadf[i, c("surv_times", "follow_list")]
    
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


###################################################
###---WRITE DATAFRAMES TO TSV FILES---###
###################################################

# sex_check dataframe to tsv
sex_check_fname <- paste("TCGA", study, "case_XIST-Y_outcomes.tsv", sep = "_")
write_tsv(
  sex_check %>% rownames_to_column(),
  sex_check_fname,
  na = "NA",
  append = FALSE,
  col_names = TRUE,
  quote = "none",
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress()
)

# counts_plus dataframe to tsv
counts_plus_fname <- paste("TCGA", study, "sample_XIST-Y_outcomes.tsv", sep = "_")

write_tsv(
  counts_plus %>% rownames_to_column(),
  counts_plus_fname,
  na = "NA",
  append = FALSE,
  col_names = TRUE,
  quote = "none",
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress()
)



###################################################
###---SURVIVAL ANALYSIS---###
###################################################
library(dplyr)
library(survival)
library(survminer)

sex_check_m <- sex_check %>% dplyr::filter(annotated_sex == "male")
df_points2 <- sex_check_m %>%
  transmute(XIST_Y = paste(status_XIST, "_", status_Y, sep = ""))
sex_check_m <- cbind(sex_check_m, df_points2)

sex_check_f <- sex_check %>% dplyr::filter(annotated_sex == "female")
df_points3 <- sex_check_f %>%
  transmute(XIST_Y = paste(status_XIST, "_", status_Y, sep = ""))
sex_check_f <- cbind(sex_check_f, df_points3)


################################################################################


km_m <- survfit(Surv(time_to_status, survival_status) ~ XIST_Y, data = sex_check_m)
km_f <- survfit(Surv(time_to_status, survival_status) ~ XIST_Y, data = sex_check_f)


km_m_plot_fname <- paste("KM", study, "Male.png", sep = "_")
png(km_m_plot_fname)

km_m %>%
  ggsurvplot(
    data = sex_check_m,
    fun = "pct",
    # linetype = "strata", # Change line type by groups
    # pval = TRUE, # Not sure if want
    # conf.int = TRUE, # Not sure if want
    risk.table = TRUE,
    fontsize = 3, # used in risk table
    surv.median.line = "hv", # median horizontal and vertical ref lines
    ggtheme = theme_light(),
    palette = c("goldenrod", "sienna", "tomato", "cadetblue", "dodgerblue","blue"),
    title = " Male - Kaplan-Meier Survival Function Estimate",
    legend.title = "",
    legend.labs = levels(sex_check_m$XIST_Y)
  )
dev.off()

km_f_plot_fname <- paste("KM", study, "Female.png", sep = "_")
png(km_f_plot_fname)

km_f %>%
  ggsurvplot(
    data = sex_check_f,
    fun = "pct",
    # linetype = "strata", # Change line type by groups
    # pval = TRUE, # Not sure if want
    # conf.int = TRUE, # Not sure if want
    risk.table = TRUE,
    fontsize = 3, # used in risk table
    surv.median.line = "hv", # median horizontal and vertical ref lines
    ggtheme = theme_light(),
    palette = c("goldenrod", "sienna", "tomato", "cadetblue", "dodgerblue","blue"),
    title = "Female - Kaplan-Meier Survival Function Estimate",
    legend.title = "",
    legend.labs = levels(sex_check_f$XIST_Y)
  )
dev.off()




counts_plus <- counts_plus %>% filter(!is.na(annotated_sex))

###################################################
###---COUNT GROUP MEMBERSHIP FOR TABLE---###
###################################################

# among those annotated "male":
n_male_nxny <- nrow(sex_check_m[sex_check_m$XIST_Y == "noXIST_noY", ])
n_male_yxyy <- nrow(sex_check_m[sex_check_m$XIST_Y == "yesXIST_yesY", ])
n_male_nxyy <- nrow(sex_check_m[sex_check_m$XIST_Y == "noXIST_yesY", ])
n_male_yxny <- nrow(sex_check_m[sex_check_m$XIST_Y == "yesXIST_noY", ])

# among those annotated "female":
n_female_nxny <- nrow(sex_check_f[sex_check_f$XIST_Y == "noXIST_noY", ])
n_female_yxyy <- nrow(sex_check_f[sex_check_f$XIST_Y == "yesXIST_yesY", ])
n_female_nxyy <- nrow(sex_check_f[sex_check_f$XIST_Y == "noXIST_yesY", ])
n_female_yxny <- nrow(sex_check_f[sex_check_f$XIST_Y == "yesXIST_noY", ])

#####################################################
###-VIOLIN PLOTS OF TPM COUNTS OF INTEREST BY M/F-###
###-TO SHOW IF THRESHOLD WAS CHOSEN APPROPRIATELY-###
#####################################################
library(patchwork)


pXIST <- ggplot(data = counts_plus, aes(factor(annotated_sex), XIST))
pXIST <- pXIST + geom_violin() + scale_y_continuous(trans = "log10") +
  geom_jitter(height = 0, width = 0.1)

pDDX3Y <- ggplot(data = counts_plus, aes(factor(annotated_sex), DDX3Y))
pDDX3Y <- pDDX3Y + geom_violin() + scale_y_continuous(trans = "log10") +
  geom_jitter(height = 0, width = 0.1)

pUSP9Y <- ggplot(data = counts_plus, aes(factor(annotated_sex), USP9Y))
pUSP9Y <- pUSP9Y + geom_violin() + scale_y_continuous(trans = "log10") +
  geom_jitter(height = 0, width = 0.1)

pUTY <- ggplot(data = counts_plus, aes(factor(annotated_sex), UTY))
pUTY <- pUTY + geom_violin() + scale_y_continuous(trans = "log10") +
  geom_jitter(height = 0, width = 0.1)

pZFY <- ggplot(data = counts_plus, aes(factor(annotated_sex), ZFY))
pZFY <- pZFY + geom_violin() + scale_y_continuous(trans = "log10") +
  geom_jitter(height = 0, width = 0.1)



p_all <-pXIST | (pDDX3Y | pUSP9Y) / (pUTY | pZFY)
p_all

# Save the composite plot
violin6_fname <- paste("Log_Violin", study, "TPM_All.png", sep = "_")
ggsave(violin6_fname, p_all, width = 10, height = 8, units = "in", dpi = 300)



pSRY <- ggplot(data = counts_plus, aes(factor(annotated_sex), SRY))
pSRY <- pSRY + geom_violin() + scale_y_continuous(trans = "log10") +
  geom_jitter(height = 0, width = 0.1)
pSRY



# Violins in plotly for interactive viewing
figDDX3Y <- counts_plus %>%
  plot_ly(
    x = ~annotated_sex, y = ~DDX3Y,
    split = ~annotated_sex,
    type = "violin",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE)
  ) %>%
  layout(
    yaxis = list(type = "log",
                 range = c(-3, 3))
  )

figUSP9Y <- counts_plus %>%
  plot_ly(
    x = ~annotated_sex, y = ~USP9Y,
    split = ~annotated_sex,
    type = "violin",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE)
  ) %>%
  layout(
    yaxis = list(type = "log",
                 range = c(-3, 3))
  )

figUTY <- counts_plus %>%
  plot_ly(
    x = ~annotated_sex, y = ~UTY,
    split = ~annotated_sex,
    type = "violin",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE)
  ) %>%
  layout(
    yaxis = list(type = "log",
                 range = c(-3, 3))
  )

figZFY <- counts_plus %>%
  plot_ly(
    x = ~annotated_sex, y = ~ZFY,
    split = ~annotated_sex,
    type = "violin",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE)
  ) %>%
  layout(
    yaxis = list(type = "log",
                 range = c(-3, 3))
  )

figXIST <- counts_plus %>%
  plot_ly(
    x = ~annotated_sex, y = ~XIST,
    split = ~annotated_sex,
    type = "violin",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE)
  ) %>%
  layout(
    yaxis = list(type = "log",
                 range = c(-3, 3))
  )

fig <- subplot(figDDX3Y, figUSP9Y, figUTY, figZFY, figXIST,
               nrows = 3, shareY = TRUE) %>%
  layout(title = "Distribution of TPM counts for genes of interest by sex",
         plot_bgcolor = "#e5ecf6",
         showlegend = FALSE
  )

fig



#####################################################
###---LINE PLOTS OF TPM COUNTS FOR THE ODD COMBOS-###
###---OF LOW XIST & LOW Y OR HIGH XIST & HIGH Y---###
###---TO ASSESS FOR ANY PATTERNS WITHIN SAMPLES---###
#####################################################

#  I've updated the code to incorporate an if-else statement.
# This modification ensures that if there are no samples meeting the criteria of 'yes XIST and yes Y,' 
# instead of encountering an error, the output will be a message indicating the absence of such samples, such as 'No samples found for yes XIST and yes Y.' This enhancement aims to handle cases where there might be a lack of data conforming to the specified conditions and provides a more informative response in such scenarios. - Ilsa 

# Two Groups of Interest:

# 1. no XIST and no Y
nxny_counts <- counts_plus %>% dplyr::filter(status_XIST == "no" & status_Y == "no")

# 2. yes XIST and yes Y
yxyy_counts <- counts_plus %>% dplyr::filter(status_XIST == "yes" & status_Y == "yes")

#Check if nxny_counts is not empty
if (nrow(nxny_counts) == 0) {
  print("No samples found for no XIST and no Y.")
} else {
  # Prepping list of gene names for plot functions below:
  gene_names <- c("XIST", "DDX3Y", "USP9Y", "UTY", "ZFY")
  
  # Reshape to work well with plot functions:
  nxny_plots <- reshape(nxny_counts,
                        varying = gene_names,
                        drop = c("AR", "AMELY", "EIF1AY", "KDM5D",
                                 "NLGN4Y", "PRKKY", "TMSB4Y", "SRY",
                                 "status_XIST", "status_Y",
                                 "survival", "followed"),
                        v.names = "TPM_counts",
                        timevar = "gene_names",
                        times = gene_names,
                        ids = row.names(nxny_counts),
                        direction = "long")
  
  
  # Proceed with yxyy_plots if nxny_counts is not empty
  if (nrow(yxyy_counts) == 0) {
    print("No samples found for yes XIST and yes Y.")
  } else {
    yxyy_plots <- reshape(yxyy_counts,
                          varying = gene_names,
                          drop = c("AR", "AMELY", "EIF1AY", "KDM5D",
                                   "NLGN4Y", "PRKKY", "TMSB4Y", "SRY",
                                   "status_XIST", "status_Y",
                                   "survival", "followed"),
                          v.names = "TPM_counts",
                          timevar = "gene_names",
                          times = gene_names,
                          ids = row.names(yxyy_counts),
                          direction = "long")
    
  }
}
# }




# Plot the lines

# Linear y-axis for nxny_plots
# if (exists("nxny_plots")) {

lineplot1_fname <- paste("TPM", study, "LowXIST-LowY_Color-by-Sex.png", sep = "_")
png(lineplot1_fname)
ggplot(data = nxny_plots, aes(x = gene_names, y = TPM_counts,
                              group = id, color = annotated_sex)) +
  scale_color_discrete(guide = "none") +
  geom_point() + geom_line() +
  ggtitle(study,"Samples with No/Low XIST No/Low Y chr, inclusive of repeats") +# Included "study" in the title for convenience.
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
# } else {
#  cat("No data available for 'yes XIST and yes Y' (nxny_plots) to plot.\n")#Output a message indicating no data for nxny_plots
# }
#Check if yxyy_plots exists before plotting, if you get an error saying yxyy_plots does not exist than run if else statement
# if (exists("yxyy_plots")) {

#Linear y-axis for yxyy_plots

lineplot2_fname <- paste("TPM", study, "HighXIST-HighY_Color-by-Sex.png", sep = "_")
png(lineplot2_fname)
ggplot(data = yxyy_plots, aes(x = gene_names, y = TPM_counts,
                              group = id, color = annotated_sex)) +
  scale_color_discrete(guide = "none") +
  geom_point() + geom_line() +
  ggtitle("Samples with XIST & Y markers TPM > 1.0, inclusive of repeats") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Log y-axis for yxyy_plots
loglineplot2_fname <- paste("TPM", study, "HighXIST-HighY_Color-by-Sex_Log-axis.png", sep = "_")
png(loglineplot2_fname)
ggplot(data = yxyy_plots, aes(x = gene_names, y = TPM_counts,
                              group = id, color = annotated_sex)) +
  scale_color_discrete(guide = "none") +
  geom_point() + geom_line() +
  scale_y_continuous(trans = "log10") +
  ggtitle("Samples with XIST & Y markers TPM > 1.0, inclusive of repeats") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
# } else {
#  cat("No data available for 'yes XIST and yes Y' (yxyy_plots) to plot.\n")#Output a message indicating no data for yxyy_plots
# }






















