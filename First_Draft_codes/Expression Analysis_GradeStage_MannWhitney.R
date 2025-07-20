#########################################################################
#########################################################################
##############################   Packages    ############################
#########################################################################
#########################################################################

# Load packages
# List required packages
packages <- c("data.table", "ggplot2")

# Load all packages
lapply(packages, library, character.only = TRUE)


# Clear all objects from the environment to avoid conflicts
rm(list = ls())
# Trigger garbage collection to free up memory
gc()

#########################################################################
#########################################################################
##############################   Variables    ###########################
#########################################################################
#########################################################################

# Define the file path for raw data
Address.DataRaw <- "/rds/projects/a/arnoldry-arnt/project_Fatemeh_ARNT/"
# Define the file path for temporary data storage
Address.Datatemp <- "/rds/projects/a/arnoldry-arnt/project_Fatemeh_ARNT/Data/temp/"
# Define the file path for cleaned data
Address.DataClean <- "/rds/projects/a/arnoldry-arnt/project_Fatemeh_ARNT/Data/Clean/"
# Define the output path for saving graphs
Address.OutputGraphs <- "/rds/projects/a/arnoldry-arnt/project_Fatemeh_ARNT/Output/Graphs/"
# Define the output path for saving tables
Address.OutputTables <- "/rds/projects/a/arnoldry-arnt/project_Fatemeh_ARNT/Output/Tables/"

#########################################################################
#########################################################################
###############################   Codes    ##############################
#########################################################################
#########################################################################

# Load sample information from the raw data directory
# This file contains metadata such as sample IDs, conditions, or patient information
dfSampleInfo <- fread(paste0(Address.DataRaw, "sample_info.txt"))

# Load normalized RNA-seq count data from the raw data directory
# This file contains gene expression values for samples in the study
dfRNACounts <- fread(paste0(Address.DataRaw, "BCPP_NMIBC_normalised_RNA-seq_counts.txt"))

#########################################################################
#########################################################################
######################   Expression Analysis    #########################
#########################################################################
#########################################################################


# Step 1: Extract ARNT expression row
arnt_expr_row <- dfRNACounts[V1 == "ENSG00000143437"]

# Step 2: Transpose the row (drop the gene ID column)
arnt_expr <- transpose(arnt_expr_row[, -"V1"])

# Step 3: Create a data.table of ARNT expression values with sample names
arnt_expr_dt <- data.table(
  Sample = names(dfRNACounts)[-1],  # all column names except V1
  ARNT_expression = as.numeric(unlist(arnt_expr))
)

# Step 4: Merge with your metadata
dfSampleInfo <- merge(dfSampleInfo, arnt_expr_dt, by = "Sample", all.x = TRUE)



# Grade and Stage Combination
dfSampleInfo[, Grade_Stage := paste0(gsub("Grade", "G", Grade), Stage)]



# Step 1: Create the new labeled variable for Grade_Stage_Class
dfSampleInfo$Grade_Stage_Class <- NA

dfSampleInfo$Grade_Stage_Class[dfSampleInfo$Grade_Stage == "G1pTa"] <- "Low"
dfSampleInfo$Grade_Stage_Class[dfSampleInfo$Grade_Stage == "G3pT1"] <- "High"

# Step 2: Convert to factor
dfSampleInfo$Grade_Stage_Class <- factor(dfSampleInfo$Grade_Stage_Class, levels = c("Low", "High"))

# Step 3: Filter only the relevant rows (G1pTa and G3pT1)
dfFiltered <- dfSampleInfo[!is.na(Grade_Stage_Class), ]

# 0. Set factor levels and labels for alignment
dfFiltered$Grade_Stage_Class <- factor(
  dfFiltered$Grade_Stage_Class,
  levels = c("Low", "High"),
  labels = c("Low (G1pTa)", "High (G3pT1)")
)

# 1. Boxplot: ARNT Expression vs Grade Stage Class
p_expr <- ggplot(dfFiltered, 
                 aes(x = Grade_Stage_Class, 
                     y = ARNT_expression, 
                     fill = Grade_Stage_Class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.8, color = "black") +
  theme_minimal() +
  labs(
    title = "ARNT Expression by Grade Stage Class",
    subtitle = paste("Wilcoxon p =", format.pval(
      wilcox.test(ARNT_expression ~ Grade_Stage_Class, data = dfFiltered, exact = FALSE)$p.value,
      digits = 3, eps = .001)),
    x     = "Grade Stage Class",
    y     = "ARNT Expression",
    fill  = "Class"
  ) +
  scale_fill_manual(
    values = c("Low (G1pTa)" = "#9ecae1", "High (G3pT1)" = "#fcae91")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )

# 2. Save the plot
ggsave(
  filename = file.path(Address.OutputGraphs,
                       "Boxplot_ARNTexpression_vs_GradeStageClass_withPval.png"),
  plot   = p_expr,
  width  = 6,
  height = 5
)
