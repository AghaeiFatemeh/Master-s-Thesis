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




# Boxplot ARNT Expression vs UROMOL_NMIBC_Class
p_expr <- ggplot(dfSampleInfo, aes(x = UROMOL_NMIBC_Class, y = ARNT_expression, fill = UROMOL_NMIBC_Class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, color = "black", size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "ARNT Expression by UROMOL_NMIBC_Class",
       x = "UROMOL_NMIBC_Class", y = "ARNT Expression") +
  scale_fill_manual(values = c("Class_1" = "#9ecae1", "Class_2a" = "#fcae91", "Class_2b" = "#a1d99b","Class_3"  = "#fc9272")) +
  theme(plot.title = element_text(hjust = 0.5))

# Kruskal–Wallis test
test_expr <- kruskal.test(ARNT_expression ~ UROMOL_NMIBC_Class, data = dfSampleInfo)
pval_expr <- format.pval(test_expr$p.value, digits = 3, eps = .001)

# Add p-value to plot
p_expr <- p_expr + labs(subtitle = paste("Kruskal–Wallis p =", pval_expr))

# Save updated plot with p-value
ggsave(paste0(Address.OutputGraphs, "Boxplot_ARNTexpression_vs_UROMOL_NMIBC_Class_withPval.png"),
       plot = p_expr, width = 6, height = 5)





# Boxplot ARNT Expression vs RNA_class
p_expr <- ggplot(dfSampleInfo, aes(x = RNA_class, y = ARNT_expression, fill = RNA_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, color = "black", size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "ARNT Expression by RNA_class",
       x = "RNA_class", y = "ARNT Expression") +
  scale_fill_manual(values = c("A" = "#9ecae1", "B" = "#fcae91", "C" = "#a1d99b")) +
  theme(plot.title = element_text(hjust = 0.5))

# Kruskal–Wallis test
test_expr <- kruskal.test(ARNT_expression ~ RNA_class, data = dfSampleInfo)
pval_expr <- format.pval(test_expr$p.value, digits = 3, eps = .001)

# Add p-value to plot
p_expr <- p_expr + labs(subtitle = paste("Kruskal–Wallis p =", pval_expr))

# Save updated plot with p-value
ggsave(paste0(Address.OutputGraphs, "Boxplot_ARNTexpression_vs_RNA_class_withPval.png"),
       plot = p_expr, width = 6, height = 5)
