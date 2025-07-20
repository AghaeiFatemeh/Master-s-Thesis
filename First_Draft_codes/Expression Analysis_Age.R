#########################################################################
#########################################################################
##############################   Packages    ############################
#########################################################################
#########################################################################

# Load packages
# List required packages
packages <- c("data.table","ggplot2", "readxl", "RColorBrewer", "nnet", "stats")

# Load all packages
lapply(packages, library, readxl, character.only = TRUE)


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


# Read the large CSV file
dfVeracyte <- fread(paste0(Address.DataRaw, "veracyte_nmibc_transcriptomes.csv"))


# Read the Clinical Data
dfClinical <- as.data.table(read_excel(paste0(Address.DataRaw, "Clinical_data.xlsx")))

cleandfclinical <- fread(paste0(Address.Datatemp, "cleandfclinical.csv"))

# Filter rows with non-finite ARNT_expression
cleandfclinical <- cleandfclinical[is.finite(ARNT_expression)]


#########################################################################
#########################################################################
#######################   ARNT Expression vs Age    #####################
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




# 1. Compute tertile cut-points based on Age
age_breaks <- quantile(dfSampleInfo$Age, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# 2. Define the plain labels
tertile_labels <- c("Young", "Middle", "Old")

# 3. Cut into tertiles
dfSampleInfo[
  , age_group := cut(
    Age,
    breaks         = age_breaks,
    labels         = tertile_labels,
    include.lowest = TRUE
  )
]

# 4. Build the “label + (min–max)” text for each tertile
lower      <- round(age_breaks[-length(age_breaks)], 0)
upper      <- round(age_breaks[-1], 0)
range_text <- paste0("(", lower, "–", upper, ")")
full_labels <- paste0(tertile_labels, "\n", range_text)
names(full_labels) <- tertile_labels

# 5. Define a colour palette
age_cols <- setNames(
  brewer.pal(n = length(tertile_labels), name = "Set3"),
  tertile_labels
)

# 6. Draw the boxplot
p_expr <- ggplot(dfSampleInfo, aes(x = age_group, y = ARNT_expression, fill = age_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  scale_x_discrete(labels = full_labels) +
  scale_fill_manual(values = age_cols, labels = full_labels) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.8, color = "black") +
  theme_minimal() +
  labs(
    title = "ARNT Expression by Age Group",
    x     = "Age Group",
    y     = "ARNT Expression"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# 7. Kruskal–Wallis test
kw_res  <- kruskal.test(ARNT_expression ~ age_group, data = dfSampleInfo)
pval_kw <- format.pval(kw_res$p.value, digits = 3, eps = .001)

# 8. Add p-value subtitle and save
p_expr <- p_expr + labs(subtitle = paste("Kruskal–Wallis p =", pval_kw))

ggsave(
  filename = file.path(Address.OutputGraphs, "Boxplot_ARNTexpression_vs_agegroup_kruskal.png"),
  plot     = p_expr,
  width    = 6,
  height   = 5
)