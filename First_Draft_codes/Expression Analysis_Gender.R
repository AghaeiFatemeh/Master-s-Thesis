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
##################   Expression DAta Preparation    #####################
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


# Step 1: Define 75th percentile cutoff
cutoff_75 <- quantile(dfSampleInfo$ARNT_expression, probs = 0.75, na.rm = TRUE)

# Step 2: Classify into "High" and "Low" expression groups
dfSampleInfo$ARNT_expr_group <- ifelse(dfSampleInfo$ARNT_expression > cutoff_75, "High", "Low")

# Step 3: Convert to factor (for consistency in plots and summaries)
dfSampleInfo$ARNT_expr_group <- factor(dfSampleInfo$ARNT_expr_group, levels = c("Low", "High"))

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

#########################################################################
#########################################################################
#############################   Analysis    #############################
#########################################################################
#########################################################################


# Step 1: Subset to complete cases for both columns
dfTemp <- dfSampleInfo[!is.na(ARNT_expr_group) & !is.na(Gender)]

# Step 2: Create contingency table
table_expr_gender <- table(dfTemp$Gender, dfTemp$ARNT_expr_group)

# Step 3: Fisher's Exact Test
test_expr_gender <- fisher.test(table_expr_gender)

# Step 4: Save LaTeX table
latex_expr_gender <- paste0(
  "\\begin{table}[h]\n",
  "\\centering\n",
  "\\begin{tabular}{lccc}\n",
  "\\hline\n",
  "Comparison & Test & P-value & 95\\% CI (Odds Ratio) \\\\\n",
  "\\hline\n",
  "Gender vs ARNT Expression Group & Fisher's Exact Test & ", 
  round(test_expr_gender$p.value, 5), " & ", 
  paste0(round(test_expr_gender$conf.int, 3), collapse = " -- "), 
  " (", round(test_expr_gender$estimate, 3), ") \\\\\n",
  "\\hline\n",
  "\\end{tabular}\n",
  "\\caption{Association between Gender and ARNT expression group (High vs Low).}\n",
  "\\label{tab:gender_vs_arnt_expr_fisher}\n",
  "\\end{table}"
)
writeLines(latex_expr_gender, con = paste0(Address.OutputTables, "FisherTest_Gender_vs_ARNTExprGroup.tex"))

# Step 5: Prepare data for bar plot
dfBar_gender <- dfTemp[, .N, by = .(Gender, ARNT_expr_group)]
dfBar_gender[, Total := sum(N), by = Gender]
dfBar_gender[, Percent := 100 * N / Total]

# Step 6: Plot with p-value annotated
library(ggplot2)
p_gender_expr <- ggplot(dfBar_gender, aes(x = Gender, y = Percent, fill = ARNT_expr_group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percent), "%")),
            position = position_stack(vjust = 0.5), size = 4, color = "white") +
  annotate("text", x = 1.5, y = 105, 
           label = paste("Fisher p =", round(test_expr_gender$p.value, 4)),
           size = 3, fontface = "italic") +
  scale_fill_manual(values = c("Low" = "#9ecae1", "High" = "#fcae91")) +
  labs(title = "Distribution of ARNT Expression Group by Gender",
       x = "Gender",
       y = "Percentage",
       fill = "ARNT Expression Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.3))

# Step 7: Save plot
ggsave(filename = paste0(Address.OutputGraphs, "StackedBar_Gender_vs_ARNTExprGroup.png"),
       plot = p_gender_expr, width = 6, height = 5)
