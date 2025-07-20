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





# Boxplot ARNT Expression vs ARNT Mutation Status
p_expr <- ggplot(dfSampleInfo, aes(x = ARNT, y = ARNT_expression, fill = ARNT)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, color = "black", size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "ARNT Expression by Mutation Status",
       x = "ARNT Status", y = "ARNT Expression") +
  scale_fill_manual(values = c("WT" = "#9ecae1", "MUT" = "#fcae91")) +
  theme(plot.title = element_text(hjust = 0.5))

# Save initial plot
ggsave(paste0(Address.OutputGraphs, "Boxplot_ARNTexpression_vs_Status.png"),
       plot = p_expr, width = 6, height = 5)

# Wilcoxon test
test_expr <- wilcox.test(ARNT_expression ~ ARNT, data = dfSampleInfo)
pval_expr <- format.pval(test_expr$p.value, digits = 3, eps = .001)

# Add p-value to plot
p_expr <- p_expr + labs(subtitle = paste("Wilcoxon p =", pval_expr))

# Save updated plot with p-value
ggsave(paste0(Address.OutputGraphs, "Boxplot_ARNTexpression_vs_Status_withPval.png"),
       plot = p_expr, width = 6, height = 5)


# Create result summary
result_expr <- data.table(
  Comparison = "ARNT Expression vs ARNT Status",
  Test = "Wilcoxon Rank-Sum Test",
  P_value = round(test_expr$p.value, 5),
  W_statistic = round(test_expr$statistic, 2)
)

# Save as TXT
writeLines(c(
  "Statistical Test: Wilcoxon Rank-Sum Test",
  "Comparison: ARNT expression levels between ARNT-MUT and ARNT-WT groups",
  paste("P-value:", result_expr$P_value),
  paste("W statistic:", result_expr$W_statistic)
), con = paste0(Address.OutputTables, "Wilcoxon_ARNTexpression_vs_Status.txt"))

# Save as LaTeX table
latex_expr <- paste0(
  "\\begin{table}[h]\n",
  "\\centering\n",
  "\\begin{tabular}{lccc}\n",
  "\\hline\n",
  "Comparison & Test & P-value & W statistic \\\\\n",
  "\\hline\n",
  "ARNT Expression vs ARNT Status & Wilcoxon Rank-Sum Test & ",
  result_expr$P_value, " & ", result_expr$W_statistic, " \\\\\n",
  "\\hline\n",
  "\\end{tabular}\n",
  "\\caption{Comparison of ARNT expression between mutation status groups using Wilcoxon test.}\n",
  "\\label{tab:arnt_expr_status}\n",
  "\\end{table}"
)

writeLines(latex_expr, con = paste0(Address.OutputTables, "Wilcoxon_ARNTexpression_vs_Status.tex"))





# Step 1: Define 75th percentile cutoff
cutoff_75 <- quantile(dfSampleInfo$ARNT_expression, probs = 0.75, na.rm = TRUE)

# Step 2: Classify into "High" and "Low" expression groups
dfSampleInfo$ARNT_expr_group <- ifelse(dfSampleInfo$ARNT_expression > cutoff_75, "High", "Low")

# Step 3: Convert to factor (for consistency in plots and summaries)
dfSampleInfo$ARNT_expr_group <- factor(dfSampleInfo$ARNT_expr_group, levels = c("Low", "High"))

# Step 4: Create contingency table: ARNT expression group vs Tumor Stage
table_expr_stage <- table(dfSampleInfo$ARNT_expr_group, dfSampleInfo$Stage)

# Step 5: Fisher's Exact Test
test_expr_stage <- fisher.test(table_expr_stage)

# Step 6: Save result as human-readable text
writeLines(c(
  "Statistical Test: Fisher's Exact Test",
  "Comparison: ARNT expression group (High/Low) vs Tumor Stage (pTa/pT1)",
  paste("P-value:", round(test_expr_stage$p.value, 5)),
  paste("Odds Ratio:", round(test_expr_stage$estimate, 3)),
  paste("95% CI:", paste(round(test_expr_stage$conf.int, 3), collapse = " - "))
), con = paste0(Address.OutputTables, "FisherTest_ARNTExprGroup_vs_Stage.txt"))

# Step 7: Save result as LaTeX table
latex_expr_stage <- paste0(
  "\\begin{table}[h]\n",
  "\\centering\n",
  "\\begin{tabular}{lccc}\n",
  "\\hline\n",
  "Comparison & Test & P-value & 95\\% CI (Odds Ratio) \\\\\n",
  "\\hline\n",
  "ARNT Expression Group vs Stage & Fisher's Exact Test & ", 
  round(test_expr_stage$p.value, 5), " & ", 
  paste0(round(test_expr_stage$conf.int, 3), collapse = " -- "), 
  " (", round(test_expr_stage$estimate, 3), ") \\\\\n",
  "\\hline\n",
  "\\end{tabular}\n",
  "\\caption{Association between ARNT expression group (High vs Low) and tumor stage (pTa vs pT1).}\n",
  "\\label{tab:arnt_expr_stage_fisher}\n",
  "\\end{table}"
)

writeLines(latex_expr_stage, con = paste0(Address.OutputTables, "FisherTest_ARNTExprGroup_vs_Stage.tex"))





# Step 1: Prepare summary data for the bar plot
dfBar_expr <- dfSampleInfo[, .N, by = .(ARNT_expr_group, Stage)]
dfBar_expr[, Total := sum(N), by = ARNT_expr_group]
dfBar_expr[, Percent := 100 * N / Total]

# Step 2: Plot
p_expr_stage <- ggplot(dfBar_expr, aes(x = ARNT_expr_group, y = Percent, fill = Stage)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percent), "%")),
            position = position_stack(vjust = 0.5), size = 4, color = "white") +
  scale_fill_manual(values = c("pTa" = "#9ecae1", "pT1" = "#fcae91")) +
  labs(title = "Distribution of Tumor Stage by ARNT Expression Group",
       x = "ARNT Expression Group",
       y = "Percentage",
       fill = "Tumor Stage") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 3: Save plot
ggsave(filename = paste0(Address.OutputGraphs, "StackedBar_ARNTExprGroup_vs_Stage.png"),
       plot = p_expr_stage, width = 6, height = 5)


