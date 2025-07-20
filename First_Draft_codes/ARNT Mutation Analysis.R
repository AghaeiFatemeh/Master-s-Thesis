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
##########################   ARNT Analysis    ###########################
#########################################################################
#########################################################################

# Create percentage summary table
dfBar <- dfSampleInfo[, .N, by = .(ARNT, Stage)]
dfBar[, Percent := N / sum(N) * 100, by = ARNT]

# Plot stacked bar chart with percentages
# ARNT vs Tumor Stage
p <- ggplot(dfBar, aes(x = ARNT, y = Percent, fill = Stage)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percent), "%")), 
            position = position_stack(vjust = 0.5), size = 4, color = "white") +
  scale_fill_manual(values = c("pTa" = "#9ecae1", "pT1" = "#fcae91")) +
  labs(title = "Distribution of Tumor Stage by ARNT Mutation Status",
       x = "ARNT Status",
       y = "Percentage",
       fill = "Stage") +
  theme_minimal()

ggsave(filename = paste0(Address.OutputGraphs, "StackedBar_ARNT_vs_Stage.png"),
       plot = p, width = 6, height = 5)

#ARNT vs RNA_class
dfBar_rna <- dfSampleInfo[, .N, by = .(ARNT, RNA_class)]
dfBar_rna[, Percent := N / sum(N) * 100, by = ARNT]

p_rna <- ggplot(dfBar_rna, aes(x = ARNT, y = Percent, fill = RNA_class)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percent), "%")), 
            position = position_stack(vjust = 0.5), size = 4, color = "white") +
  labs(title = "Distribution of RNA Class by ARNT Mutation Status",
       x = "ARNT Status", y = "Percentage", fill = "RNA Class") +
  theme_minimal()

ggsave(paste0(Address.OutputGraphs, "StackedBar_ARNT_vs_RNAclass.png"),
       plot = p_rna, width = 6, height = 5)


#ARNT vs Grade
dfBar_grade <- dfSampleInfo[, .N, by = .(ARNT, Grade)]
dfBar_grade[, Percent := N / sum(N) * 100, by = ARNT]

p_grade <- ggplot(dfBar_grade, aes(x = ARNT, y = Percent, fill = Grade)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percent), "%")), 
            position = position_stack(vjust = 0.5), size = 4, color = "white") +
  labs(title = "Distribution of Tumor Grade by ARNT Mutation Status",
       x = "ARNT Status", y = "Percentage", fill = "Grade") +
  theme_minimal()

ggsave(paste0(Address.OutputGraphs, "StackedBar_ARNT_vs_Grade.png"),
       plot = p_grade, width = 6, height = 5)


#Smoking Status vs ARNT
# Exclude rows with missing Smoking_Status
dfSmokingClean <- dfSampleInfo[!is.na(Smoking_Status)]

# Create summary table
dfBar_smoking <- dfSmokingClean[, .N, by = .(Smoking_Status, ARNT)]
dfBar_smoking[, Percent := N / sum(N) * 100, by = Smoking_Status]

# Plot stacked bar
p_smoking <- ggplot(dfBar_smoking, aes(x = Smoking_Status, y = Percent, fill = ARNT)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percent), "%")), 
            position = position_stack(vjust = 0.5), size = 4, color = "white") +
  labs(title = "ARNT Mutation Distribution by Smoking Status",
       x = "Smoking Status", y = "Percentage", fill = "ARNT") +
  theme_minimal()

# Save the plot
ggsave(paste0(Address.OutputGraphs, "StackedBar_Smoking_vs_ARNT_clean.png"),
       plot = p_smoking, width = 6, height = 5)



# Boxplot ARNT vs Age
p_age <- ggplot(dfSampleInfo, aes(x = ARNT, y = Age, fill = ARNT)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, color = "black", size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Age Distribution by ARNT Mutation Status",
       x = "ARNT Status", y = "Age") +
  scale_fill_manual(values = c("WT" = "#9ecae1", "MUT" = "#fcae91")) +
  theme(plot.title = element_text(hjust = 0.5))

# Save the plot
ggsave(paste0(Address.OutputGraphs, "Boxplot_ARNT_vs_Age.png"),
       plot = p_age, width = 6, height = 5)


test_age <- wilcox.test(Age ~ ARNT, data = dfSampleInfo)
pval_age <- format.pval(test_age$p.value, digits = 3, eps = .001)

# Add p-value as subtitle
p_age <- p_age + labs(subtitle = paste("Wilcoxon p =", pval_age))

# Save updated plot
ggsave(paste0(Address.OutputGraphs, "Boxplot_ARNT_vs_Age_withPval.png"),
       plot = p_age, width = 6, height = 5)

#########################################################################
#########################################################################
#######################   ARNT Analysis-Test    #########################
#########################################################################
#########################################################################

#ARNT vs Stage
table_stage <- table(dfSampleInfo$ARNT, dfSampleInfo$Stage)

test_stage <- fisher.test(table_stage)

result_stage <- data.table(
  Comparison = "ARNT vs Stage",
  Test = "Fisher's Exact Test",
  P_value = round(test_stage$p.value, 5),
  Odds_Ratio = ifelse(length(test_stage$estimate) == 1, round(test_stage$estimate, 3), NA),
  CI_95 = ifelse(length(test_stage$conf.int) == 2,
                 paste(round(test_stage$conf.int, 3), collapse = " - "),
                 NA)
)


# Save as TXT (human-readable explanation)
writeLines(c(
  "Statistical Test: Fisher's Exact Test",
  "Comparison: ARNT mutation status vs Tumor Stage (pTa vs pT1)",
  paste("P-value:", round(test_stage$p.value, 5)),
  paste("Odds Ratio:", round(test_stage$estimate, 3)),
  paste("95% CI:", paste(round(test_stage$conf.int, 3), collapse = " - "))
), con = paste0(Address.OutputTables, "FisherTest_ARNT_vs_Stage.txt"))

# Save as LaTeX scientific table
latex_table <- paste0(
  "\\begin{table}[h]\n",
  "\\centering\n",
  "\\begin{tabular}{lccc}\n",
  "\\hline\n",
  "Comparison & Test & P-value & 95\\% CI (Odds Ratio) \\\\\n",
  "\\hline\n",
  "ARNT vs Stage & Fisher's Exact Test & ", round(test_stage$p.value, 5), 
  " & ", paste0(round(test_stage$conf.int, 3), collapse = " -- "), 
  " (", round(test_stage$estimate, 3), ") \\\\\n",
  "\\hline\n",
  "\\end{tabular}\n",
  "\\caption{Association between ARNT mutation status and tumor stage (pTa vs pT1).}\n",
  "\\label{tab:arnt_stage_fisher}\n",
  "\\end{table}"
)

writeLines(latex_table, con = paste0(Address.OutputTables, "FisherTest_ARNT_vs_Stage.tex"))


# ARNT vs RNA Class
# Create contingency table
table_rna <- table(dfSampleInfo$ARNT, dfSampleInfo$RNA_class)

# Run Fisher's Exact Test
test_rna <- fisher.test(table_rna)

# Create summary table
result_rna <- data.table(
  Comparison = "ARNT vs RNA_class",
  Test = "Fisher's Exact Test",
  P_value = round(test_rna$p.value, 5),
  Odds_Ratio = ifelse(length(test_rna$estimate) == 1, round(test_rna$estimate, 3), NA),
  CI_95 = ifelse(length(test_rna$conf.int) == 2,
                 paste(round(test_rna$conf.int, 3), collapse = " - "),
                 NA)
)

# Save as TXT (human-readable)
writeLines(c(
  "Statistical Test: Fisher's Exact Test",
  "Comparison: ARNT mutation status vs RNA molecular class (A, B, C)",
  paste("P-value:", round(test_rna$p.value, 5)),
  paste("Odds Ratio:", ifelse(is.na(result_rna$Odds_Ratio), "Not Applicable", result_rna$Odds_Ratio)),
  paste("95% CI:", ifelse(is.na(result_rna$CI_95), "Not Applicable", result_rna$CI_95))
), con = paste0(Address.OutputTables, "FisherTest_ARNT_vs_RNAclass.txt"))

# Save as LaTeX
latex_table <- paste0(
  "\\begin{table}[h]\n",
  "\\centering\n",
  "\\begin{tabular}{lccc}\n",
  "\\hline\n",
  "Comparison & Test & P-value & 95\\% CI (Odds Ratio) \\\\\n",
  "\\hline\n",
  "ARNT vs RNA class & Fisher's Exact Test & ", round(test_rna$p.value, 5), " & ",
  ifelse(length(test_rna$conf.int) == 2,
         paste0(paste(round(test_rna$conf.int, 3), collapse = " -- "), 
                " (", round(test_rna$estimate, 3), ")"),
         "Not Applicable"),
  " \\\\\n",
  "\\hline\n",
  "\\end{tabular}\n",
  "\\caption{Association between ARNT mutation status and RNA molecular class.}\n",
  "\\label{tab:arnt_rna_fisher}\n",
  "\\end{table}"
)

writeLines(latex_table, con = paste0(Address.OutputTables, "FisherTest_ARNT_vs_RNAclass.tex"))


# ARNT vs Grade
# Create contingency table
table_grade <- table(dfSampleInfo$ARNT, dfSampleInfo$Grade)

# Run Fisher's Exact Test
test_grade <- fisher.test(table_grade)

# Create result summary table
result_grade <- data.table(
  Comparison = "ARNT vs Grade",
  Test = "Fisher's Exact Test",
  P_value = round(test_grade$p.value, 5),
  Odds_Ratio = ifelse(length(test_grade$estimate) == 1, round(test_grade$estimate, 3), NA),
  CI_95 = ifelse(length(test_grade$conf.int) == 2,
                 paste(round(test_grade$conf.int, 3), collapse = " - "),
                 NA)
)

# Save as TXT (human-readable)
writeLines(c(
  "Statistical Test: Fisher's Exact Test",
  "Comparison: ARNT mutation status vs Tumor Grade (G1, G2, G3)",
  paste("P-value:", result_grade$P_value),
  paste("Odds Ratio:", ifelse(is.na(result_grade$Odds_Ratio), "Not Applicable", result_grade$Odds_Ratio)),
  paste("95% CI:", ifelse(is.na(result_grade$CI_95), "Not Applicable", result_grade$CI_95))
), con = paste0(Address.OutputTables, "FisherTest_ARNT_vs_Grade.txt"))

# Save as LaTeX table
latex_table <- paste0(
  "\\begin{table}[h]\n",
  "\\centering\n",
  "\\begin{tabular}{lccc}\n",
  "\\hline\n",
  "Comparison & Test & P-value & 95\\% CI (Odds Ratio) \\\\\n",
  "\\hline\n",
  "ARNT vs Grade & Fisher's Exact Test & ", result_grade$P_value, " & ",
  ifelse(!is.na(result_grade$CI_95),
         paste0(result_grade$CI_95, " (", result_grade$Odds_Ratio, ")"),
         "Not Applicable"),
  " \\\\\n",
  "\\hline\n",
  "\\end{tabular}\n",
  "\\caption{Association between ARNT mutation status and tumor grade (G1, G2, G3).}\n",
  "\\label{tab:arnt_grade_fisher}\n",
  "\\end{table}"
)

writeLines(latex_table, con = paste0(Address.OutputTables, "FisherTest_ARNT_vs_Grade.tex"))



# ARNT vs Smoking Status

# Remove NA rows in Smoking_Status
dfSmokingClean <- dfSampleInfo[!is.na(Smoking_Status)]

# Contingency table
table_smoking <- table(dfSmokingClean$Smoking_Status, dfSmokingClean$ARNT)

# Run Fisher's Exact Test
test_smoking <- fisher.test(table_smoking)

# Create result summary
result_smoking <- data.table(
  Comparison = "Smoking_Status vs ARNT",
  Test = "Fisher's Exact Test",
  P_value = round(test_smoking$p.value, 5),
  Odds_Ratio = ifelse(length(test_smoking$estimate) == 1, round(test_smoking$estimate, 3), NA),
  CI_95 = ifelse(length(test_smoking$conf.int) == 2,
                 paste(round(test_smoking$conf.int, 3), collapse = " - "),
                 NA)
)

# Save as TXT (human-readable)
writeLines(c(
  "Statistical Test: Fisher's Exact Test",
  "Comparison: Smoking status (Smoker, Non-smoker, Ex-smoker) vs ARNT mutation status",
  paste("P-value:", result_smoking$P_value),
  paste("Odds Ratio:", ifelse(is.na(result_smoking$Odds_Ratio), "Not Applicable", result_smoking$Odds_Ratio)),
  paste("95% CI:", ifelse(is.na(result_smoking$CI_95), "Not Applicable", result_smoking$CI_95))
), con = paste0(Address.OutputTables, "FisherTest_Smoking_vs_ARNT.txt"))

# Save as LaTeX table
latex_table <- paste0(
  "\\begin{table}[h]\n",
  "\\centering\n",
  "\\begin{tabular}{lccc}\n",
  "\\hline\n",
  "Comparison & Test & P-value & 95\\% CI (Odds Ratio) \\\\\n",
  "\\hline\n",
  "Smoking Status vs ARNT & Fisher's Exact Test & ", result_smoking$P_value, " & ",
  ifelse(!is.na(result_smoking$CI_95),
         paste0(result_smoking$CI_95, " (", result_smoking$Odds_Ratio, ")"),
         "Not Applicable"),
  " \\\\\n",
  "\\hline\n",
  "\\end{tabular}\n",
  "\\caption{Association between smoking status and ARNT mutation status.}\n",
  "\\label{tab:smoking_arnt_fisher}\n",
  "\\end{table}"
)

writeLines(latex_table, con = paste0(Address.OutputTables, "FisherTest_Smoking_vs_ARNT.tex"))


# ARNT vs Age

# Run Wilcoxon test
test_age <- wilcox.test(Age ~ ARNT, data = dfSampleInfo)

# Create summary result
result_age <- data.table(
  Comparison = "ARNT vs Age",
  Test = "Wilcoxon Rank-Sum Test",
  P_value = round(test_age$p.value, 5),
  W_statistic = round(test_age$statistic, 2)
)

# Save as TXT (human-readable)
writeLines(c(
  "Statistical Test: Wilcoxon Rank-Sum Test (Mann–Whitney U Test)",
  "Comparison: Age distribution between ARNT-MUT and ARNT-WT groups",
  paste("P-value:", result_age$P_value),
  paste("W statistic:", result_age$W_statistic)
), con = paste0(Address.OutputTables, "WilcoxonTest_ARNT_vs_Age.txt"))

# Save as LaTeX table
latex_table <- paste0(
  "\\begin{table}[h]\n",
  "\\centering\n",
  "\\begin{tabular}{lccc}\n",
  "\\hline\n",
  "Comparison & Test & P-value & W statistic \\\\\n",
  "\\hline\n",
  "ARNT vs Age & Wilcoxon Rank-Sum Test & ", result_age$P_value, " & ", result_age$W_statistic, " \\\\\n",
  "\\hline\n",
  "\\end{tabular}\n",
  "\\caption{Comparison of patient age between ARNT mutation groups using the Wilcoxon Rank-Sum Test.}\n",
  "\\label{tab:arnt_age_wilcoxon}\n",
  "\\end{table}"
)

writeLines(latex_table, con = paste0(Address.OutputTables, "WilcoxonTest_ARNT_vs_Age.tex"))


table(dfSampleInfo$Grade)
prop.table(table(dfSampleInfo$Grade)) * 100

table(dfSampleInfo$Smoking_Status)
prop.table(table(dfSampleInfo$Smoking_Status)) * 100
