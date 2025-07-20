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




# Step 1: Contingency table
table_expr_gradestage <- table(dfSampleInfo$ARNT_expr_group, dfSampleInfo$Grade_Stage)

# Step 2: Fisher’s Exact Test
test_expr_gradestage <- fisher.test(table_expr_gradestage)

# Step 3: Save result as human-readable text
writeLines(c(
  "Statistical Test: Fisher's Exact Test",
  "Comparison: ARNT expression group (High/Low) vs Grade_Stage",
  paste("P-value:", round(test_expr_gradestage$p.value, 5)),
  if (!is.null(test_expr_gradestage$estimate))
    paste("Odds Ratio:", round(test_expr_gradestage$estimate, 3)) else "",
  if (!is.null(test_expr_gradestage$conf.int))
    paste("95% CI:", paste(round(test_expr_gradestage$conf.int, 3), collapse = " - ")) else ""
), con = paste0(Address.OutputTables, "FisherTest_ARNTExpr_vs_GradeStage.txt"))

# Step 4: Save as LaTeX
latex_table_expr_gradestage <- paste0(
  "\\begin{table}[h]\n",
  "\\centering\n",
  "\\begin{tabular}{lccc}\n",
  "\\hline\n",
  "Comparison & Test & P-value & 95\\% CI (Odds Ratio) \\\\\n",
  "\\hline\n",
  "ARNT Expression Group vs Grade+Stage & Fisher's Exact Test & ",
  round(test_expr_gradestage$p.value, 5), " & ",
  if (!is.null(test_expr_gradestage$conf.int))
    paste0(round(test_expr_gradestage$conf.int, 3), collapse = " -- ") else "NA",
  if (!is.null(test_expr_gradestage$estimate))
    paste0(" (", round(test_expr_gradestage$estimate, 3), ")") else "",
  " \\\\\n\\hline\n",
  "\\end{tabular}\n",
  "\\caption{Association between ARNT expression group and combined grade + stage.}\n",
  "\\label{tab:arnt_expr_gradestage_fisher}\n",
  "\\end{table}"
)

writeLines(latex_table_expr_gradestage, con = paste0(Address.OutputTables, "FisherTest_ARNTExpr_vs_GradeStage.tex"))

# Step 5: Prepare data for stacked bar chart
dfBar_expr_grade <- dfSampleInfo[, .N, by = .(ARNT_expr_group, Grade_Stage)]
dfBar_expr_grade[, Total := sum(N), by = ARNT_expr_group]
dfBar_expr_grade[, Percent := 100 * N / Total]
# Step 6 (updated): Plot with p-value as subtitle
p_expr_grade <- ggplot(dfBar_expr_grade, aes(x = ARNT_expr_group, y = Percent, fill = Grade_Stage)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percent), "%")),
            position = position_stack(vjust = 0.5), size = 4, color = "white") +
  labs(
    title = "Distribution of Grade+Stage by ARNT Expression Group",
    subtitle = paste("Fisher's Exact Test p =", format.pval(test_expr_gradestage$p.value, digits = 3, eps = 0.001)),
    x = "ARNT Expression Group",
    y = "Percentage",
    fill = "Grade + Stage"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

# Step 7: Save plot with subtitle
ggsave(filename = paste0(Address.OutputGraphs, "StackedBar_ARNTExpr_vs_GradeStage.png"),
       plot = p_expr_grade, width = 6, height = 5)



# Step 1: Create the new labeled variable for Grade_Stage_Class
dfSampleInfo$Grade_Stage_Class <- NA

dfSampleInfo$Grade_Stage_Class[dfSampleInfo$Grade_Stage == "G1pTa"] <- "Low (G1pTa)"
dfSampleInfo$Grade_Stage_Class[dfSampleInfo$Grade_Stage == "G3pT1"] <- "High (G3pT1)"

# Step 2: Convert to factor (fix this line!)
dfSampleInfo$Grade_Stage_Class <- factor(dfSampleInfo$Grade_Stage_Class,
                                         levels = c("Low (G1pTa)", "High (G3pT1)"))

# Step 3: Filter only the relevant rows (G1pTa and G3pT1)
dfFiltered <- dfSampleInfo[!is.na(Grade_Stage_Class), ]

# Step 4: Contingency table and Fisher test
table_expr_grade_stage_class <- table(dfFiltered$ARNT_expr_group, dfFiltered$Grade_Stage_Class)
test_expr_grade_stage_class <- fisher.test(table_expr_grade_stage_class)

# Step 5: Save as LaTeX table
latex_expr_grade_stage_class <- paste0(
  "\\begin{table}[h]\n",
  "\\centering\n",
  "\\begin{tabular}{lccc}\n",
  "\\hline\n",
  "Comparison & Test & P-value & 95\\% CI (Odds Ratio) \\\\\n",
  "\\hline\n",
  "ARNT Expression Group vs Grade-Stage Class & Fisher's Exact Test & ", 
  round(test_expr_grade_stage_class$p.value, 5), " & ", 
  paste0(round(test_expr_grade_stage_class$conf.int, 3), collapse = " -- "), 
  " (", round(test_expr_grade_stage_class$estimate, 3), ") \\\\\n",
  "\\hline\n",
  "\\end{tabular}\n",
  "\\caption{Association between ARNT expression group (High vs Low) and Grade-Stage class (G1pTa vs G3pT1).}\n",
  "\\label{tab:arnt_expr_vs_grade_stage_class}\n",
  "\\end{table}"
)

writeLines(latex_expr_grade_stage_class, con = paste0(Address.OutputTables, "FisherTest_ARNTExpr_vs_GradeStageClass.tex"))
# Step 6: Prepare data for stacked bar plot
dfBar_class <- dfFiltered[, .N, by = .(ARNT_expr_group, Grade_Stage_Class)]
dfBar_class[, Grade_Stage_Class := factor(Grade_Stage_Class,
                                          levels = c("Low (G1pTa)", "High (G3pT1)"))]
dfBar_class[, Total := sum(N), by = ARNT_expr_group]
dfBar_class[, Percent := 100 * N / Total]

# Step 7: Plot with p-value on the plot
p_class <- ggplot(dfBar_class, aes(x = ARNT_expr_group, y = Percent, fill = Grade_Stage_Class)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percent), "%")),
            position = position_stack(vjust = 0.5), size = 4, color = "white") +
  scale_fill_manual(values = c("Low (G1pTa)" = "#9ecae1", "High (G3pT1)" = "#fcae91")) +
  labs(title = "Grade-Stage Class Distribution by ARNT Expression Group",
       x = "ARNT Expression Group",
       y = "Percentage",
       fill = "Grade-Stage Class") +
  annotate("text", x = 1.5, y = 105, 
           label = paste0("Fisher p = ", round(test_expr_grade_stage_class$p.value, 4)),
           size = 4, fontface = "italic") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.3))

# Step 8: Save plot
ggsave(filename = paste0(Address.OutputGraphs, "StackedBar_ARNTExpr_vs_GradeStageClass.png"),
       plot = p_class, width = 6, height = 5)
