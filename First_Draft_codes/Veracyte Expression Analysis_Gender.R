#########################################################################
#########################################################################
##############################   Packages    ############################
#########################################################################
#########################################################################

# Load packages
# List required packages
packages <- c("data.table","ggplot2", "readxl")

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


cleandfclinical <- cleandfclinical[!is.na(Gender)]
# Remove empty strings and make it a proper factor
cleandfclinical <- cleandfclinical[Gender %in% c("male", "female")]
cleandfclinical$Gender <- factor(cleandfclinical$Gender)



cleandfclinical <- cleandfclinical[!is.na(Gender)]
# Remove empty strings and make it a proper factor
cleandfclinical <- cleandfclinical[Gender %in% c("male", "female")]
cleandfclinical$Gender <- factor(cleandfclinical$Gender)

# Boxplot: ARNT Expression vs Gender
p_expr <- ggplot(cleandfclinical, aes(x = Gender, y = ARNT_expression, fill = Gender)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.8, color = "black") +
  theme_minimal() +
  labs(title = "ARNT Expression by Gender",
       x = "Gender", y = "ARNT Expression") +
  scale_fill_manual(values = c("male" = "#fcae91", "female" = "#9ecae1")) +
  theme(plot.title = element_text(hjust = 0.5))

# Wilcoxon test (non-parametric for 2 groups)
test_expr <- wilcox.test(ARNT_expression ~ Gender, data = cleandfclinical, exact = FALSE)
pval_expr <- format.pval(test_expr$p.value, digits = 3, eps = .001)

# Add p-value to plot
p_expr <- p_expr + labs(subtitle = paste("Wilcoxon p =", pval_expr))

# Save updated plot with p-value
ggsave(paste0(Address.OutputGraphs, "Veracyte_Boxplot_ARNTexpression_vs_Gender_wilcoxon.png"),
       plot = p_expr, width = 6, height = 5)


# 1. Filter and prepare data
df_model <- cleandfclinical[!is.na(ARNT_expression) & !is.na(Gender)]
df_model$Gender <- factor(df_model$Gender)
df_model$Gender <- relevel(df_model$Gender, ref = "female")

# 2. Run logistic regression
model <- glm(Gender ~ ARNT_expression + age_at_recruit_year, data = df_model, family = binomial)

# 3. Extract summary
coef_summary <- summary(model)$coefficients
conf_int <- confint(model)

# 4. Convert to data.table
dt_table <- as.data.table(coef_summary, keep.rownames = "Variable")

# 5. Add OR and 95% CI
dt_table[, OddsRatio := exp(coef(model))]
dt_table[, CI_low := exp(conf_int[, 1])]
dt_table[, CI_high := exp(conf_int[, 2])]

# 6. Round values
dt_table[, `:=`(
  Estimate = round(Estimate, 4),
  `Std. Error` = round(`Std. Error`, 4),
  `z value` = round(`z value`, 3),
  `Pr(>|z|)` = signif(`Pr(>|z|)`, 3),
  OddsRatio = round(OddsRatio, 3),
  CI_low = round(CI_low, 3),
  CI_high = round(CI_high, 3)
)]

# 7. Create LaTeX-formatted string
latex_lines <- c(
  "\\begin{table}[h]",
  "\\centering",
  "\\begin{tabular}{lrrrrrrr}",
  "\\hline",
  "Variable & Estimate & Std. Error & z value & Pr(>|z|) & Odds Ratio & CI Low & CI High \\\\",
  "\\hline"
)

# Add each row of the table
for (i in 1:nrow(dt_table)) {
  row <- dt_table[i]
  line <- paste0(
    row$Variable, " & ",
    row$Estimate, " & ",
    row$`Std. Error`, " & ",
    row$`z value`, " & ",
    row$`Pr(>|z|)`, " & ",
    row$OddsRatio, " & ",
    row$CI_low, " & ",
    row$CI_high, " \\\\"
  )
  latex_lines <- c(latex_lines, line)
}

# Finish table
latex_lines <- c(
  latex_lines,
  "\\hline",
  "\\end{tabular}",
  "\\caption{Logistic regression results for ARNT expression predicting Gender}",
  "\\label{tab:logistic_arnt_Gender}",
  "\\end{table}"
)

# 8. Save as .tex file
writeLines(latex_lines, con = paste0(Address.OutputTables, "Veracyte_LogisticRegression_ARNT_Gender_Table.tex"))
