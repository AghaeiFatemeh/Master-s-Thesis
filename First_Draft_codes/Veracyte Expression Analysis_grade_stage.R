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
###################   ARNT Expression vs Mortality    ###################
#########################################################################
#########################################################################

# 1. Combine grade + stage into e.g. "G1pTa" or "G3pT1"
#    - sub() will replace "Grade_" (or "Grade ") with "G"
cleandfclinical[
  , grade_stage := paste0(
    sub("Grade[_ ]+", "G", grade),
    stage
  )
]

# 2. Create the new class variable and label only the two combos
cleandfclinical[
  , grade_stage_class := NA_character_
][
  grade_stage == "G1pTa", grade_stage_class := "Low"
][
  grade_stage == "G3pT1", grade_stage_class := "High"
]

# 3. Turn into a factor with the correct order
cleandfclinical[
  , grade_stage_class := factor(grade_stage_class, levels = c("Low", "High"))
]

# 4. Filter to just those two groups
df_labelled <- cleandfclinical[!is.na(grade_stage_class)]

# 0. Set factor levels and labels
df_labelled$grade_stage_class <- factor(
  df_labelled$grade_stage_class,
  levels = c("Low", "High"),
  labels = c("Low (G1pTa)", "High (G3pT1)")
)

# 1. Draw the boxplot with jittered sample points and reversed colors
p_expr <- ggplot(df_labelled,
                 aes(x = grade_stage_class,
                     y = ARNT_expression,
                     fill = grade_stage_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.8, color = "black") +
  theme_minimal() +
  labs(
    title = "ARNT Expression by Grade–Stage Class",
    subtitle = paste("Wilcoxon p =", format.pval(
      wilcox.test(ARNT_expression ~ grade_stage_class, data = df_labelled, exact = FALSE)$p.value,
      digits = 3, eps = .001)),
    x    = "Grade–Stage Class",
    y    = "ARNT Expression",
    fill = "Class"
  ) +
  scale_fill_manual(
    values = c("Low (G1pTa)" = "#9ecae1", "High (G3pT1)" = "#fcae91")  # reversed colors
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )

# 2. Save the plot
ggsave(
  filename = file.path(Address.OutputGraphs,
                       "Veracyte_Boxplot_ARNTexpression_vs_GradeStageClass_withPval.png"),
  plot   = p_expr,
  width  = 6,
  height = 5
)


# 1. Filter and prepare data
df_model <- df_labelled[
  !is.na(ARNT_expression) & !is.na(grade_stage_class)
][
  , grade_stage_class := factor(grade_stage_class, levels = c("Low","High"))
][
  , grade_stage_class := relevel(grade_stage_class, ref = "Low")
]

# 2. Run logistic regression
model <- glm(
  grade_stage_class ~ ARNT_expression + Gender + age_at_recruit_year,
  data   = df_model,
  family = binomial
)

# 3. Extract summary and CIs
coef_summary <- summary(model)$coefficients; conf_int <- confint(model)

# 4. Convert to data.table
dt_table <- as.data.table(coef_summary, keep.rownames = "Variable")

# 5. Add OR and 95% CI
dt_table[
  , OddsRatio := exp(coef(model))
][
  , CI_low    := exp(conf_int[,1])
][
  , CI_high   := exp(conf_int[,2])
]

# 6. Round values
dt_table[, `:=`(
  Estimate     = round(Estimate,    4),
  `Std. Error` = round(`Std. Error`,4),
  `z value`    = round(`z value`,   3),
  `Pr(>|z|)`   = signif(`Pr(>|z|)`, 3),
  OddsRatio    = round(OddsRatio,   3),
  CI_low       = round(CI_low,      3),
  CI_high      = round(CI_high,     3)
)]

# 7. Build LaTeX-formatted lines
latex_lines <- c(
  "\\begin{table}[h]",
  "\\centering",
  "\\begin{tabular}{lrrrrrrr}",
  "\\hline",
  "Variable & Estimate & Std. Error & z value & Pr(>|z|) & Odds Ratio & CI Low & CI High \\\\",
  "\\hline",
  sapply(
    seq_len(nrow(dt_table)),
    function(i) paste0(
      dt_table$Variable[i], " & ",
      dt_table$Estimate[i], " & ",
      dt_table$`Std. Error`[i], " & ",
      dt_table$`z value`[i], " & ",
      dt_table$`Pr(>|z|)`[i], " & ",
      dt_table$OddsRatio[i], " & ",
      dt_table$CI_low[i], " & ",
      dt_table$CI_high[i], " \\\\\\\\"
    )
  ),
  "\\hline",
  "\\end{tabular}",
  "\\caption{Logistic regression results for ARNT expression predicting grade–stage class (Low vs High).}",
  "\\label{tab:logistic_arnt_grade_stage_class}",
  "\\end{table}"
)

# 8. Write to .tex file in one line
writeLines(latex_lines, con = paste0(Address.OutputTables, "VeracyteLogisticRegression_ARNT_grade_stage_class_Table.tex"))
