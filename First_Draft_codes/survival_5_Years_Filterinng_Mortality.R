#########################################################################
#########################################################################
##############################   Packages    ############################
#########################################################################
#########################################################################

# Load packages
# List required packages
install.packages("survminer")
packages <- c("data.table","ggplot2", "readxl", "survival", "survminer","coin")

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
########################   Survival Analysis    #########################
#########################################################################
#########################################################################

library(data.table)
setDT(cleandfclinical)

# 1) Create binary mortality status (1 = dead, 0 = alive)
cleandfclinical[
  , mortality_binary := fifelse(
    Mortality_status %in% c("Dead", "dead", "1"), 1L,
    fifelse(Mortality_status %in% c("Alive", "alive", "0"), 0L, NA_integer_)
  )
]

# 2) Clean up smoking status (blank → NA)
cleandfclinical[
  , Smoking_status := fifelse(Smoking_status == "", NA_character_, Smoking_status)
]

# 3) Convert Gender and Smoking_status to factors
cleandfclinical[
  , `:=`(
    Gender         = factor(Gender),
    Smoking_status = factor(Smoking_status)
  )
]

# 4) Subset to rows with no missing values in the model variables
vars <- c(
  "time_to_overall_survival",
  "mortality_binary",
  "ARNT_expression",
  "age_at_recruit_year",
  "Gender",
  "Smoking_status"
)
cleandfclinical_clean <- cleandfclinical[
  complete.cases(cleandfclinical[, ..vars]),
]

library(survival)


# Step 1: Filter data to include only patients who died within 5 years (1825 days) or are censored before that
cleandfclinical_os <- cleandfclinical_clean[
  time_to_overall_survival <= 1825
]

# Step 2: Dichotomize ARNT expression at the 75th percentile
q3 <- quantile(cleandfclinical_os$ARNT_expression, 0.75, na.rm = TRUE)
cleandfclinical_os[
  , ARNT_group := factor(
    ifelse(ARNT_expression >= q3, "High", "Low"),
    levels = c("Low", "High")
  )
]

# Step 3: Subset complete cases for OS analysis
vars_os <- c("time_to_overall_survival", "mortality_binary", "ARNT_group")
cleandfclinical_os <- cleandfclinical_os[
  complete.cases(cleandfclinical_os[, ..vars_os])
]

# Step 4: Standard log-rank test
lr_test_os <- survdiff(
  Surv(time_to_overall_survival, mortality_binary) ~ ARNT_group,
  data = cleandfclinical_os
)
print(lr_test_os)

# Step 5: Permutation-based log-rank test (small sample alternative)
approx_lr_os <- logrank_test(
  Surv(time_to_overall_survival, mortality_binary) ~ ARNT_group,
  data = cleandfclinical_os,
  distribution = approximate(B = 10000)
)
print(approx_lr_os)
pvalue(approx_lr_os)

# Step 6: Cox proportional hazards model
cox_os_model <- coxph(
  Surv(time_to_overall_survival, mortality_binary) ~ ARNT_group + age_at_recruit_year + Gender + Smoking_status,
  data = cleandfclinical_os
)
cox_os_sum <- summary(cox_os_model)
p_adj_os <- signif(cox_os_sum$coefficients["ARNT_groupHigh", "Pr(>|z|)"], 2)
HR_os    <- signif(cox_os_sum$coefficients["ARNT_groupHigh", "exp(coef)"], 2)

# Step 7: Kaplan-Meier plot
km_fit_os <- survfit(
  Surv(time_to_overall_survival, mortality_binary) ~ ARNT_group,
  data = cleandfclinical_os
)

pval_log_os <- surv_pvalue(
  km_fit_os,
  data = cleandfclinical_os
)$pval

km_plot_os <- ggsurvplot(
  km_fit_os,
  data         = cleandfclinical_os,
  legend.title = "ARNT Expression",
  legend.labs  = c("Low", "High"),
  risk.table   = TRUE,
  conf.int     = TRUE,
  palette      = c("steelblue", "firebrick")
)

# Annotate only log-rank p
max_time_os <- max(cleandfclinical_os$time_to_overall_survival, na.rm = TRUE)
annot_x_os  <- max_time_os * 0.25

km_plot_os$plot <- km_plot_os$plot +
  annotate("text",
           x     = annot_x_os, y = 0.20,
           label = paste0("log-rank p = ", signif(pval_log_os, 2)),
           hjust = 0)

# Save the plot
ggsave(
  filename = paste0(Address.OutputGraphs, "KM_ARNT_OS_High_Low_only_logrank_filtered_Mortalityonly.png"),
  plot     = km_plot_os$plot,
  width    = 6,
  height   = 5,
  dpi      = 300
)





# 1. Combine grade + stage into e.g. "G1pTa" or "G3pT1"
#    - sub() will replace "Grade_" (or "Grade ") with "G"
cleandfclinical_os[
  , grade_stage := paste0(
    sub("Grade[_ ]+", "G", grade),
    stage
  )
]

# 2. Create the new class variable and label only the two combos
cleandfclinical_os[
  , grade_stage_class := NA_character_
][
  grade_stage == "G1pTa", grade_stage_class := "Low"
][
  grade_stage == "G3pT1", grade_stage_class := "High"
]

# 3. Turn into a factor with the correct order
cleandfclinical_os[
  , grade_stage_class := factor(grade_stage_class, levels = c("Low", "High"))
]

# 4. Filter to just those two groups
cleandfclinical_os <- cleandfclinical_os[!is.na(grade_stage_class)]

# 0. Set factor levels and labels
cleandfclinical_os$grade_stage_class <- factor(
  cleandfclinical_os$grade_stage_class,
  levels = c("Low", "High"),
  labels = c("Low (G1pTa)", "High (G3pT1)")
)








# 1. Fit the Cox model for overall survival
cox_os_model <- coxph(
  Surv(time_to_overall_survival, mortality_binary) ~ ARNT_group + Gender + Smoking_status + grade_stage_class,
  data = cleandfclinical_os
)

# 2. Get the summary
cox_os_summary <- summary(cox_os_model)

# 3. Create a data frame for the LaTeX table
cox_os_table <- data.frame(
  Variable   = rownames(cox_os_summary$coefficients),
  HR         = signif(cox_os_summary$coefficients[, "exp(coef)"], 2),
  `p-value`  = signif(cox_os_summary$coefficients[, "Pr(>|z|)"], 2),
  check.names = FALSE
)

# 4. Convert to LaTeX format
latex_os <- toLatex(
  xtable::xtable(cox_os_table,
                 caption = "Cox proportional hazards model for overall survival",
                 label = "tab:cox_os"),
  include.rownames = FALSE
)

# 5. Save the LaTeX table to the output folder
writeLines(latex_os, con = paste0(Address.OutputTables, "cox_os_summary_5_Years_Mortality_table_GS.txt"))
