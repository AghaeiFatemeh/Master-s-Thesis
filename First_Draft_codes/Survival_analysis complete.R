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



# 5) Dichotomize ARNT expression at the 75th percentile
q3 <- quantile(cleandfclinical_clean$ARNT_expression, 0.75, na.rm = TRUE)
cleandfclinical_clean[
  , ARNT_group := factor(
    ifelse(ARNT_expression >= q3, "High", "Low"),
    levels = c("Low","High")
  )
]




# 6) Create binary progression variable (1 if TRUE, 0 if FALSE or NA)
cleandfclinical_clean[
  , progression_binary := fifelse(Prog_to_MIBC == TRUE, 1L, 0L)
]

# 7) Define PFS_time
cleandfclinical_clean[
  , PFS_time := fifelse(
    progression_binary == 1 & !is.na(time_to_prog),
    time_to_prog,
    time_to_overall_survival
  )
]

# 8) Define PFS_event
cleandfclinical_clean[
  , PFS_event := fifelse(
    (progression_binary == 1 & !is.na(time_to_prog)) | mortality_binary == 1,
    1L,
    0L
  )
]

# 9) Subset to complete cases for PFS analysis
vars_pfs <- c("PFS_time", "PFS_event", "ARNT_group")
cleandfclinical_clean_pfs <- cleandfclinical_clean[
  complete.cases(cleandfclinical_clean[, ..vars_pfs])
]


# 10) Run the log-rank test comparing PFS by ARNT_group
library(survival)
lr_test <- survdiff(
  Surv(PFS_time, PFS_event) ~ ARNT_group,
  data = cleandfclinical_clean_pfs
)
print(lr_test)

# 11) Small sample log-rank test (permutation-based)
library(coin)
approx_lr <- logrank_test(
  Surv(PFS_time, PFS_event) ~ ARNT_group,
  data         = cleandfclinical_clean_pfs,
  distribution = approximate(B = 10000)
)
print(approx_lr)
pvalue(approx_lr)





# 1) Fit multivariate Cox for PFS to get the adjusted p-value and HR for ARNT
cox_pfs_model <- coxph(
  Surv(PFS_time, PFS_event) ~ ARNT_group + age_at_recruit_year + Gender + Smoking_status,
  data = cleandfclinical_clean_pfs
)
cox_sum <- summary(cox_pfs_model)

# Extract the adjusted p-value and HR for ARNT_groupHigh
p_adj <- signif(cox_sum$coefficients["ARNT_groupHigh", "Pr(>|z|)"], 2)
HR    <- signif(cox_sum$coefficients["ARNT_groupHigh", "exp(coef)"], 2)

# 2) Fit Kaplan-Meier and get log-rank p-value
km_fit <- survfit(
  Surv(PFS_time, PFS_event) ~ ARNT_group,
  data = cleandfclinical_clean_pfs
)

pval_log <- surv_pvalue(
  km_fit,
  data = cleandfclinical_clean_pfs
)$pval  # log-rank p-value

# 3) Create the Kaplan-Meier plot
km_plot <- ggsurvplot(
  km_fit,
  data         = cleandfclinical_clean_pfs,
  legend.title = "ARNT Expression",
  legend.labs  = c("Low", "High"),
  risk.table   = TRUE,
  conf.int     = TRUE,
  palette      = c("steelblue", "firebrick")
)

# 4) Annotate both p-values and HR onto the plot
max_time <- max(cleandfclinical_clean_pfs$PFS_time, na.rm = TRUE)
annot_x  <- max_time * 0.25

km_plot$plot <- km_plot$plot +
  annotate("text",
           x     = annot_x, y = 0.20,
           label = paste0("log-rank p = ", signif(pval_log, 2)),
           hjust = 0)

# Save the plot
ggsave(
  filename = paste0(Address.OutputGraphs, "KM_ARNT_PFS_High_Low_with_logrank_only_all age group.png"),
  plot     = km_plot$plot,
  width    = 6,
  height   = 5,
  dpi      = 300
)