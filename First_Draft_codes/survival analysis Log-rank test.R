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

# 5) Fit Cox model for overall survival
cox_model <- coxph(
  Surv(time_to_overall_survival, mortality_binary) ~ 
    ARNT_expression + age_at_recruit_year + Gender,
  data = cleandfclinical_clean
)
summary(cox_model)

# 5) Dichotomize ARNT expression at the 75th percentile
q3 <- quantile(cleandfclinical$ARNT_expression, 0.75, na.rm = TRUE)
cleandfclinical[
  , ARNT_group := factor(
    ifelse(ARNT_expression >= q3, "High", "Low"),
    levels = c("Low","High")
  )
]


# 6) Create a progression flag from Prog_to_MIBC (1 if progressed, 0 otherwise)
cleandfclinical[
  , progression_binary := fifelse(Prog_to_MIBC %in% c("Yes","1"), 1L, 0L)
]

# 7) Define PFS_time and PFS_event
cleandfclinical[
  , PFS_time := fifelse(progression_binary == 1, time_to_prog, time_to_overall_survival)
]
cleandfclinical[
  , PFS_event := fifelse(progression_binary == 1 | mortality_binary == 1, 1L, 0L)
]

# 8) Subset to complete cases for PFS analysis
vars_pfs <- c("PFS_time", "PFS_event", "ARNT_group")

cleandfclinical_clean_pfs <- cleandfclinical[
  complete.cases(cleandfclinical[, ..vars_pfs]),
]

# 9) Run the log-rank test comparing PFS by ARNT_group
lr_test <- survdiff(
  Surv(PFS_time, PFS_event) ~ ARNT_group,
  data = cleandfclinical_clean_pfs
)
print(lr_test)

##Small sample exact test
# time    = PFS_time (numeric)
# status  = PFS_event (0/1)
# group   = ARNT_group (“Low”/“High” factor)
approx_lr <- coin::logrank_test(
  Surv(PFS_time, PFS_event) ~ ARNT_group,
  data         = cleandfclinical_clean_pfs,
  distribution = approximate(B = 10000)  # 10 000 reps
)
print(approx_lr)
pvalue(approx_lr)