#########################################################################
#########################################################################
##############################   Packages    ############################
#########################################################################
#########################################################################

# Load packages
# List required packages
install.packages("survminer")
packages <- c("data.table","ggplot2", "readxl", "survival", "survminer")

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


cleandfclinical[,
                mortality_binary := fifelse(
                  Mortality_status %in% c("Dead","dead","1"), 1L,
                  fifelse(Mortality_status %in% c("Alive","alive","0"), 0L, NA_integer_)
                )
]



cleandfclinical[Smoking_status == "", Smoking_status := NA_character_]



cleandfclinical[, `:=`(
  Gender         = factor(Gender),
  Smoking_status = factor(Smoking_status)
)]




vars <- c("time_to_overall_survival",
          "mortality_binary",
          "ARNT_expression",
          "age_at_recruit_year",
          "Gender",
          "Smoking_status")

cleandfclinical_clean <- cleandfclinical[complete.cases(cleandfclinical[, ..vars])]




cox_model <- coxph(
  Surv(time_to_overall_survival, mortality_binary) ~ 
    ARNT_expression + age_at_recruit_year + Gender,
  data = cleandfclinical_clean
)

summary(cox_model)





# Add ARNT expression group: High if above 75th percentile, Low otherwise
q3 <- quantile(cleandfclinical_clean$ARNT_expression, 0.75, na.rm = TRUE)

cleandfclinical_clean[, ARNT_group := ifelse(ARNT_expression >= q3, "High", "Low")]
cleandfclinical_clean[, ARNT_group := factor(ARNT_group, levels = c("Low", "High"))]


# Fit survival object
surv_obj <- Surv(cleandfclinical_clean$time_to_overall_survival,
                 cleandfclinical_clean$mortality_binary)

# Plot
km_fit <- survfit(surv_obj ~ ARNT_group, data = cleandfclinical_clean)



ggsave(
  filename = paste0(Address.OutputGraphs, "KM_ARNT_High_Low_survival.png"),
  plot = ggsurvplot(km_fit,
                    data = cleandfclinical_clean,
                    pval = TRUE,
                    conf.int = TRUE,
                    risk.table = TRUE,
                    legend.title = "ARNT Expression",
                    legend.labs = c("Low", "High"),
                    palette = c("steelblue", "firebrick"))$plot,
  width = 6,
  height = 5,
  dpi = 300
)




cox_cat <- coxph(Surv(time_to_overall_survival, mortality_binary) ~ 
                   ARNT_group + age_at_recruit_year + Gender + Smoking_status,
                 data = cleandfclinical_clean)
summary(cox_cat)
