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


#########################################################################
#########################################################################
#############################   Cleaning    #############################
#########################################################################
#########################################################################

for (col in names(dfClinical)) {
  if (is.factor(dfClinical[[col]]) || is.character(dfClinical[[col]])) {
    cat("\nColumn:", col, "\n")
    print(table(dfClinical[[col]], useNA = "ifany"))
  }
}


# Copy the original dataset
cleandfclinical <- dfClinical

# Replace 'Undeterminable' with NA in the grade column
cleandfclinical$grade[cleandfclinical$grade == "Undeterminable"] <- NA

# Replace 'pTis' and 'pTx' with NA in the stage column
cleandfclinical$stage[cleandfclinical$stage %in% c("pTis", "pTx")] <- NA

# Replace 'Unknown' with NA in the concomitant_in_situ column
cleandfclinical$concomitant_in_situ[cleandfclinical$concomitant_in_situ %in% c("Unknown")] <- NA

# Replace '99' with NA in the Smoking_status column
cleandfclinical$Smoking_status[cleandfclinical$Smoking_status %in% c("99")] <- NA

# Replace 'No' with NA in the bladder_ca column
cleandfclinical$bladder_ca[cleandfclinical$bladder_ca %in% c("No")] <- NA

table(cleandfclinical$bladder_ca, useNA = "ifany")

cleandfclinical$Recurrence[is.na(cleandfclinical$Recurrence)] <- "False"

## Convert age_at_recruit from days to years and round to 1 decimal
cleandfclinical[, age_at_recruit_year := age_at_recruit / 365.25 ]

summary(cleandfclinical[, c("number_tumours", "size_largest_tumour(cm)")])


#########################################################################
#########################################################################
######################   Merging Exprssion Data   #######################
#########################################################################
#########################################################################


# Step 1: Extract ARNT expression
dfARNT <- dfVeracyte[V1 == "ARNT"]

# Step 2: Get sample names and fix them
sample_names <- colnames(dfARNT)[-1]
fixed_names <- sub("\\.CEL$", "", sample_names)
fixed_names <- gsub("\\.", "-", fixed_names)
fixed_names <- paste0(fixed_names, ".CEL")

# Step 3: Transpose ARNT row
arnt_values <- as.numeric(dfARNT[, -1, with = FALSE])
dfARNT_t <- data.table(celfile_name = fixed_names, ARNT_expression = arnt_values)

# Merge on celfile_name
cleandfclinical <- merge(cleandfclinical, dfARNT_t, by = "celfile_name", all.x = TRUE)

##Save the dataset as .rds
#saveRDS(cleandfclinical, file = paste0(Address.Datatemp, "cleandfclinical.rds"))

fwrite(cleandfclinical, file = paste0(Address.Datatemp, "cleandfclinical.csv"))
