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

#########################################################################
#########################################################################
##########################   Distribution    ############################
#########################################################################
#########################################################################

#Plot distribution of each column

#Gender
ggplot(cleandfclinical, aes(x = Gender)) +
  geom_bar(fill = "skyblue", color = "black", na.rm = FALSE) +
  theme_minimal() +
  ggtitle("Gender Distribution") +
  xlab("Gender") +
  ylab("Count") +
  ggsave(paste0(Address.OutputGraphs, "GenderBarPlot2.png"), width = 10, height = 6)

#Smoking_status
ggplot(cleandfclinical, aes(x = Smoking_status)) +
  geom_bar(fill = "skyblue", color = "black", na.rm = FALSE) +
  theme_minimal() +
  ggtitle("Smoking Status Distribution") +
  xlab("Smoking Status") +
  ylab("Count") +
  ggsave(paste0(Address.OutputGraphs, "SmokingStatusBarPlot2.png"), width = 10, height = 6)

#grade
ggplot(cleandfclinical, aes(x = grade)) +
  geom_bar(fill = "skyblue", color = "black", na.rm = FALSE) +
  theme_minimal() +
  ggtitle("Grade Distribution") +
  xlab("Grade") +
  ylab("Count") +
  ggsave(paste0(Address.OutputGraphs, "GradeBarPlot2.png"), width = 10, height = 6)

#stage
ggplot(cleandfclinical, aes(x = stage)) +
  geom_bar(fill = "skyblue", color = "black", na.rm = FALSE) +
  theme_minimal() +
  ggtitle("Stage Distribution") +
  xlab("Stage") +
  ylab("Count") +
  ggsave(paste0(Address.OutputGraphs, "StageBarPlot2.png"), width = 10, height = 6)

#concomitant_in_situ
ggplot(cleandfclinical, aes(x = concomitant_in_situ)) +
  geom_bar(fill = "skyblue", color = "black", na.rm = FALSE) +
  theme_minimal() +
  ggtitle("Concomitant In Situ Distribution") +
  xlab("Concomitant In Situ") +
  ylab("Count") +
  ggsave(paste0(Address.OutputGraphs, "ConcomitantInSituBarPlot.png"), width = 10, height = 6)

#Mortality_status
ggplot(cleandfclinical, aes(x = Mortality_status)) +
  geom_bar(fill = "skyblue", color = "black", na.rm = FALSE) +
  theme_minimal() +
  ggtitle("Mortality Status Distribution") +
  xlab("Mortality Status") +
  ylab("Count") +
  ggsave(paste0(Address.OutputGraphs, "MortalityStatusBarPlot2.png"), width = 10, height = 6)

#Recurrence
ggplot(cleandfclinical, aes(x = Recurrence)) +
  geom_bar(fill = "skyblue", color = "black", na.rm = FALSE) +
  theme_minimal() +
  ggtitle("Recurrence Distribution") +
  xlab("Recurrence") +
  ylab("Count") +
  ggsave(paste0(Address.OutputGraphs, "RecurrenceBarPlot2.png"), width = 10, height = 6)

#age_at_recruit_year
ggplot(cleandfclinical, aes(x = age_at_recruit_year)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  theme_minimal() +
  ggtitle("Age Distribution") +
  ggsave(paste0(Address.OutputGraphs, "AgeHistVerocyte2.png"), width = 10, height = 6)


#number_tumours
ggplot(cleandfclinical, aes(x = number_tumours)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  theme_minimal() +
  ggtitle("number_tumours Distribution") +
  ggsave(paste0(Address.OutputGraphs, "number_tumoursHistVerocyte2.png"), width = 10, height = 6)

#size_largest_tumour(cm)
ggplot(cleandfclinical, aes(x = `size_largest_tumour(cm)`)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  theme_minimal() +
  ggtitle("Size of Largest Tumour Distribution") +
  ggsave(paste0(Address.OutputGraphs, "size_largest_tumour_Hist.png"), width = 10, height = 6)

