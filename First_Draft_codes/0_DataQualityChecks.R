#########################################################################
#########################################################################
##############################   Packages    ############################
#########################################################################
#########################################################################

# Load packages
# List required packages
packages <- c("data.table","ggplot2")

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
###############################   Reading    ##############################
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
########################   CheckNormalization    ########################
#########################################################################
#########################################################################

#Reshape the data for visualization

# Rename first column to "Gene"
setnames(dfRNACounts, 1, "Gene")

# Melt the data into long format for plotting
dfRNACountsLong <- melt(dfRNACounts, id.vars = "Gene",
                        variable.name = "Sample",
                        value.name = "Expression")

# Quick check
head(dfRNACountsLong)

#Density plot of expression distributions per sample

ggplot(dfRNACountsLong, aes(x = Expression, color = Sample)) +
  geom_density() +
  theme_minimal() +
  ggtitle("Density plot of normalized expression values") +
  xlab("Expression") +
  ylab("Density") +
  ggsave(paste0(Address.OutputGraphs,"ExpressionDistribution.png"), width = 10, height = 6)

#Density plot of expression distributions all toegther
# Plot all expression values together (no sample-wise coloring)
p <- ggplot(dfRNACountsLong, aes(x = Expression)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  theme_minimal() +
  ggtitle("Overall Density of Normalized Expression Values") +
  xlab("Expression") +
  ylab("Density") +
  ggsave(paste0(Address.OutputGraphs, "ExpressionDistributionTogether.png"), width = 10, height = 6)

summary(dfRNACountsLong[,.(Expression)])

#########################################################################
#########################################################################
########################   CheckDistribution    #########################
#########################################################################
#########################################################################

#Checking unique values for each column
sapply(dfSampleInfo, function(x) length(unique(x)))
#Cheking missing values
colSums(is.na(dfSampleInfo))

#Plot distribution of each column

ggplot(dfSampleInfo, aes(x = Age)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  theme_minimal() +
  ggtitle("Age Distribution") +
  ggsave(paste0(Address.OutputGraphs, "AgeHist.png"), width = 10, height = 6)

ggplot(dfSampleInfo, aes(y = Age)) +
  geom_boxplot(fill = "orange") +
  theme_minimal() +
  ggtitle("Age Boxplot (Check for Outliers)") +
  ggsave(paste0(Address.OutputGraphs, "AgeBoxPlot.png"), width = 10, height = 6)


ggplot(dfSampleInfo, aes(x = Stage)) +
  geom_bar(fill = "lightgreen") +
  theme_minimal() +
  ggtitle("Stage Counts") +
  ggsave(paste0(Address.OutputGraphs, "StageBar.png"), width = 10, height = 6)

ggplot(dfSampleInfo, aes(x = Grade)) +
  geom_bar(fill = "plum") +
  theme_minimal() +
  ggtitle("Grade Counts") +
  ggsave(paste0(Address.OutputGraphs, "GradeBar.png"), width = 10, height = 6)

ggplot(dfSampleInfo, aes(x = ARNT)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  ggtitle("ARNT Mutation Status") +
  ggsave(paste0(Address.OutputGraphs, "ARNTBar.png"), width = 10, height = 6)

ggplot(dfSampleInfo, aes(x = Gender)) +
  geom_bar(fill = "lightgreen") +
  theme_minimal() +
  ggtitle("Stage Counts") +
  ggsave(paste0(Address.OutputGraphs, "GenderBar.png"), width = 10, height = 6)

#########################################################################
#########################################################################
##############################   Packages    ############################
#########################################################################
#########################################################################

# Load packages
# List required packages
packages <- c("data.table","ggplot2")

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
###############################   Reading    ##############################
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
########################   CheckNormalization    ########################
#########################################################################
#########################################################################

#Reshape the data for visualization

# Rename first column to "Gene"
setnames(dfRNACounts, 1, "Gene")

# Melt the data into long format for plotting
dfRNACountsLong <- melt(dfRNACounts, id.vars = "Gene",
                        variable.name = "Sample",
                        value.name = "Expression")

# Quick check
head(dfRNACountsLong)

#Density plot of expression distributions per sample

ggplot(dfRNACountsLong, aes(x = Expression, color = Sample)) +
  geom_density() +
  theme_minimal() +
  ggtitle("Density plot of normalized expression values") +
  xlab("Expression") +
  ylab("Density") +
  ggsave(paste0(Address.OutputGraphs,"ExpressionDistribution.png"), width = 10, height = 6)

#Density plot of expression distributions all toegther
# Plot all expression values together (no sample-wise coloring)
p <- ggplot(dfRNACountsLong, aes(x = Expression)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  theme_minimal() +
  ggtitle("Overall Density of Normalized Expression Values") +
  xlab("Expression") +
  ylab("Density") +
  ggsave(paste0(Address.OutputGraphs, "ExpressionDistributionTogether.png"), width = 10, height = 6)

summary(dfRNACountsLong[,.(Expression)])

#########################################################################
#########################################################################
########################   CheckDistribution    #########################
#########################################################################
#########################################################################

#Checking unique values for each column
sapply(dfSampleInfo, function(x) length(unique(x)))
#Cheking missing values
colSums(is.na(dfSampleInfo))

#Plot distribution of each column

ggplot(dfSampleInfo, aes(x = Age)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  theme_minimal() +
  ggtitle("Age Distribution") +
  ggsave(paste0(Address.OutputGraphs, "AgeHist.png"), width = 10, height = 6)

ggplot(dfSampleInfo, aes(y = Age)) +
  geom_boxplot(fill = "orange") +
  theme_minimal() +
  ggtitle("Age Boxplot (Check for Outliers)") +
  ggsave(paste0(Address.OutputGraphs, "AgeBoxPlot.png"), width = 10, height = 6)


ggplot(dfSampleInfo, aes(x = Stage)) +
  geom_bar(fill = "lightgreen") +
  theme_minimal() +
  ggtitle("Stage Counts") +
  ggsave(paste0(Address.OutputGraphs, "StageBar.png"), width = 10, height = 6)

ggplot(dfSampleInfo, aes(x = Grade)) +
  geom_bar(fill = "plum") +
  theme_minimal() +
  ggtitle("Grade Counts") +
  ggsave(paste0(Address.OutputGraphs, "GradeBar.png"), width = 10, height = 6)

ggplot(dfSampleInfo, aes(x = ARNT)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  ggtitle("ARNT Mutation Status") +
  ggsave(paste0(Address.OutputGraphs, "ARNTBar.png"), width = 10, height = 6)

ggplot(dfSampleInfo, aes(x = Smoking_Status)) +
  geom_bar(fill = "lightgreen") +
  theme_minimal() +
  ggtitle("Stage Counts") +
  ggsave(paste0(Address.OutputGraphs, "SmokingStatusBar.png"), width = 10, height = 6)

ggplot(dfSampleInfo, aes(x = RNA_class)) +
  geom_bar(fill = "lightgreen") +
  theme_minimal() +
  ggtitle("Stage Counts") +
  ggsave(paste0(Address.OutputGraphs, "RNAClassBar.png"), width = 10, height = 6)

ggplot(dfSampleInfo, aes(x = UROMOL_NMIBC_Class)) +
  geom_bar(fill = "lightgreen") +
  theme_minimal() +
  ggtitle("Stage Counts") +
  ggsave(paste0(Address.OutputGraphs, "UROMOLNMIBCClassBar.png"), width = 10, height = 6)

ggplot(dfSampleInfo, aes(x = Gender)) +
  geom_bar(fill = "lightgreen") +
  theme_minimal() +
  ggtitle("Stage Counts") +
  ggsave(paste0(Address.OutputGraphs, "GenderBar.png"), width = 10, height = 6)

#Compare ARNT groups
ggplot(dfSampleInfo, aes(x = ARNT, y = Age)) +
  geom_boxplot() +
  ggtitle("Age by ARNT Mutation Status") +
  ggsave(paste0(Address.OutputGraphs, "AgeARNTBoxPlot.png"), width = 10, height = 6)
