#########################################################################
#########################################################################
##############################   Packages    ############################
#########################################################################
#########################################################################

# Load packages
# List required packages
packages <- c("data.table")

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



# Read the large CSV file
dfVeracyte <- fread(paste0(Address.DataRaw, "veracyte_nmibc_transcriptomes.csv"))

# Check dimensions and column names
dim(dfVeracyte)
head(dfVeracyte)




# Filter the ARNT row
arnt_row <- dfVeracyte[dfVeracyte$V1 == "ARNT", ]

# Transpose the expression values (excluding the first column)
arnt_expr <- as.numeric(unlist(arnt_row[, -1, with = FALSE]))

# Create a data frame for plotting
dfARNT <- data.frame(Expression = arnt_expr)


#########################################################################
#########################################################################
##########################   Distribution    ############################
#########################################################################
#########################################################################



#Density plot of expression distributions all toegther
# Plot all expression values together (no sample-wise coloring)
p <- ggplot(dfARNT, aes(x = Expression)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  theme_minimal() +
  ggtitle("Overall Density of Veracyte Expression Values") +
  xlab("Expression") +
  ylab("Density") +
  ggsave(paste0(Address.OutputGraphs, "VeracyteExpressionDistributionTogether.png"), width = 10, height = 6)
