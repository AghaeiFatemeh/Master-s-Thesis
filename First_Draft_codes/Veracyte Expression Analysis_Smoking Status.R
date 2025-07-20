#########################################################################
#########################################################################
##############################   Packages    ############################
#########################################################################
#########################################################################

# Load packages
# List required packages
packages <- c("data.table","ggplot2", "readxl", "RColorBrewer", "nnet", "stats")

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

# Remove rows where Smoking_status is NA or blank
cleandfclinical <- cleandfclinical[!is.na(Smoking_status) & Smoking_status != ""]



# 1. Define smoking‐status levels, labels and colours
smoke_levels <- c("Yes", "No, but I used to smoke", "No, I have never smoked", "")
smoke_labels <- c("Yes", "Former smoker", "Never smoked", "Unknown")
smoke_cols   <- setNames(
  brewer.pal(n = length(smoke_levels), name = "Set3"),
  smoke_levels
)

# 2. Draw the boxplot by Smoking_status
p_expr <- ggplot(cleandfclinical, 
                 aes(x = Smoking_status, 
                     y = ARNT_expression, 
                     fill = Smoking_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  scale_x_discrete(labels = smoke_labels) +
  scale_fill_manual(values = smoke_cols, labels = smoke_labels) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.8, color = "black") +
  theme_minimal() +
  labs(
    title = "ARNT Expression by Smoking Status",
    x     = "Smoking Status",
    y     = "ARNT Expression"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# 3. Kruskal–Wallis test
kw_res  <- kruskal.test(ARNT_expression ~ Smoking_status,
                        data = cleandfclinical)
pval_kw <- format.pval(kw_res$p.value, digits = 3, eps = .001)

# 4. Add subtitle and save
p_expr <- p_expr + labs(subtitle = paste("Kruskal–Wallis p =", pval_kw))

ggsave(
  filename = paste0(Address.OutputGraphs, 
                    "Veracyte_Boxplot_ARNTexpression_vs_SmokingStatus_Kruskal.png"),
  plot   = p_expr,
  width  = 6,
  height = 5
)




# 5. Prepare data for linear model
df_model <- cleandfclinical[
  !is.na(ARNT_expression) & 
    Smoking_status != "" & 
    !is.na(Smoking_status)
]
df_model$Smoking_status <- factor(df_model$Smoking_status, levels = smoke_levels)
df_model$Smoking_status <- relevel(df_model$Smoking_status, ref = "No, I have never smoked")

# 6. Fit linear regression
model <- lm(ARNT_expression ~ Smoking_status + Gender + age_at_recruit_year, data = df_model)

# 7. Extract coefficients and confidence intervals
coef_summary <- summary(model)$coefficients
conf_int     <- confint(model)

# 8. Tidy into a data.table
library(data.table)
dt <- as.data.table(coef_summary, keep.rownames = "Variable")
dt[, `:=`(
  CI_low  = conf_int[,1],
  CI_high = conf_int[,2]
)]
dt[, `:=`(
  Estimate     = round(Estimate,    4),
  `Std. Error` = round(`Std. Error`,4),
  `t value`    = round(`t value`,   3),
  `Pr(>|t|)`   = signif(`Pr(>|t|)`, 3),
  CI_low       = round(CI_low,      3),
  CI_high      = round(CI_high,     3)
)]

# 9. Build LaTeX table
latex_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Linear regression of ARNT expression on Smoking status and Gender}",
  "\\label{tab:lm_arnt_smoking_gender}",
  "\\begin{tabular}{lrrrrrrr}",
  "\\toprule",
  "Variable & Estimate & Std. Error & t value & Pr(>|t|) & CI low & CI high \\\\",
  "\\midrule",
  sapply(seq_len(nrow(dt)), function(i) {
    r <- dt[i]
    sprintf(
      "%s & %.4f & %.4f & %.3f & %.3g & %.3f & %.3f \\\\",
      r$Variable, r$Estimate, r$`Std. Error`,
      r$`t value`,  r$`Pr(>|t|)`, r$CI_low, r$CI_high
    )
  }),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

# 10. Save as .tex file
writeLines(
  latex_lines, 
  con = paste0(Address.OutputTables, 
               "Veracyte_LinearRegression_ARNT_SmokingStatus_Gender_Table.tex")
)
