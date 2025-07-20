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

# 1. Compute tertile cut‐points
age_breaks <- quantile(
  cleandfclinical$age_at_recruit_year,
  probs = c(0, 1/3, 2/3, 1),
  na.rm = TRUE
)

# 2. Define the plain labels
tertile_labels <- c("Young", "Middle", "Old")

# 3. Cut into tertiles
cleandfclinical[
  , age_group := cut(
    age_at_recruit_year,
    breaks         = age_breaks,
    labels         = tertile_labels,
    include.lowest = TRUE
  )
]

# 4. Build the “label + (min–max)” text for each tertile
lower <- round(age_breaks[-length(age_breaks)], 0)
upper <- round(age_breaks[-1], 0)
range_text <- paste0("(", lower, "–", upper, ")")
# e.g. range_text = c("(25–45)", "(46–60)", "(61–80)")

full_labels <- paste0(tertile_labels, "\n", range_text)
# named vector so ggplot knows which label goes with which factor level
names(full_labels) <- tertile_labels

# 5. Define a colour palette
age_cols <- setNames(
  brewer.pal(n = length(tertile_labels), name = "Set3"),
  tertile_labels
)

# 6. Draw the boxplot
p_expr <- ggplot(cleandfclinical,
                 aes(x = age_group,
                     y = ARNT_expression,
                     fill = age_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  # apply our new labels and colours
  scale_x_discrete(labels = full_labels) +
  scale_fill_manual(values = age_cols,
                    labels = full_labels) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.8, color = "black") +
  theme_minimal() +
  labs(
    title = "ARNT Expression by Age Group",
    x     = "Age Group",
    y     = "ARNT Expression"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# 7. Kruskal–Wallis test
kw_res  <- kruskal.test(ARNT_expression ~ age_group,
                        data = cleandfclinical)
pval_kw <- format.pval(kw_res$p.value, digits = 3, eps = .001)

# 8. Add subtitle and show (or save)
p_expr <- p_expr +
  labs(subtitle = paste("Kruskal–Wallis p =", pval_kw))


# Save updated plot with p-value
ggsave(paste0(Address.OutputGraphs, "Veracyte_Boxplot_ARNTexpression_vs_age_group_Kruskal–Wallis.png"),
       plot = p_expr, width = 6, height = 5)





# 1. Prepare the data
df_model <- cleandfclinical[
  !is.na(ARNT_expression) & !is.na(age_group)
][
  , age_group := factor(age_group, levels = c("Young","Middle","Old"))
]

# 2. Fit a linear model predicting ARNT_expression by age_group and Gender
model <- lm(ARNT_expression ~ age_group + Gender, data = df_model)

# 3. Extract coefs and CIs
coef_summary <- summary(model)$coefficients; conf_int <- confint(model)

# 4. Tidy into data.table
dt <- as.data.table(coef_summary, keep.rownames = "Variable")
dt[, `:=`(
  CI_low   = conf_int[,1],
  CI_high  = conf_int[,2]
)][, `:=`(
  Estimate     = round(Estimate,    4),
  `Std. Error` = round(`Std. Error`,4),
  `t value`    = round(`t value`,   3),
  `Pr(>|t|)`   = signif(`Pr(>|t|)`, 3),
  CI_low       = round(CI_low,      3),
  CI_high      = round(CI_high,     3)
)]

# 5. Build LaTeX table
latex_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Linear regression of ARNT expression on age group and Gender}",
  "\\label{tab:lm_arnt_age_gender}",
  "\\begin{tabular}{lrrrrrrr}",
  "\\toprule",
  "Variable & Estimate & Std. Error & t value & Pr(>|t|) & CI low & CI high \\\\",
  "\\midrule",
  sapply(seq_len(nrow(dt)), function(i){
    r <- dt[i]
    sprintf(
      "%s & %.4f & %.4f & %.3f & %.3g & %.3f & %.3f \\\\",
      r$Variable, r$Estimate, r$`Std. Error`,
      r$`t value`, r$`Pr(>|t|)`, r$CI_low, r$CI_high
    )
  }),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

# 6. Write out in one line
writeLines(latex_lines, con = paste0(Address.OutputTables, "Veracyte_LinearRegression_ARNT_age_gender_Table.tex"))
