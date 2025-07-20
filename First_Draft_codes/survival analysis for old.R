# 2) Fit Kaplan-Meier and get log-rank p-value
km_fit <- survfit(
  Surv(PFS_time, PFS_event) ~ ARNT_group,
  data = cleandfclinical_clean_pfs[age_at_recruit_year>=75]
)

pval_log <- surv_pvalue(
  km_fit,
  data = cleandfclinical_clean_pfs[age_at_recruit_year>=75]
)$pval  # log-rank p-value

# 3) Create the Kaplan-Meier plot
km_plot <- ggsurvplot(
  km_fit,
  data         = cleandfclinical_clean_pfs[age_at_recruit_year>=75],
  legend.title = "ARNT Expression",
  legend.labs  = c("Low", "High"),
  risk.table   = TRUE,
  conf.int     = TRUE,
  palette      = c("steelblue", "firebrick")
)

# 4) Annotate both p-values and HR onto the plot
max_time <- max(cleandfclinical_clean_pfs[age_at_recruit_year>=75]$PFS_time, 
                na.rm = TRUE)
annot_x  <- max_time * 0.25

km_plot$plot <- km_plot$plot +
  annotate("text",
           x     = annot_x, y = 0.20,
           label = paste0("log-rank p = ", signif(pval_log, 2)),
           hjust = 0)

# Save the plot
ggsave(
  filename = paste0(Address.OutputGraphs, "KM_ARNT_PFS_High_Low_with_logrank_only_old.png"),
  plot     = km_plot$plot,
  width    = 6,
  height   = 5,
  dpi      = 300
)