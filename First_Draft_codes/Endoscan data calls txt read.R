# load data.table
library(data.table)

# define your Samples directory
sample_dir <- "/rds/projects/a/arnoldry-arnt/project_Fatemeh_ARNT/oncoscan_data/Samples"

# find all calls.txt files (recursive search, full paths)
file_paths <- list.files(
  path       = sample_dir,
  pattern    = "^calls\\.txt$",
  recursive  = TRUE,
  full.names = TRUE
)

# read and combine
combined_dt <- rbindlist(
  lapply(file_paths, function(fp) {
    dt <- fread(fp, sep = "\t")                     # tab-delimited read
    dt[, patient_id := basename(dirname(fp))]       # parent folder as ID
    dt                                              # return augmented dt
  }),
  use.names = TRUE,                                 # match columns by name
  fill      = TRUE                                  # in case some files differ
)

# (optional) move patient_id to front
setcolorder(
  combined_dt,
  c("patient_id", setdiff(names(combined_dt), "patient_id"))
)

# combined_dt now holds all your calls.txt with a patient_id column
