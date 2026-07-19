#!/usr/bin/env Rscript
# Purpose: Convert .rds files in data/ to .rda format so that data() can load them.
#
# Background:
#   - Bioconductor reviewer (mireia-bioinfo) requested migrating rds from
#     inst/extdata/ to data/ so users can load via data(object_name).
#   - R's data() does NOT support .rds directly; it requires .rda (or .RData).
#   - This script converts each data/<name>.rds to data/<name>.rda, ensuring
#     the saved object name matches the file basename (so that
#     data(<name>) loads an object called <name>).
#
# Run this BEFORE R CMD build, e.g.:
#   Rscript inst/scripts/convert_data_to_rda.R
#
# After running, data/ will contain only .rda files (the .rds versions are
# removed to keep the directory clean).

data_dir <- file.path("data")
if (!dir.exists(data_dir)) {
  stop("data/ directory not found. Run from package root.")
}

rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) {
  message("No .rds files in data/. Nothing to convert.")
  quit(status = 0)
}

message(sprintf("Found %d .rds files to convert:", length(rds_files)))

for (f in rds_files) {
  obj_name <- tools::file_path_sans_ext(basename(f))
  obj <- readRDS(f)
  # Ensure the object in the .rda is named after the file basename
  assign(obj_name, obj)
  rda_path <- file.path(data_dir, paste0(obj_name, ".rda"))
  save(list = obj_name, file = rda_path, compress = "xz", version = 2)
  message(sprintf("  %s -> %s", basename(f), basename(rda_path)))
  # Remove the original .rds to avoid confusion
  unlink(f)
}

# Final listing
message("\nFinal data/ contents:")
for (f in list.files(data_dir, full.names = TRUE)) {
  message(sprintf("  %s (%.1f KB)", basename(f), file.info(f)$size / 1024))
}
