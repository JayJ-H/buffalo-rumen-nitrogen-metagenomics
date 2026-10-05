packages <- c("sandwich", "lmtest", "jsonlite", "digest", "ggplot2", "patchwork",
              "ggrepel", "svglite", "ragg", "png", "scales")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
cat("Dependencies available. See environment/ for the verified versions.\n")
