args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/")
local_library <- Sys.getenv("TAHP_R_LIBRARY")
if (nzchar(local_library)) .libPaths(c(local_library, .libPaths()))
packages <- c("sandwich", "lmtest", "jsonlite", "digest", "ggplot2", "patchwork",
              "ggrepel", "svglite", "ragg", "png", "scales")
versions <- data.frame(Package = c("R", packages),
                       Version = c(as.character(getRversion()),
                                   vapply(packages, function(p) as.character(packageVersion(p)),
                                          character(1))))
dir.create(file.path(root, "environment"), showWarnings = FALSE)
write.csv(versions, file.path(root, "environment/R_versions.csv"), row.names = FALSE)
