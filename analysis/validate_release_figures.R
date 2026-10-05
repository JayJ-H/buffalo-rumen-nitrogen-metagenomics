args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
Sys.setenv(TAHP_FIGURE_STYLE = "journal")
source(file.path(dirname(script), "validate_nitrogen_figures.R"))
