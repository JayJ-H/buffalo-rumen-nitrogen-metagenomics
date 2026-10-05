args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
local_library <- Sys.getenv("TAHP_R_LIBRARY")
if (nzchar(local_library)) .libPaths(c(local_library, .libPaths()))
input <- file.path(root, "results/nitrogen_function")
journal_style <- identical(Sys.getenv("TAHP_FIGURE_STYLE"), "journal")
src <- file.path(root, if (journal_style) "source_data/nitrogen_tahp_revised_20261004"
                 else "source_data/nitrogen_figures_20261004")
fig <- file.path(root, if (journal_style) "figures/nitrogen_tahp_revised"
                 else "figures/nitrogen_function")
read <- function(folder, name) read.csv(file.path(folder, name), check.names = FALSE)
validation <- jsonlite::fromJSON(file.path(fig, "figure_validation.json"))
for (name in names(validation$source_input_hashes)) {
  stopifnot(digest::digest(file = file.path(input, name), algo = "sha256") ==
              validation$source_input_hashes[[name]])
}
scores <- read(input, "scores_per_animal.csv")
long <- read(src, "Figure1a_all20_individual_abundance.csv")
stopifnot(nrow(long) == 600, length(unique(long$animal_id)) == 30,
          length(unique(long$Feature_ID)) == 20,
          all(table(long$Feature_ID, long$diet_group) == 6))
for (f in unique(long$Feature_ID)) {
  rows <- long[long$Feature_ID == f, ]
  ii <- match(rows$sample_id, scores$sample_id)
  stopifnot(max(abs(rows$percent_all_genes -
                    scores[[paste0(f, "__fraction")]][ii] * 100)) < 1e-12)
}
composition <- read(src, "Figure1c_KO_composition.csv")
stopifnot(nrow(composition) == 4, abs(sum(composition$Share) - 1) < 1e-12,
          abs(sum(composition$Share[composition$KO != "Other KOs"]) -
                .5619849140196044) < 1e-12)
carbon <- read(input, "carbon_function_associations.csv")
displayed <- read(src, "Figure2a_complete80_associations.csv")
stopifnot(nrow(displayed) == 80, identical(carbon$Pair_ID, displayed$Pair_ID),
          max(abs(carbon$Partial_r - displayed$Partial_r)) < 1e-12,
          all(displayed$q_HC3 > .05), sum(displayed$Permutation_support) == 24)
residuals <- read(src, "Figure2b_GS_starch_adjusted_residuals.csv")
target <- carbon[carbon$Pair_ID == "glutamine_synthetase__starch_alpha_glucan", ]
stopifnot(abs(cor(residuals$C_residual, residuals$N_residual) - target$Partial_r) < 1e-10,
          abs(coef(lm(N_residual ~ C_residual, residuals))[2] - target$Beta) < 1e-10)
volcano <- read(src, "FigureS2_all230_volcano_coordinates.csv")
original <- read(input, "diet_all_contrasts.csv")
stopifnot(nrow(volcano) == 230, all(table(volcano$Comparison) == 23),
          identical(volcano$Feature_ID, original$Feature_ID),
          max(abs(volcano$log2_ratio - original$log_ratio / log(2))) < 1e-12,
          max(abs(volcano$minus_log10_q + log10(original$q_all230))) < 1e-12)
heat <- read(src, "FigureS3_all600_profile_heatmap_values.csv")
stopifnot(nrow(heat) == 600,
          max(abs(tapply(heat$Z, heat$Feature_ID, mean))) < 1e-10,
          max(abs(tapply(heat$Z, heat$Feature_ID, sd) - 1)) < 1e-10)
pca_scores <- read(src, "FigureS1_PCA_scores_all30.csv")
pca_loadings <- read(src, "FigureS1_PCA_loadings_all20.csv")
pca_variance <- read(src, "FigureS1_PCA_variance_all20.csv")
stopifnot(nrow(pca_scores) == 30, nrow(pca_loadings) == 20,
          nrow(pca_variance) == 20,
          abs(sum(pca_variance$Variance_percent) - 100) < 1e-10)

tiff_metadata <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  bytes <- readBin(con, "raw", 2)
  endian <- if (identical(bytes, charToRaw("II"))) "little" else "big"
  stopifnot(readBin(con, "integer", 1, size = 2, signed = FALSE, endian = endian) == 42)
  offset <- readBin(con, "integer", 1, size = 4, endian = endian)
  seek(con, offset)
  count <- readBin(con, "integer", 1, size = 2, signed = FALSE, endian = endian)
  tags <- list()
  for (i in seq_len(count)) {
    tag <- readBin(con, "integer", 1, size = 2, signed = FALSE, endian = endian)
    type <- readBin(con, "integer", 1, size = 2, signed = FALSE, endian = endian)
    n <- readBin(con, "integer", 1, size = 4, endian = endian)
    raw <- readBin(con, "raw", 4)
    if (tag %in% c(256, 257, 282, 283, 296)) {
      rc <- rawConnection(raw)
      value <- readBin(rc, "integer", 1, size = if (type == 3) 2 else 4,
                       signed = type != 3, endian = endian)
      close(rc)
      tags[[as.character(tag)]] <- list(type = type, n = n, value = value)
    }
  }
  rational <- function(tag) {
    seek(con, tags[[as.character(tag)]]$value)
    value <- readBin(con, "integer", 2, size = 4, endian = endian)
    value[1] / value[2]
  }
  c(width = tags[["256"]]$value, height = tags[["257"]]$value,
    x_dpi = rational(282), y_dpi = rational(283),
    unit = tags[["296"]]$value)
}
exports <- read(fig, "export_validation.csv")
stopifnot(nrow(exports) == 6, all(exports$Width_mm == 174))
details <- list()
for (i in seq_len(nrow(exports))) {
  row <- exports[i, ]
  base <- file.path(fig, row$Figure)
  stopifnot(all(file.exists(paste0(base, c(".svg", ".pdf", ".tiff", ".png", "_preview.png")))))
  meta <- tiff_metadata(paste0(base, ".tiff"))
  stopifnot(meta["x_dpi"] == 600, meta["y_dpi"] == 600, meta["unit"] == 2,
            abs(meta["width"] - row$Width_mm / 25.4 * 600) <= 1,
            abs(meta["height"] - row$Height_mm / 25.4 * 600) <= 1)
  image <- png::readPNG(paste0(base, ".png"))
  stopifnot(abs(dim(image)[2] - row$Width_mm / 25.4 * 300) <= 1,
            abs(dim(image)[1] - row$Height_mm / 25.4 * 300) <= 1,
            sd(as.vector(image[, , 1:3])) > .03)
  con <- file(paste0(base, ".pdf"), "rb")
  signature <- rawToChar(readBin(con, "raw", 5))
  close(con)
  stopifnot(signature == "%PDF-")
  svg <- paste(readLines(paste0(base, ".svg"), warn = FALSE), collapse = "\n")
  stopifnot(grepl("<text", svg, fixed = TRUE),
            grepl("Arial", svg, fixed = TRUE),
            !grepl("NaN|Infinity", svg))
  details[[i]] <- c(Figure = row$Figure, meta)
}
jsonlite::write_json(list(
  status = "PASS", input_hashes_unchanged = TRUE,
  all_600_profile_values_verified = TRUE,
  all_230_volcano_values_verified = TRUE,
  all_80_heatmap_pairs_verified = TRUE,
  residuals_match_original_model = TRUE,
  PCA_counts_and_variance_verified = TRUE,
  all_6_TIFF_dimensions_and_600dpi_verified = TRUE,
  all_6_PNG_dimensions_and_nonblank_verified = TRUE,
  PDF_signatures_verified = TRUE, SVG_editable_text_verified = TRUE,
  export_metadata = details),
  file.path(fig, "independent_figure_validation.json"),
  pretty = TRUE, auto_unbox = TRUE)
cat("All numerical, dimensional and export checks PASS.\n")
