args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
local_library <- Sys.getenv("TAHP_R_LIBRARY")
if (nzchar(local_library)) .libPaths(c(local_library, .libPaths()))
needed <- c("ggplot2", "patchwork", "ggrepel", "svglite", "ragg",
            "png", "jsonlite", "digest", "scales")
stopifnot(all(vapply(needed, requireNamespace, logical(1), quietly = TRUE)))
library(ggplot2)
library(patchwork)
options(stringsAsFactors = FALSE, digits = 7)
set.seed(20261004)

input <- file.path(root, "results/nitrogen_function")
journal_style <- identical(Sys.getenv("TAHP_FIGURE_STYLE"), "journal")
out <- file.path(root, if (journal_style) "figures/nitrogen_tahp_revised" else
                   "figures/nitrogen_function")
source_out <- file.path(root, if (journal_style)
  "source_data/nitrogen_tahp_revised_20261004" else "source_data/nitrogen_figures_20261004")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(source_out, recursive = TRUE, showWarnings = FALSE)
input_paths <- list.files(input, full.names = TRUE)
input_hashes <- setNames(vapply(input_paths, function(p) {
  digest::digest(file = p, algo = "sha256")
}, character(1)), basename(input_paths))
read_input <- function(name) {
  read.csv(file.path(input, name), check.names = FALSE,
           colClasses = if (name %in% c("scores_per_animal.csv",
             "carbon_pair_scores.csv")) c(animal_id = "character") else NA)
}
save_source <- function(data, name) {
  previous <- options(digits = 17)
  on.exit(options(previous))
  write.csv(data, file.path(source_out, name), row.names = FALSE, na = "")
}
compact_number <- function(x) trimws(formatC(x, format = "g", digits = 3))
cfg <- jsonlite::fromJSON(file.path(input, "definition_freeze.json"),
                         simplifyVector = FALSE)
features <- names(cfg$base_sets)
stopifnot(length(features) == 20)
scores <- read_input("scores_per_animal.csv")
desc <- read_input("functional_description.csv")
ko <- read_input("KO_contribution.csv")
diet <- read_input("diet_omnibus.csv")
means <- read_input("diet_adjusted_means.csv")
carbon <- read_input("carbon_function_associations.csv")
pair_scores <- read_input("carbon_pair_scores.csv")
loo <- read_input("carbon_LOO.csv")
contrasts <- read_input("diet_all_contrasts.csv")
stopifnot(nrow(scores) == 30, nrow(diet) == 23,
          nrow(carbon) == 80, nrow(contrasts) == 230,
          all(table(scores$diet_group) == 6),
          all(diet$q_HC3 > .05), all(contrasts$q_all230 > .05),
          all(carbon$q_HC3 > .05), sum(carbon$q_permutation < .05) == 24)

labels <- c(
  Arginine_biosynthesis = "Arginine biosynthesis",
  "Alanine,_aspartate_and_glutamate_metabolism" = "Ala / Asp / Glu metabolism",
  "Glycine,_serine_and_threonine_metabolism" = "Gly / Ser / Thr metabolism",
  Cysteine_and_methionine_metabolism = "Cys / Met metabolism",
  "Valine,_leucine_and_isoleucine_degradation" = "BCAA degradation",
  "Valine,_leucine_and_isoleucine_biosynthesis" = "BCAA biosynthesis",
  Lysine_biosynthesis = "Lysine biosynthesis",
  Lysine_degradation = "Lysine degradation",
  Arginine_and_proline_metabolism = "Arg / Pro metabolism",
  Histidine_metabolism = "Histidine metabolism",
  Tyrosine_metabolism = "Tyrosine metabolism",
  Phenylalanine_metabolism = "Phenylalanine metabolism",
  Tryptophan_metabolism = "Tryptophan metabolism",
  "Phenylalanine,_tyrosine_and_tryptophan_biosynthesis" = "Phe / Tyr / Trp biosynthesis",
  Nitrogen_metabolism = "Nitrogen metabolism",
  Biosynthesis_of_amino_acids = "Amino-acid biosynthesis",
  glutamine_synthetase = "Glutamine synthetase (GS)",
  glutamate_synthase = "Glutamate synthase (GOGAT)",
  glutamate_dehydrogenase = "Glutamate dehydrogenase (GDH)",
  urease_structural_genes = "Urease structural genes",
  BCAA_synthesis_to_degradation = "BCAA synthesis / degradation",
  lysine_synthesis_to_degradation = "Lysine synthesis / degradation",
  GS_to_GDH = "GS / GDH"
)
short <- c(Nitrogen_metabolism = "Nitrogen metabolism",
           glutamine_synthetase = "GS", glutamate_synthase = "GOGAT",
           glutamate_dehydrogenase = "GDH", urease_structural_genes = "Urease")
enzymes <- c("glutamine_synthetase", "glutamate_synthase",
             "glutamate_dehydrogenase", "urease_structural_genes")
focus <- c("Nitrogen_metabolism", "glutamine_synthetase",
           "glutamate_synthase", "glutamate_dehydrogenase")
diet_colors <- c(A = "#477D9F", B = "#DF9A63", C = "#63A7A0",
                 D = "#C77F95", E = "#9CAF68")
enzyme_colors <- c(glutamine_synthetase = "#477D9F",
                   glutamate_synthase = "#63A7A0",
                   glutamate_dehydrogenase = "#C77F95",
                   urease_structural_genes = "#DF9A63")
focus_colors <- c(Nitrogen_metabolism = "#687983", enzyme_colors[1:3])
class_colors <- c(Pathway = "#78AFA9", Enzyme = "#BA819D", Ratio = "#DEAD74")
theme_set(theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.text = element_text(size = 6.2, colour = "#34434B"),
        axis.title = element_text(size = 6.8),
        axis.line = element_line(linewidth = .3, colour = "#56646B"),
        axis.ticks = element_line(linewidth = .3, colour = "#56646B"),
        axis.ticks.length = grid::unit(1.2, "mm"),
        plot.title = element_text(size = 7.8, face = "bold", colour = "#25343C"),
        plot.subtitle = element_text(size = 6.3, colour = "#657179"),
        plot.caption = element_text(size = 6.1, hjust = 0, colour = "#657179"),
        plot.tag = element_text(size = 10, face = "bold", colour = "#25343C"),
        plot.margin = margin(5, 6, 5, 5),
        strip.text = element_text(size = 6.6, face = "bold"),
        strip.background = element_rect(fill = "#F0F3F3", colour = NA),
        legend.key.size = grid::unit(3, "mm"),
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 6.2),
        legend.background = element_blank(),
        panel.grid = element_blank()))
median_iqr <- function(x) {
  data.frame(y = median(x), ymin = quantile(x, .25), ymax = quantile(x, .75))
}
feature_long <- do.call(rbind, lapply(features, function(f) {
  data.frame(sample_id = scores$sample_id, animal_id = scores$animal_id,
             diet_group = scores$diet_group, Feature_ID = f,
             log_fraction = scores[[f]], percent_all_genes =
               scores[[paste0(f, "__fraction")]] * 100)
}))
feature_long$diet_group <- factor(feature_long$diet_group, levels = LETTERS[1:5])
feature_long$Family <- ifelse(feature_long$Feature_ID %in% enzymes, "Enzyme", "Pathway")
stopifnot(nrow(feature_long) == 600,
          max(abs(exp(feature_long$log_fraction) * 100 -
                    feature_long$percent_all_genes)) < 1e-12)
order_mean <- names(sort(tapply(feature_long$percent_all_genes,
                               feature_long$Feature_ID, mean)))
feature_long$Feature <- factor(labels[feature_long$Feature_ID],
                               levels = unname(labels[order_mean]))
save_source(feature_long, "Figure1a_all20_individual_abundance.csv")
enzyme_long <- feature_long[feature_long$Feature_ID %in% enzymes, ]
enzyme_long$Enzyme <- factor(short[enzyme_long$Feature_ID],
                             levels = unname(short[enzymes]))
save_source(enzyme_long, "Figure1b_enzyme_diet_individuals.csv")

p1a <- ggplot(feature_long, aes(x = percent_all_genes, y = Feature)) +
  geom_violin(aes(fill = Family), orientation = "y", scale = "width",
              trim = TRUE, alpha = .48, linewidth = .25, colour = "#849399") +
  geom_point(aes(colour = Family),
             position = position_jitter(height = .12, width = 0, seed = 20261004),
             size = .6, alpha = .42) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 1.35,
               fill = "white", colour = "#263940", stroke = .3) +
  scale_x_log10(breaks = c(.00001, .0001, .001, .01, .1, 1),
                labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1", "1"),
                expand = expansion(mult = c(.035, .1))) +
  scale_fill_manual(values = class_colors) +
  scale_colour_manual(values = class_colors) +
  labs(title = "Functional gene-abundance distributions",
       subtitle = "20 fixed collections; 30 individual buffaloes",
       x = "Relative abundance (% of all genes; log scale)", y = NULL) +
  theme(legend.position = "none", axis.text.y = element_text(size = 6.2))

p1b <- ggplot(enzyme_long, aes(diet_group, percent_all_genes,
                               fill = diet_group, colour = diet_group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = .35, linewidth = .25) +
  geom_boxplot(width = .16, outlier.shape = NA, fill = "white",
               colour = "#435159", linewidth = .3) +
  geom_point(position = position_jitter(width = .1, seed = 20261004),
             size = 1, alpha = .8) +
  facet_wrap(~ Enzyme, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = diet_colors) +
  scale_colour_manual(values = diet_colors) +
  scale_y_continuous(labels = compact_number,
                     expand = expansion(mult = c(.1, .15))) +
  labs(title = "Enzyme-defined\ncollections", subtitle = "n = 6 per diet",
       x = "Diet", y = "Relative abundance (%)") +
  theme(legend.position = "none", panel.spacing = grid::unit(2.6, "mm"),
        axis.text = element_text(size = 5.6), axis.title.y = element_text(size = 6.1))

ko_n <- ko[ko$Feature_ID == "Nitrogen_metabolism", ]
top_ko <- c("K00262", "K01915", "K00266")
shares <- ko_n$Cohort_set_mass_share[match(top_ko, ko_n$KO)]
stopifnot(abs(sum(shares) - .5619849140196044) < 1e-12)
composition <- data.frame(KO = c(top_ko, "Other KOs"),
                          Share = c(shares, 1 - sum(shares)))
composition$KO <- factor(composition$KO, levels = composition$KO)
composition$ymax <- cumsum(composition$Share)
composition$ymin <- c(0, head(composition$ymax, -1))
composition$Label <- sprintf("%s  %.2f%%", composition$KO, composition$Share * 100)
ko_colors <- c(K00262 = "#C77F95", K01915 = "#477D9F",
               K00266 = "#63A7A0", "Other KOs" = "#DDE4E5")
save_source(composition, "Figure1c_KO_composition.csv")
save_source(ko_n, "Figure1c_all_nitrogen_KO_contributions.csv")
p1c <- ggplot(composition) +
  geom_rect(aes(ymin = ymin, ymax = ymax, xmin = 1.08, xmax = 2, fill = KO),
            colour = "white", linewidth = .55) +
  annotate("text", x = 0, y = 0, label = "56.20%\nThree major KOs",
           family = "Arial", size = 3, colour = "#263940", lineheight = 1.15) +
  coord_polar(theta = "y", start = -.2) +
  scale_x_continuous(limits = c(0, 2.15)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = ko_colors, labels = composition$Label) +
  labs(title = "Major nitrogen-metabolism KOs", fill = NULL) +
  theme_void(base_size = 7, base_family = "Arial") +
  theme(plot.title = element_text(size = 7.8, face = "bold", colour = "#25343C"),
        plot.margin = margin(5, 6, 5, 5), legend.position = "right",
        legend.text = element_text(size = 6.4),
        legend.key.size = grid::unit(3.1, "mm"),
        plot.tag = element_text(size = 10, face = "bold"))

nitrogen <- feature_long[feature_long$Feature_ID == "Nitrogen_metabolism", ]
nitrogen_test <- diet[diet$Feature_ID == "Nitrogen_metabolism", ]
nitrogen_means <- means[means$Feature_ID == "Nitrogen_metabolism", ]
save_source(nitrogen, "Figure1d_nitrogen_individuals.csv")
save_source(nitrogen_means, "Figure1d_nitrogen_adjusted_means.csv")
p1d <- ggplot(nitrogen, aes(diet_group, percent_all_genes,
                           fill = diet_group, colour = diet_group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = .35, linewidth = .3) +
  geom_boxplot(width = .18, outlier.shape = NA, fill = "white",
               colour = "#435159", linewidth = .35) +
  geom_point(position = position_jitter(width = .1, seed = 20261004), size = 1.4) +
  scale_fill_manual(values = diet_colors) + scale_colour_manual(values = diet_colors) +
  scale_y_continuous(labels = compact_number, expand = expansion(mult = c(.1, .1))) +
  labs(title = "Nitrogen-metabolism\nprofiles",
       subtitle = sprintf("HC3 P = %.4f; BH q = %.4f",
                          nitrogen_test$P_HC3, nitrogen_test$q_HC3),
       x = "Diet", y = "Relative abundance (%)") +
  theme(legend.position = "none")

figure1 <- p1a + p1b + p1c + p1d +
  plot_layout(design = "AAB\nAAB\nCCD", heights = c(1, 1, 1.35)) +
  plot_annotation(tag_levels = "a",
    caption = "Each point is one buffalo. Abundances use the all-gene denominator; KO shares use cohort collection mass.")

carbon_labels <- c(starch_alpha_glucan = "Starch /\nalpha-glucan",
                   cellulose_chain = "Cellulose", xylan = "Xylan",
                   pectin_backbone = "Pectin")
carbon$Feature <- factor(labels[carbon$Nitrogen_feature],
                         levels = rev(unname(labels[features])))
carbon$Carbon <- factor(carbon$Carbon_feature, levels = names(carbon_labels))
carbon$Permutation_support <- carbon$q_permutation < .05
save_source(carbon, "Figure2a_complete80_associations.csv")
p2a <- ggplot(carbon, aes(Carbon, Feature, fill = Partial_r)) +
  geom_tile(colour = "white", linewidth = .45) +
  geom_text(aes(label = sprintf("%.2f", Partial_r)), size = 2.05, colour = "#243A43") +
  geom_point(data = carbon[carbon$Permutation_support, ],
             aes(x = as.numeric(Carbon) + .28, y = as.numeric(Feature) + .27),
             inherit.aes = FALSE, shape = 21, size = .65, stroke = .18,
             colour = "#26383F", fill = "white") +
  scale_fill_gradient2(low = "#6DADA9", mid = "#F8FAFA", high = "#CE8299",
                       midpoint = 0, limits = c(-1, 1), name = "Partial r",
                       breaks = c(-1, -.5, 0, .5, 1)) +
  scale_x_discrete(labels = carbon_labels, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Carbon-nitrogen conditional covariation",
       subtitle = "80 fixed pairs; circles: permutation q < 0.05",
       x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 5.8),
        axis.text.y = element_text(size = 6),
        legend.position = "bottom",
        legend.key.width = grid::unit(5.5, "mm"),
        legend.key.height = grid::unit(2.5, "mm"))

gs_pair <- pair_scores[pair_scores$Pair_ID ==
                        "glutamine_synthetase__starch_alpha_glucan", ]
Z <- model.matrix(~ factor(diet_group, levels = LETTERS[1:5]) +
                    initial_BW_kg + log_KO_annotation_load, gs_pair)
gs_pair$N_residual <- lm.fit(Z, gs_pair$response)$residuals
gs_pair$C_residual <- lm.fit(Z, gs_pair$carbon)$residuals
gs_row <- carbon[carbon$Pair_ID == unique(gs_pair$Pair_ID), ]
residual_beta <- coef(lm(N_residual ~ C_residual, gs_pair))[2]
stopifnot(nrow(gs_pair) == 30,
          abs(cor(gs_pair$N_residual, gs_pair$C_residual) - gs_row$Partial_r) < 1e-10,
          abs(residual_beta - gs_row$Beta) < 1e-10)
save_source(gs_pair, "Figure2b_GS_starch_adjusted_residuals.csv")
p2b <- ggplot(gs_pair, aes(C_residual, N_residual)) +
  geom_hline(yintercept = 0, linewidth = .25, colour = "#DCE2E3") +
  geom_vline(xintercept = 0, linewidth = .25, colour = "#DCE2E3") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linewidth = .65, colour = "#465A65") +
  geom_point(aes(fill = diet_group), shape = 21, size = 2.1,
             colour = "white", stroke = .3) +
  scale_fill_manual(values = diet_colors) +
  scale_x_continuous(labels = compact_number) +
  scale_y_continuous(labels = compact_number) +
  labs(title = "GS and starch-processing profiles",
       subtitle = sprintf("r = %.3f; HC3 q = %.4f\nPermutation q = %.4f",
                          gs_row$Partial_r, gs_row$q_HC3, gs_row$q_permutation),
       x = "Adjusted starch / alpha-glucan log abundance",
       y = "Adjusted GS log abundance", fill = "Diet") +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title = element_text(size = 6.3),
        legend.margin = margin(0, 0, 0, 0))

forest <- carbon[carbon$Carbon_feature == "starch_alpha_glucan" &
                   carbon$Nitrogen_feature %in% focus, ]
forest <- forest[match(focus, forest$Nitrogen_feature), ]
forest$Feature_short <- factor(short[forest$Nitrogen_feature],
                               levels = rev(unname(short[focus])))
save_source(forest, "Figure2c_four_prespecified_HC3_effects.csv")
p2c <- ggplot(forest, aes(Beta, Feature_short, colour = Nitrogen_feature)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = .35, colour = "#A3AFB3") +
  geom_segment(aes(x = CI_lower, xend = CI_upper,
                   y = Feature_short, yend = Feature_short), linewidth = .7) +
  geom_point(size = 2.3) +
  geom_text(aes(label = sprintf("q = %.3f", q_HC3)), x = .91,
            hjust = 1, size = 2.05, colour = "#647179") +
  scale_colour_manual(values = focus_colors) +
  scale_x_continuous(limits = c(min(forest$CI_lower) - .05, .95),
                     breaks = c(-1, -.5, 0), labels = compact_number, expand = c(0, 0)) +
  labs(title = "Adjusted association estimates",
       subtitle = "HC3 95% confidence intervals\nBH across 80 pairs",
       x = "Log-log regression coefficient", y = NULL) +
  theme(legend.position = "none", axis.text.y = element_text(size = 6))

loo_focus <- loo[loo$Pair_ID %in% forest$Pair_ID, ]
loo_focus$Nitrogen_feature <- sub("__starch_alpha_glucan$", "", loo_focus$Pair_ID)
loo_focus$Feature_short <- factor(short[loo_focus$Nitrogen_feature],
                                  levels = levels(forest$Feature_short))
stopifnot(nrow(loo_focus) == 120, all(loo_focus$Beta < 0))
save_source(loo_focus, "Figure2d_all120_leave_one_out_effects.csv")
p2d <- ggplot(loo_focus, aes(Beta, Feature_short, fill = Nitrogen_feature)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = .3, colour = "#A3AFB3") +
  geom_violin(orientation = "y", scale = "width", trim = TRUE,
              alpha = .35, linewidth = .25, colour = "#819198") +
  geom_point(position = position_jitter(width = 0, height = .1, seed = 20261004),
             size = .7, colour = "#63747D", alpha = .65) +
  geom_point(data = forest, aes(Beta, Feature_short), inherit.aes = FALSE,
             shape = 23, size = 2.1, fill = "white", colour = "#263940", stroke = .4) +
  scale_fill_manual(values = focus_colors) +
  scale_x_continuous(labels = compact_number) +
  labs(title = "Leave-one-buffalo-out stability",
       subtitle = "30 omissions per pair\nDiamonds: full-cohort estimates",
       x = "Log-log regression coefficient", y = NULL) +
  theme(legend.position = "none", axis.text.y = element_text(size = 6))

figure2 <- p2a + p2b + p2c + p2d +
  plot_layout(design = "AB\nAC\nAD", widths = c(1.15, 1), heights = c(1.1, .8, .8)) +
  plot_annotation(tag_levels = "a",
    caption = "Models adjust for diet, initial body weight and KO annotation load. HC3 is primary; permutation is sensitivity.")

matrix_log <- as.matrix(scores[, features])
rownames(matrix_log) <- scores$sample_id
pca <- prcomp(matrix_log, center = TRUE, scale. = TRUE)
variance <- pca$sdev^2 / sum(pca$sdev^2)
pca_scores <- data.frame(sample_id = scores$sample_id, animal_id = scores$animal_id,
                         diet_group = scores$diet_group,
                         PC1 = pca$x[, 1], PC2 = pca$x[, 2])
centroids <- aggregate(cbind(PC1, PC2) ~ diet_group, pca_scores, mean)
pca_segments <- merge(pca_scores, centroids, by = "diet_group",
                       suffixes = c("", "_centroid"))
scree <- data.frame(PC = seq_along(variance), Variance_percent = variance * 100,
                     Cumulative_percent = cumsum(variance) * 100)
loadings <- data.frame(Feature_ID = features, Label = labels[features],
                       PC1 = pca$rotation[, 1], PC2 = pca$rotation[, 2])
save_source(pca_scores, "FigureS1_PCA_scores_all30.csv")
save_source(centroids, "FigureS1_PCA_diet_centroids.csv")
save_source(scree, "FigureS1_PCA_variance_all20.csv")
save_source(loadings, "FigureS1_PCA_loadings_all20.csv")
stopifnot(abs(sum(variance) - 1) < 1e-12,
          max(abs(cov(pca$x)[row(cov(pca$x)) != col(cov(pca$x))])) < 1e-9)
ps1a <- ggplot(pca_scores, aes(PC1, PC2, fill = diet_group)) +
  geom_hline(yintercept = 0, colour = "#E1E6E7", linewidth = .3) +
  geom_vline(xintercept = 0, colour = "#E1E6E7", linewidth = .3) +
  geom_segment(data = pca_segments,
               aes(xend = PC1_centroid, yend = PC2_centroid, colour = diet_group),
               linewidth = .3, alpha = .4) +
  geom_point(aes(shape = diet_group), size = 2.7, colour = "white", stroke = .45) +
  geom_point(data = centroids, size = 4.2, shape = 23, colour = "#40545E", stroke = .5) +
  scale_fill_manual(values = diet_colors) + scale_colour_manual(values = diet_colors) +
  scale_shape_manual(values = c(A = 21, B = 22, C = 24, D = 23, E = 25)) +
  scale_x_continuous(labels = compact_number) +
  scale_y_continuous(labels = compact_number) +
  labs(title = "PCA of nitrogen-related functional profiles",
       subtitle = "20 centered and scaled log-abundance variables",
       x = sprintf("PC1 (%.2f%%)", variance[1] * 100),
       y = sprintf("PC2 (%.2f%%)", variance[2] * 100),
       fill = "Diet", colour = "Diet", shape = "Diet") +
  theme(legend.position = "bottom")
ps1b <- ggplot(scree[1:8, ], aes(factor(PC), Variance_percent)) +
  geom_col(fill = "#78AFA9", width = .68) +
  geom_text(aes(label = sprintf("%.1f", Variance_percent)),
            vjust = -.4, size = 2) +
  scale_y_continuous(labels = compact_number, expand = expansion(mult = c(0, .15))) +
  labs(title = "Variance explained", x = "Principal component", y = "Variance (%)") +
  theme(axis.text = element_text(size = 5.7))
top_load <- loadings[order(abs(loadings$PC1), decreasing = TRUE)[1:12], ]
top_load$Label <- gsub(" metabolism", "", top_load$Label, fixed = TRUE)
top_load$Label[top_load$Feature_ID %in% enzymes] <-
  short[top_load$Feature_ID[top_load$Feature_ID %in% enzymes]]
top_load$Label <- factor(top_load$Label, levels = rev(top_load$Label))
ps1c <- ggplot(top_load, aes(PC1, Label)) +
  geom_segment(aes(x = 0, xend = PC1, yend = Label), linewidth = .65, colour = "#B5C9CA") +
  geom_point(size = 1.9, colour = "#BA819D") +
  geom_vline(xintercept = 0, linewidth = .25, colour = "#D2DDDF") +
  scale_x_continuous(labels = compact_number) +
  labs(title = "Largest PC1 loadings", x = "Loading", y = NULL) +
  theme(axis.text.y = element_text(size = 5.5))
figure_s1 <- ps1a + ps1b + ps1c +
  plot_layout(design = "AB\nAC", widths = c(1.35, 1), heights = c(.7, 1.3)) +
  plot_annotation(tag_levels = "a",
    caption = "All 30 animals are shown. Diamonds denote diet centroids; connecting lines are descriptive, not confidence intervals.")

volcano <- contrasts
volcano$log2_ratio <- volcano$log_ratio / log(2)
volcano$minus_log10_q <- -log10(volcano$q_all230)
volcano$Comparison <- paste(volcano$comparison, "vs", volcano$reference)
volcano$Family <- ifelse(volcano$Feature_ID %in% enzymes, "Enzyme",
                         ifelse(volcano$Feature_ID %in% names(cfg$ratios), "Ratio", "Pathway"))
volcano$Short_label <- ifelse(volcano$Feature_ID %in% enzymes[1:3],
                              short[volcano$Feature_ID], "")
save_source(volcano, "FigureS2_all230_volcano_coordinates.csv")
ps2 <- ggplot(volcano, aes(log2_ratio, minus_log10_q, colour = Family)) +
  geom_hline(yintercept = -log10(.05), linetype = "dashed",
              colour = "#6D7E86", linewidth = .35) +
  geom_vline(xintercept = 0, colour = "#DDE4E5", linewidth = .25) +
  geom_point(aes(shape = Family), size = 1.85, alpha = .8) +
  ggrepel::geom_text_repel(data = volcano[volcano$Short_label != "", ],
                           aes(label = Short_label), size = 1.95,
                           seed = 20261004, max.overlaps = Inf,
                           min.segment.length = 0, segment.size = .18,
                           box.padding = .1, point.padding = .1) +
  facet_wrap(~ Comparison, ncol = 2) +
  scale_colour_manual(values = class_colors) +
  scale_shape_manual(values = c(Pathway = 16, Enzyme = 17, Ratio = 15)) +
  scale_x_continuous(labels = compact_number) +
  scale_y_continuous(limits = c(0, 1.48), breaks = c(0, .5, 1, 1.3), labels = compact_number,
                     expand = c(0, 0)) +
  labs(title = "Complete diet-contrast volcano atlas",
       subtitle = "All 23 responses in each contrast; horizontal line: BH q = 0.05",
       x = "Adjusted abundance ratio (log2 scale)",
       y = "-log10(BH q across all 230 contrasts)",
       colour = "Collection", shape = "Collection") +
  theme(legend.position = "bottom", panel.spacing = grid::unit(3.2, "mm"))

standardized <- scale(matrix_log)
heat <- do.call(rbind, lapply(seq_along(features), function(j) {
  data.frame(sample_id = scores$sample_id, animal_id = scores$animal_id,
             diet_group = scores$diet_group, Feature_ID = features[j],
             Log_abundance = matrix_log[, j], Z = standardized[, j])
}))
animal_order <- scores$animal_id[order(scores$diet_group, as.integer(scores$animal_id))]
heat$Animal <- factor(heat$animal_id, levels = animal_order)
heat$Feature <- factor(labels[heat$Feature_ID], levels = rev(unname(labels[features])))
heat$diet_group <- factor(heat$diet_group, levels = LETTERS[1:5])
save_source(heat, "FigureS3_all600_profile_heatmap_values.csv")
stopifnot(max(abs(tapply(heat$Z, heat$Feature_ID, mean))) < 1e-10,
          max(abs(tapply(heat$Z, heat$Feature_ID, sd) - 1)) < 1e-10)
ps3 <- ggplot(heat, aes(Animal, Feature, fill = Z)) +
  geom_tile(colour = "white", linewidth = .15) +
  facet_grid(~ diet_group, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#5D9F9B", mid = "#F8FAFA", high = "#C77590",
                       midpoint = 0, name = "Within-feature\nz score") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Individual nitrogen-related functional profiles",
       subtitle = "20 fixed collections x 30 animals; standardized log abundances",
       x = "Individual buffalo", y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 5.8),
        legend.position = "right", panel.spacing.x = grid::unit(1.6, "mm"))

atlas <- feature_long
atlas$Feature <- factor(labels[atlas$Feature_ID], levels = unname(labels[features]))
atlas$Facet <- factor(vapply(as.character(atlas$Feature), function(label) {
  f <- names(labels)[match(label, labels)]
  q <- diet$q_HC3[match(f, diet$Feature_ID)]
  paste0(paste(strwrap(label, width = 23), collapse = "\n"),
         sprintf("\nBH q = %.3f", q))
}, character(1)), levels = unique(vapply(unname(labels[features]), function(label) {
  f <- names(labels)[match(label, labels)]
  q <- diet$q_HC3[match(f, diet$Feature_ID)]
  paste0(paste(strwrap(label, width = 23), collapse = "\n"),
         sprintf("\nBH q = %.3f", q))
}, character(1))))
save_source(atlas, "FigureS4_all600_diet_violin_values.csv")
ps4 <- ggplot(atlas, aes(diet_group, percent_all_genes,
                         fill = diet_group, colour = diet_group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = .35, linewidth = .25) +
  geom_boxplot(width = .16, outlier.shape = NA, fill = "white",
               colour = "#435159", linewidth = .25) +
  geom_point(position = position_jitter(width = .1, seed = 20261004),
             size = .9, alpha = .8) +
  facet_wrap(~ Facet, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = diet_colors) + scale_colour_manual(values = diet_colors) +
  scale_y_continuous(labels = compact_number, expand = expansion(mult = c(.1, .1))) +
  labs(title = "Complete diet-by-function violin atlas",
       subtitle = "All 20 base collections; HC3 omnibus q values use the complete 23-response family",
       x = "Diet", y = "Relative abundance (% of all genes)") +
  theme(legend.position = "none", strip.text = element_text(size = 5.8),
        axis.text = element_text(size = 5.1),
        panel.spacing = grid::unit(2.2, "mm"))

if (journal_style) {
  source(file.path(root, "analysis/style_nitrogen_figures_tahp.R"), local = TRUE)
}
exports <- list()
save_figure <- function(plot, id, width_mm, height_mm, panels) {
  base <- file.path(out, id)
  svglite::svglite(paste0(base, ".svg"), width = width_mm / 25.4,
                   height = height_mm / 25.4, bg = "white")
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_mm / 25.4,
                       height = height_mm / 25.4, family = "Arial", bg = "white")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(base, ".tiff"), width = width_mm, height = height_mm,
                 units = "mm", res = 600, compression = "lzw", background = "white")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(base, ".png"), width = width_mm, height = height_mm,
                units = "mm", res = 300, background = "white")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(base, "_preview.png"), width = width_mm, height = height_mm,
                units = "mm", res = 240, background = "white")
  print(plot)
  dev.off()
  image <- png::readPNG(paste0(base, "_preview.png"))
  rgb <- image[, , 1:3]
  border <- c(rgb[1:3, , ], rgb[(dim(rgb)[1] - 2):dim(rgb)[1], , ],
              rgb[, 1:3, ], rgb[, (dim(rgb)[2] - 2):dim(rgb)[2], ])
  svg <- paste(readLines(paste0(base, ".svg"), warn = FALSE), collapse = "\n")
  text_count <- length(gregexpr("<text", svg, fixed = TRUE)[[1]])
  stopifnot(sd(as.vector(rgb)) > .03, mean(border < .9) < .005,
            text_count > 5, grepl("Arial", svg, fixed = TRUE))
  exports[[length(exports) + 1]] <<- data.frame(
    Figure = id, Width_mm = width_mm, Height_mm = height_mm,
    TIFF_dpi = 600, PNG_dpi = 300, Preview_dpi = 240,
    Preview_width_px = dim(image)[2], Preview_height_px = dim(image)[1],
    Pixel_SD = sd(as.vector(rgb)), Edge_ink_fraction = mean(border < .9),
    SVG_text_elements = text_count, Panels = panels)
}

save_figure(figure1, "Figure1_nitrogen_composition", 174, if (journal_style) 195 else 185, 4)
save_figure(figure2, "Figure2_carbon_nitrogen_associations", 174, if (journal_style) 200 else 190, 4)
save_figure(figure_s1, "FigureS1_functional_PCA", 174, if (journal_style) 155 else 125, 3)
save_figure(ps2, "FigureS2_complete_volcano_atlas", 174, 205, 10)
save_figure(ps3, "FigureS3_individual_function_heatmap", 174, 145, 1)
save_figure(ps4, "FigureS4_complete_diet_violin_atlas", 174, 225, 20)
export_table <- do.call(rbind, exports)
write.csv(export_table, file.path(out, "export_validation.csv"), row.names = FALSE)

images <- lapply(export_table$Figure, function(id) {
  png::readPNG(file.path(out, paste0(id, "_preview.png")))
})
ragg::agg_png(file.path(out, "contact_sheet.png"), width = 2700, height = 3000,
              units = "px", res = 220, background = "white")
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 3)))
for (i in seq_along(images)) {
  grid::pushViewport(grid::viewport(layout.pos.row = (i - 1) %/% 3 + 1,
                                    layout.pos.col = (i - 1) %% 3 + 1))
  img <- images[[i]]
  fitted <- min(900 * .93 / ncol(img), 1500 * .88 / nrow(img))
  grid::grid.raster(img, x = .5, y = .48,
                    width = grid::unit(ncol(img) * fitted / 900, "npc"),
                    height = grid::unit(nrow(img) * fitted / 1500, "npc"))
  grid::grid.text(sub("_.*$", "", export_table$Figure[i]), x = .04, y = .985,
                  just = c("left", "top"),
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 12, fontface = "bold"))
  grid::popViewport()
}
grid::popViewport()
dev.off()
stopifnot(identical(input_hashes,
                   setNames(vapply(input_paths, function(p) {
                     digest::digest(file = p, algo = "sha256")
                   }, character(1)), basename(input_paths))))
jsonlite::write_json(list(
  status = "PASS", backend = "R", figure_pages = 6,
  style = if (journal_style) "TAHP revised" else "original",
  source_input_hashes = input_hashes,
  historical_statistics_unchanged = TRUE,
  all_30_animals_included = TRUE, base_feature_count = 20,
  full_diet_contrasts = 230, full_carbon_pairs = 80,
  primary_HC3_FDR_discoveries = 0, permutation_sensitivity_discoveries = 24,
  KO_top3_mass_share = sum(shares),
  residual_GS_r = cor(gs_pair$N_residual, gs_pair$C_residual),
  residual_GS_beta = unname(residual_beta),
  PCA_transform = "centered and unit-SD-scaled natural-log all-gene fractions",
  exports = export_table),
  file.path(out, "figure_validation.json"), auto_unbox = TRUE,
  pretty = TRUE, digits = 16)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Six R figure pages and complete panel source tables exported.\n")
