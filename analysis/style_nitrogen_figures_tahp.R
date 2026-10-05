# Presentation-only overrides. All analysis objects come from the frozen renderer.
ink <- "#222222"
diet_colors <- c(A = "#4477AA", B = "#EE7733", C = "#228833",
                 D = "#CC6677", E = "#AA3377")
diet_shapes <- c(A = 21, B = 22, C = 24, D = 23, E = 25)
class_colors <- c(Pathway = "#4477AA", Enzyme = "#CC6677", Ratio = "#777777")
enzyme_colors <- c(glutamine_synthetase = "#4477AA",
                   glutamate_synthase = "#228833",
                   glutamate_dehydrogenase = "#CC6677",
                   urease_structural_genes = "#EE7733")
focus_colors <- c(Nitrogen_metabolism = "#555555", enzyme_colors[1:3])
base_theme <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(text = element_text(colour = ink),
        axis.text = element_text(size = 8, colour = ink),
        axis.title = element_text(size = 8.5, colour = ink),
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.title.y = element_text(margin = margin(r = 4)),
        axis.line = element_line(linewidth = .35, colour = ink),
        axis.ticks = element_line(linewidth = .3, colour = ink),
        axis.ticks.length = grid::unit(1, "mm"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold", colour = ink),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = grid::unit(3.3, "mm"),
        legend.background = element_blank(),
        plot.title = element_blank(), plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        plot.tag = element_text(size = 10, face = "bold", colour = ink),
        plot.margin = margin(5, 5, 5, 5),
        panel.spacing = grid::unit(3, "mm"))
theme_set(base_theme)
pt8 <- 8 / ggplot2::.pt
plain_feature <- function(x) {
  x <- sub("Glutamine synthetase \\(GS\\)", "GS", x)
  x <- sub("Glutamate synthase \\(GOGAT\\)", "GOGAT", x)
  x <- sub("Glutamate dehydrogenase \\(GDH\\)", "GDH", x)
  sub("Urease structural genes", "Urease", x)
}
small_number <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    if (v == 0) return("0")
    if (abs(v) < .001) return(formatC(v, format = "e", digits = 0))
    format(signif(v, 2), trim = TRUE, scientific = FALSE)
  }, character(1))
}
diet_boxes <- function(data, mapping) {
  ggplot(data, mapping) +
    geom_boxplot(width = .54, outlier.shape = NA, linewidth = .32,
                 colour = ink, alpha = .65) +
    geom_point(position = position_jitter(width = .1, height = 0, seed = 20261004),
               size = 1.1, shape = 16, colour = ink, alpha = .85) +
    scale_fill_manual(values = diet_colors) +
    scale_y_continuous(labels = small_number, breaks = scales::breaks_pretty(n = 3),
                       expand = expansion(mult = c(.08, .1))) +
    theme(legend.position = "none")
}

p1a <- ggplot(feature_long, aes(percent_all_genes, Feature, fill = Family)) +
  geom_boxplot(orientation = "y", width = .48, outlier.shape = NA,
               linewidth = .3, colour = ink, alpha = .65) +
  geom_point(position = position_jitter(height = .11, width = 0, seed = 20261004),
             size = .62, alpha = .4, colour = ink) +
  scale_fill_manual(values = class_colors) +
  scale_x_log10(breaks = 10^(-5:0), labels = scales::label_log(),
                expand = expansion(mult = c(.035, .08))) +
  scale_y_discrete(labels = plain_feature) +
  labs(x = "Relative abundance (%)", y = NULL) +
  theme(legend.position = "none")
p1b <- diet_boxes(enzyme_long, aes(diet_group, percent_all_genes, fill = diet_group)) +
  facet_wrap(~ Enzyme, ncol = 1, scales = "free_y") +
  labs(x = "Diet", y = "Relative abundance (%)") +
  theme(panel.spacing = grid::unit(1.5, "mm"))
ko_colors <- c(K00262 = "#CC6677", K01915 = "#4477AA",
               K00266 = "#228833", "Other KOs" = "#BBBBBB")
p1c <- ggplot(composition, aes(Share * 100, factor(KO, levels = rev(levels(KO))),
                               fill = KO)) +
  geom_col(width = .62, colour = ink, linewidth = .25) +
  geom_text(aes(label = sprintf("%.2f", Share * 100)), hjust = -.15,
            size = pt8, colour = ink, family = "Arial") +
  scale_fill_manual(values = ko_colors) +
  scale_x_continuous(limits = c(0, 55), breaks = c(0, 20, 40),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Share of nitrogen-metabolism KOs (%)", y = NULL) +
  theme(legend.position = "none")
p1d <- diet_boxes(nitrogen, aes(diet_group, percent_all_genes, fill = diet_group)) +
  annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.15,
           label = sprintf("P = %.3f; q = %.3f", nitrogen_test$P_HC3,
                           nitrogen_test$q_HC3), family = "Arial", size = pt8) +
  scale_y_continuous(labels = small_number, breaks = scales::breaks_pretty(n = 4),
                     expand = expansion(mult = c(.08, .22))) +
  labs(x = "Diet", y = "Nitrogen-metabolism abundance (%)")
figure1 <- p1a + p1b + p1c + p1d +
  plot_layout(design = "AAAABB\nAAAABB\nCCCDDD", heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a")

# Short symbols in cells identify the permutation sensitivity test only.
p2a <- ggplot(carbon, aes(Carbon, Feature, fill = Partial_r)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = sprintf("%.2f", Partial_r),
                colour = ifelse(abs(Partial_r) > .72, "white", ink)),
            family = "Arial", size = pt8) +
  geom_point(data = carbon[carbon$Permutation_support, ],
             aes(x = as.numeric(Carbon) + .33, y = as.numeric(Feature) + .3),
             inherit.aes = FALSE, shape = 21, size = .8, stroke = .25,
             colour = ink, fill = "white") +
  scale_colour_identity() +
  scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
                       limits = c(-1, 1), midpoint = 0,
                       breaks = c(-1, 0, 1), name = "Partial r") +
  scale_x_discrete(labels = c(starch_alpha_glucan = "Starch",
                              cellulose_chain = "Cellulose", xylan = "Xylan",
                              pectin_backbone = "Pectin"), expand = c(0, 0)) +
  scale_y_discrete(labels = plain_feature, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.key.height = grid::unit(2.5, "mm"),
        legend.key.width = grid::unit(8, "mm"))
p2b <- ggplot(gs_pair, aes(C_residual, N_residual)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = ink, linewidth = .45) +
  geom_point(aes(fill = diet_group, shape = diet_group),
             size = 2, stroke = .35, colour = ink) +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
           label = sprintf("r = %.3f\nHC3 q = %.3f\nPermutation q = %.3f",
                           gs_row$Partial_r, gs_row$q_HC3, gs_row$q_permutation),
           size = pt8, family = "Arial", lineheight = 1.1) +
  scale_fill_manual(values = diet_colors) +
  scale_shape_manual(values = diet_shapes) +
  scale_x_continuous(labels = small_number, breaks = scales::breaks_pretty(n = 3),
                     expand = expansion(mult = c(.07, .12))) +
  scale_y_continuous(labels = small_number, breaks = scales::breaks_pretty(n = 3),
                     expand = expansion(mult = c(.06, .32))) +
  labs(x = "Starch-processing residual", y = "GS residual",
       fill = "Diet", shape = "Diet") +
  guides(fill = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = grid::unit(0, "mm"),
        legend.key.width = grid::unit(3, "mm"))
p2c <- ggplot(forest, aes(Beta, Feature_short)) +
  geom_vline(xintercept = 0, linewidth = .3, linetype = "dashed", colour = "#777777") +
  geom_segment(aes(x = CI_lower, xend = CI_upper, yend = Feature_short),
                linewidth = .5, colour = ink) +
  geom_point(aes(fill = Nitrogen_feature), shape = 22, size = 2.2,
             colour = ink, stroke = .3) +
  scale_fill_manual(values = focus_colors) +
  scale_x_continuous(labels = small_number, n.breaks = 4,
                     expand = expansion(mult = c(.05, .1))) +
  scale_y_discrete(labels = function(x) sub("Nitrogen metabolism", "N metabolism", x)) +
  labs(x = "Coefficient (95% CI)", y = NULL) +
  theme(legend.position = "none")
p2d <- ggplot(loo_focus, aes(Beta, Feature_short, fill = Nitrogen_feature)) +
  geom_vline(xintercept = 0, linewidth = .3, linetype = "dashed", colour = "#777777") +
  geom_boxplot(orientation = "y", width = .42, outlier.shape = NA,
               linewidth = .3, colour = ink, alpha = .65) +
  geom_point(position = position_jitter(width = 0, height = .1, seed = 20261004),
             size = .75, alpha = .55, colour = ink) +
  geom_point(data = forest, aes(Beta, Feature_short), inherit.aes = FALSE,
             shape = 23, size = 2.2, fill = "white", colour = ink, stroke = .4) +
  scale_fill_manual(values = focus_colors) +
  scale_y_discrete(labels = function(x) sub("Nitrogen metabolism", "N metabolism", x)) +
  scale_x_continuous(labels = small_number, n.breaks = 4) +
  labs(x = "Leave-one-out coefficient", y = NULL) +
  theme(legend.position = "none")
figure2 <- p2a + p2b + p2c + p2d +
  plot_layout(design = "AB\nAC\nAD", widths = c(1.35, 1), heights = c(1.25, .8, .85)) +
  plot_annotation(tag_levels = "a")

ps1a <- ggplot(pca_scores, aes(PC1, PC2, fill = diet_group, shape = diet_group)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = .3, colour = "#BBBBBB") +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = .3, colour = "#BBBBBB") +
  geom_point(size = 2.6, colour = ink, stroke = .4) +
  scale_fill_manual(values = diet_colors) +
  scale_shape_manual(values = diet_shapes) +
  scale_x_continuous(n.breaks = 5) + scale_y_continuous(n.breaks = 5) +
  labs(x = sprintf("PC1 (%.2f%%)", variance[1] * 100),
       y = sprintf("PC2 (%.2f%%)", variance[2] * 100),
       fill = "Diet", shape = "Diet") +
  theme(legend.position = "right")
ps1b <- ggplot(scree[1:8, ], aes(factor(PC), Variance_percent)) +
  geom_col(fill = "#4477AA", colour = ink, linewidth = .25, width = .65) +
  geom_text(aes(label = sprintf("%.1f", Variance_percent)), vjust = -.4,
            family = "Arial", size = pt8) +
  scale_y_continuous(n.breaks = 4, expand = expansion(mult = c(0, .15))) +
  labs(x = "Principal component", y = "Variance explained (%)")
ps1c <- ggplot(top_load, aes(PC1, Label)) +
  geom_vline(xintercept = 0, linewidth = .3, colour = "#777777") +
  geom_col(width = .62, fill = "#777777") +
  scale_x_continuous(n.breaks = 3) +
  scale_y_discrete(labels = function(x) {
    x <- gsub(" biosynthesis", " synth.", x, fixed = TRUE)
    x
  }) +
  labs(x = "PC1 loading", y = NULL)
figure_s1 <- ps1a + ps1b + ps1c +
  plot_layout(design = "AA\nBC", heights = c(1, 1), widths = c(.85, 1.15)) +
  plot_annotation(tag_levels = "a")

ps2 <- ggplot(volcano, aes(log2_ratio, minus_log10_q, colour = Family, shape = Family)) +
  geom_hline(yintercept = -log10(.05), linetype = "dashed", colour = ink, linewidth = .3) +
  geom_vline(xintercept = 0, colour = "#CCCCCC", linewidth = .25) +
  geom_point(size = 1.8, alpha = .85) +
  ggrepel::geom_text_repel(data = volcano[volcano$Short_label != "", ],
                           aes(label = Short_label), size = pt8, colour = ink,
                           family = "Arial", seed = 20261004, max.overlaps = Inf,
                           min.segment.length = 0, segment.size = .25,
                           box.padding = .25, point.padding = .15) +
  facet_wrap(~ Comparison, ncol = 2) +
  scale_colour_manual(values = class_colors) +
  scale_shape_manual(values = c(Pathway = 16, Enzyme = 17, Ratio = 15)) +
  scale_x_continuous(labels = small_number, n.breaks = 4) +
  scale_y_continuous(limits = c(0, 1.5), breaks = c(0, .5, 1, 1.3),
                     expand = c(0, 0)) +
  labs(x = expression(log[2]~"(adjusted abundance ratio)"),
       y = expression(-log[10]~"(BH q)"), colour = NULL, shape = NULL) +
  theme(legend.position = "bottom", panel.spacing = grid::unit(3, "mm"))
ps3 <- ggplot(heat, aes(Animal, Feature, fill = Z)) +
  geom_tile(colour = "white", linewidth = .12) +
  facet_grid(~ diet_group, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
                       midpoint = 0, name = "z score") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(labels = plain_feature, expand = c(0, 0)) +
  labs(x = "Individual buffalo", y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
        legend.position = "right", panel.spacing.x = grid::unit(1.2, "mm"),
        legend.key.height = grid::unit(7, "mm"),
        legend.key.width = grid::unit(2.5, "mm"))

# Multiplying tiny percentages by a displayed power of ten improves tick legibility.
atlas_style <- atlas
maxima <- tapply(atlas_style$percent_all_genes, atlas_style$Feature_ID, max)
exponents <- floor(log10(maxima))
atlas_style$Display_exponent <- exponents[atlas_style$Feature_ID]
atlas_style$Display_abundance <- atlas_style$percent_all_genes /
  10^atlas_style$Display_exponent
superscript_digits <- c("0" = "\u2070", "1" = "\u00b9", "2" = "\u00b2",
                        "3" = "\u00b3", "4" = "\u2074", "5" = "\u2075",
                        "6" = "\u2076", "7" = "\u2077", "8" = "\u2078",
                        "9" = "\u2079", "-" = "\u207b")
exponent_text <- function(x) {
  paste0("\u00d7 10", paste0(superscript_digits[strsplit(as.character(x), "")[[1]]],
                             collapse = ""))
}
facet_labels <- vapply(features, function(f) {
  label <- plain_feature(labels[[f]])
  label <- gsub(" / ", "/", label, fixed = TRUE)
  label <- gsub(" metabolism", " metab.", label, fixed = TRUE)
  label <- gsub(" biosynthesis", " synth.", label, fixed = TRUE)
  paste0(paste(strwrap(label, width = 22), collapse = "\n"),
         sprintf("\n%s; q = %.3f", exponent_text(exponents[[f]]),
                 diet$q_HC3[match(f, diet$Feature_ID)]))
}, character(1))
atlas_style$Display_facet <- factor(atlas_style$Feature_ID, levels = features)
save_source(atlas_style, "FigureS4_display_scaling.csv")
ps4 <- diet_boxes(atlas_style, aes(diet_group, Display_abundance, fill = diet_group)) +
  facet_wrap(~ Display_facet, ncol = 4, scales = "free_y",
             labeller = labeller(Display_facet = setNames(facet_labels, features))) +
  labs(x = "Diet", y = "Relative abundance (%)") +
  theme(panel.spacing = grid::unit(2.5, "mm"))
