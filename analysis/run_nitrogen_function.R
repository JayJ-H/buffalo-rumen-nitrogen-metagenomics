args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
local_library <- Sys.getenv("TAHP_R_LIBRARY")
if (nzchar(local_library)) .libPaths(c(local_library, .libPaths()))
library(sandwich)
source(file.path(root, "analysis/association_helpers.R"))
out <- file.path(root, "results/nitrogen_function")
cfg <- jsonlite::fromJSON(file.path(out, "definition_freeze.json"), simplifyVector = FALSE)
check_hashes <- function(paths, prefix = root) {
  for (p in names(paths)) {
    stopifnot(digest::digest(file = file.path(prefix, p), algo = "sha256") == paths[[p]])
  }
}
check_hashes(cfg$inputs)
check_hashes(cfg$historical_hashes)
check_hashes(cfg$prepared_files, out)
write_csv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE, na = "")
features <- unlist(cfg$features, use.names = FALSE)
d <- read.csv(file.path(out, "scores_per_animal.csv"), check.names = FALSE,
              colClasses = c(animal_id = "character"))
w <- read.csv(file.path(out, "within_KO_scores_per_animal.csv"), check.names = FALSE)
d$diet <- factor(d$diet_group, levels = LETTERS[1:5])
d$BW010 <- (d$initial_BW_kg - mean(d$initial_BW_kg)) / 10
stopifnot(nrow(d) == 30, length(features) == 23, all(table(d$diet) == 6),
          identical(d$sample_id, w$sample_id))
X <- model.matrix(~ diet + BW010, d)
X0 <- model.matrix(~ BW010, d)
df <- nrow(X) - ncol(X)
stopifnot(df == 24L)
q0 <- qr(X0)
q1 <- qr(X)
A <- qr.solve(X, diag(nrow(X)))
h <- rowSums(qr.Q(q1)^2)
fit_diet <- function(y, design = X, reduced = X0) {
  estimable <- sd(y) > 1e-12
  a <- qr.solve(design, diag(nrow(design)))
  hh <- rowSums(qr.Q(qr(design))^2)
  fit <- lm.fit(design, y)
  b <- fit$coefficients
  e <- fit$residuals
  V <- a %*% diag((e / (1 - hh))^2) %*% t(a)
  j <- 2:5
  estimable <- estimable && rcond(V[j, j]) > 1e-12
  Fstat <- if (estimable) as.numeric(t(b[j]) %*% solve(V[j, j], b[j])) / 4 else NA_real_
  p <- if (estimable) pf(Fstat, 4, nrow(design)-ncol(design), lower.tail = FALSE) else 1
  rss0 <- sum(lm.fit(reduced, y)$residuals^2)
  list(b = b, V = V, F = Fstat, P = p, Estimable = estimable,
       partial_R2 = if (rss0 > 1e-20) (rss0 - sum(e^2))/rss0 else NA_real_)
}
set.seed(as.integer(cfg$seed))
orders <- replicate(as.integer(cfg$permutations), sample.int(30))
within_orders <- replicate(as.integer(cfg$permutations),
                           unlist(lapply(split(seq_len(30), d$diet), sample)))
write_csv(orders, "diet_permutation_indices.csv")
write_csv(within_orders, "within_diet_permutation_indices.csv")
permutation_diet <- function(y) {
  null <- lm.fit(X0, y)
  observed_rss0 <- sum(null$residuals^2)
  observed_rss1 <- sum(lm.fit(X, y)$residuals^2)
  if (observed_rss0 < 1e-20 || observed_rss1 < 1e-20) return(1)
  observed <- ((observed_rss0 - observed_rss1)/4)/(observed_rss1/df)
  ystar <- matrix(null$residuals[orders], 30) + null$fitted.values
  s0 <- colSums(qr.resid(q0, ystar)^2)
  s1 <- colSums(qr.resid(q1, ystar)^2)
  stats <- ((s0-s1)/4)/(s1/df)
  (1 + sum(stats >= observed - 1e-12))/(ncol(orders)+1)
}
checks <- list()
diet_family <- function(input, denominator) {
  rows <- lapply(features, function(name) {
    y <- input[[name]]
    f <- fit_diet(y)
    if (f$Estimable) {
      model <- lm(y ~ diet + BW010, d)
      err <- max(abs(f$V - vcovHC(model, type = "HC3")))
      stopifnot(err < 1e-8)
      checks[[length(checks)+1]] <<- data.frame(Test = paste(denominator, name), Max_error = err)
    }
    data.frame(Feature_ID = name, N = 30, df1 = 4, df2 = 24, Denominator = denominator,
               Estimable = f$Estimable, F_HC3 = f$F, P_HC3 = f$P,
               P_permutation = if (f$Estimable) permutation_diet(y) else 1,
               partial_R2 = f$partial_R2)
  })
  result <- do.call(rbind, rows)
  result$q_HC3 <- p.adjust(result$P_HC3, "BH", n = 23)
  result$q_permutation <- p.adjust(result$P_permutation, "BH", n = 23)
  result
}
primary <- diet_family(d, "all_gene_mass")
sensitivity <- diet_family(w, "single_KO_annotated_gene_mass")
write_csv(primary, "diet_omnibus.csv")
write_csv(sensitivity, "diet_within_KO_sensitivity.csv")
grids <- data.frame(diet = factor(LETTERS[1:5], levels = LETTERS[1:5]), BW010 = 0)
M <- model.matrix(~ diet + BW010, grids)
pairs <- combn(LETTERS[1:5], 2, simplify = FALSE)
contrast_rows <- adjusted <- distributions <- loo_diet <- list()
for (name in features) {
  y <- d[[name]]
  f <- fit_diet(y)
  for (pair in pairs) {
    cvec <- M[match(pair[2], LETTERS[1:5]), ] - M[match(pair[1], LETTERS[1:5]), ]
    b <- sum(cvec * f$b)
    se <- sqrt(max(0, as.numeric(t(cvec) %*% f$V %*% cvec)))
    contrast_rows[[length(contrast_rows)+1]] <- data.frame(
      Feature_ID = name, reference = pair[1], comparison = pair[2],
      log_ratio = b, SE_HC3 = se, ratio = exp(b),
      ratio_CI_low = exp(b-qt(.975, df)*se), ratio_CI_high = exp(b+qt(.975, df)*se),
      P = if (f$Estimable && se > 0) 2*pt(abs(b/se), df, lower.tail = FALSE) else 1)
  }
  unit <- if (name %in% names(cfg$ratios)) "annotation_mass_ratio" else "percent_all_genes"
  multiplier <- if (unit == "annotation_mass_ratio") 1 else 100
  for (i in seq_len(5)) {
    b <- sum(M[i, ]*f$b)
    se <- sqrt(max(0, as.numeric(t(M[i, ]) %*% f$V %*% M[i, ])))
    adjusted[[length(adjusted)+1]] <- data.frame(
      Feature_ID = name, diet_group = LETTERS[i], unit = unit,
      geometric_mean = exp(b)*multiplier,
      CI_lower = exp(b-qt(.975, df)*se)*multiplier,
      CI_upper = exp(b+qt(.975, df)*se)*multiplier)
    values <- exp(y[d$diet == LETTERS[i]])*multiplier
    distributions[[length(distributions)+1]] <- data.frame(
      Feature_ID = name, diet_group = LETTERS[i], n = 6, unit = unit,
      mean = mean(values), SD = sd(values), median = median(values),
      min = min(values), max = max(values))
  }
  for (i in seq_len(30)) {
    z <- fit_diet(y[-i], X[-i, ], X0[-i, ])
    loo_diet[[length(loo_diet)+1]] <- data.frame(
      Feature_ID = name, excluded_animal = d$animal_id[i], F_HC3 = z$F, P_HC3 = z$P,
      B_log_ratio = z$b[2], C_log_ratio = z$b[3], D_log_ratio = z$b[4], E_log_ratio = z$b[5])
  }
}
contrasts <- do.call(rbind, contrast_rows)
contrasts$q_all230 <- p.adjust(contrasts$P, "BH", n = 230)
write_csv(contrasts, "diet_all_contrasts.csv")
write_csv(do.call(rbind, adjusted), "diet_adjusted_means.csv")
write_csv(do.call(rbind, distributions), "diet_group_summary.csv")
write_csv(do.call(rbind, loo_diet), "diet_LOO.csv")

pair_members <- read.csv(file.path(out, "carbon_pair_membership.csv"), check.names = FALSE)
pair_members$Estimable <- tolower(as.character(pair_members$Estimable)) == "true"
pair_scores <- read.csv(file.path(out, "carbon_pair_scores.csv"), check.names = FALSE,
                       colClasses = c(animal_id = "character"))
stopifnot(nrow(pair_members) == 80)
set.seed(as.integer(cfg$seed) + 1L)
boot_indices <- t(replicate(as.integer(cfg$bootstrap_draws),
                   unlist(lapply(split(seq_len(30), d$diet), function(v) sample(v, 6, TRUE)))))
write_csv(boot_indices, "bootstrap_indices.csv")
carbon_rows <- carbon_loo <- list()
for (pair_id in pair_members$Pair_ID) {
  z <- pair_scores[pair_scores$Pair_ID == pair_id, ]
  stopifnot(identical(z$sample_id, d$sample_id))
  z$diet <- d$diet
  z$BW010 <- d$BW010
  full <- model.matrix(~ carbon + diet + BW010 + log_KO_annotation_load, z)
  null <- model.matrix(~ diet + BW010 + log_KO_annotation_load, z)
  without_load <- model.matrix(~ carbon + diet + BW010, z)
  estimable <- pair_members$Estimable[match(pair_id, pair_members$Pair_ID)] &&
    qr(full)$rank == ncol(full) && sd(qr.resid(qr(null), z$carbon)) > 1e-12
  if (!estimable) {
    carbon_rows[[length(carbon_rows)+1]] <- data.frame(
      Pair_ID = pair_id, Estimable = FALSE, N = 30, Residual_df = 22,
      Beta = NA, SE_HC3 = NA, CI_lower = NA, CI_upper = NA, P_HC3 = 1,
      Beta_standardized = NA, Partial_r = NA, P_permutation = 1,
      Beta_no_load = NA, P_no_load = 1, LOO_sign_agreement = NA,
      LOO_beta_min = NA, LOO_beta_max = NA, Bootstrap_success = 0,
      Bootstrap_positive_fraction = NA, Bootstrap_CI_lower = NA, Bootstrap_CI_upper = NA)
    next
  }
  f <- robust_matrix_fit(full, matrix(z$response, ncol = 1, dimnames = list(NULL, pair_id)),
                         term = "carbon")$results
  ff <- robust_matrix_fit(without_load,
                          matrix(z$response, ncol = 1, dimnames = list(NULL, pair_id)),
                          term = "carbon")$results
  model <- lm(response ~ carbon + diet + BW010 + log_KO_annotation_load, z)
  ct <- lmtest::coeftest(model, vcov. = vcovHC(model, type = "HC3"), df = 22)
  err <- max(abs(c(f$Beta_per_100g_ADG - ct["carbon",1],
                   f$SE_HC3 - ct["carbon",2], f$P - ct["carbon",4])))
  stopifnot(err < 1e-7)
  checks[[length(checks)+1]] <- data.frame(Test = pair_id, Max_error = err)
  nullfit <- lm.fit(null, z$response)
  ystar <- matrix(nullfit$residuals[within_orders], 30) + nullfit$fitted.values
  inv <- solve(crossprod(full))
  aa <- inv %*% t(full)
  bperm <- as.numeric(aa[2, ] %*% ystar)
  rss <- colSums(qr.resid(qr(full), ystar)^2)
  tperm <- bperm / sqrt((rss/22)*inv[2,2])
  obs <- coef(model)["carbon"] / sqrt(sum(residuals(model)^2)/22*inv[2,2])
  pperm <- (1 + sum(abs(tperm) >= abs(obs)-1e-12))/(ncol(within_orders)+1)
  loo <- sapply(seq_len(30), function(i) lm.fit(full[-i, ], z$response[-i])$coefficients[2])
  carbon_loo[[length(carbon_loo)+1]] <- data.frame(
    Pair_ID = pair_id, excluded_animal = d$animal_id, Beta = loo)
  boot <- apply(boot_indices, 1, function(ids) {
    if (qr(full[ids, ])$rank < ncol(full)) return(NA_real_)
    lm.fit(full[ids, ], z$response[ids])$coefficients[2]
  })
  yr <- qr.resid(qr(null), z$response)
  xr <- qr.resid(qr(null), z$carbon)
  carbon_rows[[length(carbon_rows)+1]] <- data.frame(
    Pair_ID = pair_id, Estimable = TRUE, N = 30, Residual_df = 22,
    Beta = f$Beta_per_100g_ADG, SE_HC3 = f$SE_HC3, CI_lower = f$CI_lower, CI_upper = f$CI_upper,
    P_HC3 = f$P, Beta_standardized = f$Beta_per_100g_ADG * sd(z$carbon)/sd(z$response),
    Partial_r = cor(xr, yr), P_permutation = pperm,
    Beta_no_load = ff$Beta_per_100g_ADG, P_no_load = ff$P,
    LOO_sign_agreement = mean(sign(loo) == sign(f$Beta_per_100g_ADG)),
    LOO_beta_min = min(loo), LOO_beta_max = max(loo),
    Bootstrap_success = sum(is.finite(boot)),
    Bootstrap_positive_fraction = mean(boot > 0, na.rm = TRUE),
    Bootstrap_CI_lower = unname(quantile(boot, .025, na.rm = TRUE)),
    Bootstrap_CI_upper = unname(quantile(boot, .975, na.rm = TRUE)))
}
carbon <- do.call(rbind, carbon_rows)
carbon$q_HC3 <- p.adjust(carbon$P_HC3, "BH", n = 80)
carbon$q_permutation <- p.adjust(carbon$P_permutation, "BH", n = 80)
carbon$q_no_load <- p.adjust(carbon$P_no_load, "BH", n = 80)
carbon$Supported_HC3_and_permutation <- carbon$q_HC3 < .05 & carbon$q_permutation < .05
carbon <- merge(pair_members[, c("Pair_ID", "Nitrogen_feature", "Carbon_feature", "Shared_KOs")],
                carbon, by = "Pair_ID", sort = FALSE)
carbon <- carbon[match(pair_members$Pair_ID, carbon$Pair_ID), ]
write_csv(carbon, "carbon_function_associations.csv")
write_csv(do.call(rbind, carbon_loo), "carbon_LOO.csv")
write_csv(do.call(rbind, checks), "R_package_checks.csv")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
check_hashes(cfg$historical_hashes)
check_hashes(cfg$inputs)
summary <- list(
  status = "PASS", n_animals = 30, diet_responses = 23, contrasts = 230,
  carbon_pairs = 80, diet_HC3_supported = sum(primary$q_HC3 < .05),
  diet_permutation_supported = sum(primary$q_permutation < .05),
  contrasts_supported = sum(contrasts$q_all230 < .05),
  carbon_HC3_supported = sum(carbon$q_HC3 < .05),
  carbon_both_supported = sum(carbon$Supported_HC3_and_permutation),
  max_package_error = max(vapply(checks, function(z) z$Max_error, numeric(1))),
  historical_hashes_unchanged = TRUE, membership_redefined = FALSE)
jsonlite::write_json(summary, file.path(out, "run_validation.json"), pretty = TRUE, auto_unbox = TRUE)
cat("\nDiet: all 23 responses\n")
print(primary[, c("Feature_ID", "P_HC3", "q_HC3", "q_permutation")], row.names = FALSE)
cat("\nCarbon-function: first 12 by HC3 q, all 80 saved\n")
print(head(carbon[order(carbon$q_HC3), c("Nitrogen_feature", "Carbon_feature",
      "Partial_r", "q_HC3", "q_permutation", "LOO_sign_agreement")], 12), row.names = FALSE)
cat("\nRun status PASS; independent Python validation remains required.\n")
