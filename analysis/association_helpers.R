# Multivariate QR fits share the same design across features. HC3 is checked
# against sandwich::vcovHC and lmtest::coeftest before project models run.
robust_matrix_fit <- function(X, Y, term = "ADG100") {
  Y <- as.matrix(Y)
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(nrow(Y) == n, qr(X)$rank == p, n > p + 1L)
  fitted <- lm.fit(X, Y)
  B <- matrix(fitted$coefficients, nrow = p,
              dimnames = list(colnames(X), colnames(Y)))
  E <- matrix(fitted$residuals, nrow = n, dimnames = dimnames(Y))
  h <- rowSums(qr.Q(fitted$qr)[, seq_len(p), drop = FALSE]^2)
  stopifnot(all(h < 1 - 1e-10))
  A <- qr.solve(X, diag(n))
  j <- match(term, colnames(X))
  stopifnot(!is.na(j))
  se <- sqrt(colSums(sweep(E / (1 - h), 1, A[j, ], "*")^2))
  stopifnot(all(is.finite(se)), all(se > 0))
  beta <- B[j, ]
  df <- n - p
  P <- 2 * pt(abs(beta / se), df = df, lower.tail = FALSE)
  mse <- colSums(E^2) / df
  cooks <- sweep(E^2, 2, p * mse, "/") * h / (1 - h)^2
  deleted_var <- sweep(-E^2 / (1 - h), 2, colSums(E^2), "+") / (df - 1)
  stopifnot(all(deleted_var > 0))
  student <- E / sqrt(deleted_var * (1 - h))
  result <- data.frame(
    Feature_ID = colnames(Y), N = n, Residual_df = df,
    Beta_per_100g_ADG = as.numeric(beta), SE_HC3 = as.numeric(se),
    CI_lower = as.numeric(beta - qt(.975, df) * se),
    CI_upper = as.numeric(beta + qt(.975, df) * se), P = as.numeric(P),
    Max_Cooks_distance = apply(cooks, 2, max),
    N_Cooks_above_4_n = colSums(cooks > 4 / n),
    Max_abs_external_studentized_residual = apply(abs(student), 2, max),
    stringsAsFactors = FALSE, row.names = NULL
  )
  list(results = result, B = B, E = E, leverage = h, cooks = cooks)
}

multiplicative_clr <- function(P, delta) {
  stopifnot(all(P >= 0), max(abs(colSums(P) - 1)) < 1e-10)
  zeros <- P == 0
  nz <- colSums(zeros)
  stopifnot(all(nz * delta < 1))
  replaced <- sweep(P, 2, 1 - nz * delta, "*")
  replaced[zeros] <- delta
  stopifnot(all(replaced > 0), max(abs(colSums(replaced) - 1)) < 1e-10)
  logged <- log(replaced)
  sweep(logged, 2, colMeans(logged), "-")
}

validate_helpers <- function() {
  set.seed(20261002)
  d <- data.frame(ADG100 = rnorm(30), diet = factor(rep(LETTERS[1:5], each = 6)),
                  BW010 = rnorm(30))
  X <- model.matrix(~ ADG100 + diet + BW010, d)
  Y <- matrix(rnorm(30 * 8), 30, 8)
  Y <- Y * (1 + abs(d$ADG100)) + d$ADG100
  colnames(Y) <- paste0("fixture", seq_len(ncol(Y)))
  fit <- robust_matrix_fit(X, Y)
  errors <- numeric()
  for (k in seq_len(ncol(Y))) {
    m <- lm(Y[, k] ~ ADG100 + diet + BW010, data = d)
    ct <- lmtest::coeftest(m, vcov. = sandwich::vcovHC(m, type = "HC3"),
                          df = df.residual(m))
    errors <- c(errors,
                abs(fit$results$Beta_per_100g_ADG[k] - coef(m)["ADG100"]),
                abs(fit$results$SE_HC3[k] - ct["ADG100", 2]),
                abs(fit$results$P[k] - ct["ADG100", 4]),
                abs(fit$results$Max_Cooks_distance[k] - max(cooks.distance(m))),
                abs(fit$results$Max_abs_external_studentized_residual[k] -
                      max(abs(rstudent(m)))))
  }
  stopifnot(max(errors) < 1e-9)
  P <- matrix(c(0, .2, .8, .1, 0, .9), 3, 2)
  transformed <- multiplicative_clr(P, .01)
  expected <- log(c(.01, .198, .792))
  stopifnot(max(abs(transformed[, 1] - (expected - mean(expected)))) < 1e-12)
  stopifnot(max(abs(colMeans(transformed))) < 1e-12)
  for (i in c(1L, 12L, 30L)) {
    reduced <- lm.fit(X[-i, ], Y[-i, ])$coefficients["ADG100", ]
    for (k in seq_len(ncol(Y))) {
      m <- lm(Y[-i, k] ~ ADG100 + diet + BW010, data = d[-i, ])
      stopifnot(abs(reduced[k] - coef(m)["ADG100"]) < 1e-10)
    }
  }
  list(status = "PASS", synthetic_features = ncol(Y),
       max_error_against_lm_sandwich_lmtest = max(errors),
       zero_replacement_fixture = "PASS", LOO_fixture = "PASS")
}
