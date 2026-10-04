## =============================================================================
## ME317 Project 2026 -- Task 3: Portfolio risk and backtesting
##
## Assumes in the environment (from Tasks 1-2):
##   prices_clean    : data.frame, date + 10 adjusted price columns
##   log_returns     : data.frame, date + 10 log-return columns
##   simple_returns  : data.frame, date + 10 simple-return columns
##
## Produces:
##   portfolio       : date, value, loss, loss_rel
##   risk_table      : VaR and ES at 95% / 99% for three methods
##   backtest_table  : exceedance counts and Kupiec POF results
##
## NOTE: not executed in this environment. Run section by section and check the
##       printed diagnostics before using any output.
## =============================================================================

## ---- 0. Setup ----------------------------------------------------------------

set.seed(12345)

if (!requireNamespace("MASS", quietly = TRUE))
  stop("Install MASS: install.packages('MASS')")
library(MASS)   # fitdistr()

stopifnot(exists("prices_clean"), "date" %in% names(prices_clean))

dir.create("data",    showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

## ---- 1. Parameters -----------------------------------------------------------

INITIAL   <- 1000                       # dollars invested per stock at t = 0
ALPHAS    <- c(0.95, 0.99)              # VaR/ES confidence levels
IN_END    <- as.Date("2020-12-31")      # end of estimation window
OUT_START <- as.Date("2021-01-01")      # start of backtest window

## ---- 2. Portfolio construction -----------------------------------------------

dates   <- as.Date(prices_clean$date)
P       <- as.matrix(prices_clean[, -1, drop = FALSE])
TICKERS <- colnames(P)
n       <- nrow(P)

## $1000 per stock buys a FIXED number of shares. "No rebalancing" means the
## share counts are constant; the weights are therefore equal only on day one
## and drift thereafter with relative performance.
shares <- INITIAL / P[1, ]
names(shares) <- TICKERS

## Initial weights: $1000 / $10,000 = 0.10 each, by construction.
w0 <- rep(INITIAL, length(TICKERS)) / (INITIAL * length(TICKERS))
names(w0) <- TICKERS
cat("Initial weights:\n"); print(round(w0, 4))
cat("Initial portfolio value: $", format(sum(INITIAL * length(TICKERS))), "\n\n", sep = "")

## Portfolio value path.
V <- as.numeric(P %*% shares)

## Weight drift: comparing day-one to final weights shows how far the portfolio
## moved from equal weighting without any trading. Relevant because it changes
## the portfolio's effective sector exposure over the sample.
W <- sweep(P, 2, shares, "*") / V       # each row divided by that day's value
cat("Weights on first vs last day:\n")
print(round(rbind(first = W[1, ], last = W[n, ]), 4))

## Daily loss, in dollars: L_t = -(V_t - V_{t-1}) = V_{t-1} - V_t.
## Positive values are losses, which is the course's sign convention.
V_prev <- V[-n]
V_curr <- V[-1]
loss   <- V_prev - V_curr
loss_dates <- dates[-1]

## Relative loss, used for the scale-consistent robustness check in Section 7.
loss_rel <- loss / V_prev

portfolio <- data.frame(date = loss_dates, value = V_curr,
                        loss = loss, loss_rel = loss_rel, row.names = NULL)

## ---- 3. Split into estimation and backtest windows ---------------------------

in_idx  <- loss_dates <= IN_END
out_idx <- loss_dates >= OUT_START

L_in  <- loss[in_idx];      L_out  <- loss[out_idx]
Lr_in <- loss_rel[in_idx];  Lr_out <- loss_rel[out_idx]
V_out_prev <- V_prev[out_idx]

T_in  <- length(L_in)
T_out <- length(L_out)

cat("\nEstimation window: ", format(min(loss_dates[in_idx])), " to ",
    format(max(loss_dates[in_idx])), " (", T_in, " days)\n", sep = "")
cat("Backtest window:   ", format(min(loss_dates[out_idx])), " to ",
    format(max(loss_dates[out_idx])), " (", T_out, " days)\n", sep = "")
cat("Expected exceedances out-of-sample: ",
    paste(sprintf("%.0f%%: %.1f", 100 * ALPHAS, (1 - ALPHAS) * T_out),
          collapse = " | "), "\n\n", sep = "")

## ---- 4. Method (i): Empirical / historical simulation ------------------------
## VaR is the empirical quantile of the in-sample loss distribution; ES is the
## mean of losses at or beyond it.
##
## CAVEAT to report: at 99% with ~250 observations the exceedance set contains
## only about 2-3 losses, so empirical ES at 99% is an average over a handful of
## points and is very imprecise. The count is printed so you can say so.

empirical_var_es <- function(L, alphas) {
  v <- as.numeric(quantile(L, probs = alphas, type = 7))
  e <- sapply(seq_along(alphas), function(k) mean(L[L >= v[k]]))
  n_tail <- sapply(seq_along(alphas), function(k) sum(L >= v[k]))
  list(VaR = v, ES = as.numeric(e), n_tail = n_tail)
}

emp <- empirical_var_es(L_in, ALPHAS)
cat("Observations in the empirical tail (drives ES precision): ",
    paste(sprintf("%.0f%%: %d", 100 * ALPHAS, emp$n_tail), collapse = " | "),
    "\n", sep = "")

## ---- 5. Method (ii): Normal model -------------------------------------------
## L ~ N(mu, sigma^2), fitted by moment matching (equal to MLE here).
##   VaR_a = mu + sigma * qnorm(a)
##   ES_a  = mu + sigma * dnorm(qnorm(a)) / (1 - a)

mu_n    <- mean(L_in)
sigma_n <- sd(L_in)

var_norm <- mu_n + sigma_n * qnorm(ALPHAS)
es_norm  <- mu_n + sigma_n * dnorm(qnorm(ALPHAS)) / (1 - ALPHAS)

cat("\nNormal fit:  mu = ", round(mu_n, 3), " | sigma = ", round(sigma_n, 3),
    "\n", sep = "")

## ---- 6. Method (iii): Student-t model ---------------------------------------
## fitdistr(densfun = "t") estimates a location-scale t: (x - m)/s ~ t_df.
##
## Two implementation points:
##  (a) The losses are in dollars (magnitudes in the hundreds), which makes the
##      optimiser poorly scaled. Fitting on standardised data and transforming
##      back is numerically far more stable: if y = x/c and y ~ t(m, s, df),
##      then x ~ t(c*m, c*s, df). The df is scale-invariant.
##  (b) `s` is a SCALE, not a standard deviation. For df > 2,
##      sd = s * sqrt(df / (df - 2)). Comparing s directly against sigma_n
##      would understate the fitted dispersion.

c_scale <- sd(L_in)
y       <- L_in / c_scale

fit_t <- tryCatch(
  fitdistr(y, densfun = "t",
           start = list(m = median(y), s = IQR(y) / 1.35, df = 5),
           lower = c(m = -Inf, s = 1e-6, df = 2.1)),   # df > 2 for finite variance
  error = function(e) {
    warning("fitdistr failed (", conditionMessage(e),
            "); retrying with a profile over df.", call. = FALSE)
    NULL
  }
)

if (is.null(fit_t)) {
  ## Fallback: profile the log-likelihood over a grid of df, optimising m and s
  ## at each value, then take the maximum. Slower but robust.
  grid <- seq(2.5, 30, by = 0.1)
  ll <- sapply(grid, function(d) {
    f <- try(fitdistr(y, densfun = function(x, m, s) dt((x - m) / s, df = d) / s,
                      start = list(m = median(y), s = IQR(y) / 1.35),
                      lower = c(-Inf, 1e-6)), silent = TRUE)
    if (inherits(f, "try-error")) -Inf else f$loglik
  })
  d_hat <- grid[which.max(ll)]
  f_hat <- fitdistr(y, densfun = function(x, m, s) dt((x - m) / s, df = d_hat) / s,
                    start = list(m = median(y), s = IQR(y) / 1.35),
                    lower = c(-Inf, 1e-6))
  est <- c(m = f_hat$estimate["m"], s = f_hat$estimate["s"], df = d_hat)
} else {
  est <- fit_t$estimate
}

## Transform back to the dollar scale.
mu_t <- as.numeric(est["m"])  * c_scale
s_t  <- as.numeric(est["s"])  * c_scale
nu_t <- as.numeric(est["df"])

sd_t <- if (nu_t > 2) s_t * sqrt(nu_t / (nu_t - 2)) else NA_real_

cat("Student-t fit:  mu = ", round(mu_t, 3),
    " | scale s = ", round(s_t, 3),
    " | nu = ", round(nu_t, 3),
    " | implied sd = ", round(sd_t, 3), "\n", sep = "")
if (nu_t <= 4)
  message("Note: nu <= 4 implies infinite kurtosis; the fourth moment does not ",
          "exist, so sample kurtosis is not a valid check on this fit.")

## VaR and ES for a location-scale t (McNeil-Frey-Embrechts):
##   VaR_a = mu + s * t_nu^{-1}(a)
##   ES_a  = mu + s * [ g_nu(t_a) / (1-a) ] * (nu + t_a^2) / (nu - 1)
## where g_nu is the standard t density and t_a = qt(a, nu). Requires nu > 1.

t_var_es <- function(mu, s, nu, alphas) {
  q  <- qt(alphas, df = nu)
  v  <- mu + s * q
  e  <- mu + s * (dt(q, df = nu) / (1 - alphas)) * ((nu + q^2) / (nu - 1))
  list(VaR = v, ES = e)
}
tt <- t_var_es(mu_t, s_t, nu_t, ALPHAS)

## ---- 7. Risk comparison table -----------------------------------------------

risk_table <- data.frame(
  method = rep(c("Empirical", "Normal", "Student-t"), each = length(ALPHAS)),
  alpha  = rep(ALPHAS, times = 3),
  VaR    = c(emp$VaR,  var_norm, tt$VaR),
  ES     = c(emp$ES,   es_norm,  tt$ES),
  row.names = NULL
)
risk_table$ES_VaR_ratio <- risk_table$ES / risk_table$VaR

## Express as a percentage of the portfolio value on the last in-sample day,
## which makes the magnitudes interpretable.
V_ref <- V[max(which(dates <= IN_END))]
risk_table$VaR_pct <- 100 * risk_table$VaR / V_ref
risk_table$ES_pct  <- 100 * risk_table$ES  / V_ref

cat("\n--- VaR and ES, estimated on 2020 (dollars; V_ref = $",
    round(V_ref, 0), ") ---\n", sep = "")
print(transform(risk_table,
                VaR = round(VaR, 2), ES = round(ES, 2),
                ES_VaR_ratio = round(ES_VaR_ratio, 3),
                VaR_pct = round(VaR_pct, 3), ES_pct = round(ES_pct, 3)),
      row.names = FALSE)

## ---- 8. Kupiec Proportion of Failures test ----------------------------------
## Null: the true exceedance probability equals p = 1 - alpha.
##   LR = 2 * [ x*log(p_hat/p) + (T-x)*log((1-p_hat)/(1-p)) ],   p_hat = x/T
## Asymptotically chi-squared with 1 df.
##
## The x = 0 and x = T cases must be handled separately: p_hat = 0 makes
## log(p_hat/p) equal -Inf, and 0 * -Inf evaluates to NaN in R rather than to
## its limit of 0. Taking limits gives LR = -2*T*log(1-p) when x = 0 and
## LR = -2*T*log(p) when x = T. With a calm backtest year and models calibrated
## on 2020, x = 0 at the 99% level is a realistic outcome, so this matters.

kupiec_pof <- function(x, T, alpha) {
  p    <- 1 - alpha
  phat <- x / T
  if (x == 0) {
    LR <- -2 * T * log(1 - p)
  } else if (x == T) {
    LR <- -2 * T * log(p)
  } else {
    LR <- 2 * (x * log(phat / p) + (T - x) * log((1 - phat) / (1 - p)))
  }
  data.frame(violations = x, T = T, p_expected = p, p_hat = phat,
             expected_violations = p * T,
             LR_stat = LR,
             p_value = pchisq(LR, df = 1, lower.tail = FALSE),
             reject_5pct = pchisq(LR, df = 1, lower.tail = FALSE) < 0.05,
             row.names = NULL)
}

## ---- 9. Backtest: 2020 parameters applied to 2021 ---------------------------
## The VaR forecast is UNCONDITIONAL: one number per method per level, held
## constant across 2021. There is no volatility updating, which is precisely the
## limitation the exceedance pattern should reveal.

backtest_table <- do.call(rbind, lapply(seq_len(nrow(risk_table)), function(i) {
  m <- risk_table$method[i]; a <- risk_table$alpha[i]; v <- risk_table$VaR[i]
  x <- sum(L_out > v)
  cbind(data.frame(method = m, alpha = a, VaR_forecast = v),
        kupiec_pof(x, T_out, a))
}))
rownames(backtest_table) <- NULL

cat("\n--- Backtest on 2021 (dollar losses vs fixed 2020 VaR) ---\n")
print(transform(backtest_table,
                VaR_forecast = round(VaR_forecast, 2),
                p_hat = round(p_hat, 4),
                expected_violations = round(expected_violations, 1),
                LR_stat = round(LR_stat, 3),
                p_value = signif(p_value, 4)),
      row.names = FALSE)

## Dates of the exceedances, per method and level. Useful for the write-up: if
## they cluster, the unconditional model is missing volatility dynamics, which
## Kupiec alone cannot detect (that is what Christoffersen's independence test
## adds).
cat("\n--- Exceedance dates ---\n")
for (i in seq_len(nrow(risk_table))) {
  v  <- risk_table$VaR[i]
  hit <- loss_dates[out_idx][L_out > v]
  cat(sprintf("%-10s %.0f%%: %d exceedance(s)%s\n",
              risk_table$method[i], 100 * risk_table$alpha[i], length(hit),
              if (length(hit)) paste0(" -- ", paste(format(hit), collapse = ", ")) else ""))
}

## ---- 10. Robustness check: scale-consistent (relative) VaR ------------------
## The dollar VaR above is estimated on a ~$10,000 portfolio in 2020 and applied
## in 2021, when the portfolio is worth a different amount. A fixed dollar
## threshold is therefore not a fixed risk threshold. Estimating VaR as a
## fraction of portfolio value and rescaling by the previous day's value each
## day removes that mismatch. Report this alongside the primary result; if the
## exceedance counts differ materially, the relative version is the defensible
## one.

emp_r  <- empirical_var_es(Lr_in, ALPHAS)
mu_rn  <- mean(Lr_in); sd_rn <- sd(Lr_in)
var_rn <- mu_rn + sd_rn * qnorm(ALPHAS)

cr    <- sd(Lr_in)
fit_r <- tryCatch(
  fitdistr(Lr_in / cr, densfun = "t",
           start = list(m = median(Lr_in / cr), s = IQR(Lr_in / cr) / 1.35, df = 5),
           lower = c(m = -Inf, s = 1e-6, df = 2.1)),
  error = function(e) NULL)

if (!is.null(fit_r)) {
  var_rt <- as.numeric(fit_r$estimate["m"]) * cr +
    as.numeric(fit_r$estimate["s"]) * cr *
    qt(ALPHAS, df = as.numeric(fit_r$estimate["df"]))
} else {
  var_rt <- rep(NA_real_, length(ALPHAS))
  warning("Relative-scale t fit failed; robustness check incomplete.", call. = FALSE)
}

rel_var <- data.frame(method = rep(c("Empirical", "Normal", "Student-t"),
                                   each = length(ALPHAS)),
                      alpha  = rep(ALPHAS, times = 3),
                      VaR_rel = c(emp_r$VaR, var_rn, var_rt))

backtest_rel <- do.call(rbind, lapply(seq_len(nrow(rel_var)), function(i) {
  vr <- rel_var$VaR_rel[i]
  if (is.na(vr)) return(NULL)
  ## Daily dollar threshold = relative VaR x previous day's portfolio value.
  x <- sum(L_out > vr * V_out_prev)
  cbind(data.frame(method = rel_var$method[i], alpha = rel_var$alpha[i],
                   VaR_rel_pct = 100 * vr),
        kupiec_pof(x, T_out, rel_var$alpha[i]))
}))
rownames(backtest_rel) <- NULL

cat("\n--- Robustness: relative (percentage-of-value) VaR, rescaled daily ---\n")
print(transform(backtest_rel,
                VaR_rel_pct = round(VaR_rel_pct, 3),
                p_hat = round(p_hat, 4),
                expected_violations = round(expected_violations, 1),
                LR_stat = round(LR_stat, 3),
                p_value = signif(p_value, 4)),
      row.names = FALSE)

## ---- 11. Diagnostic: volatility across the two windows ----------------------
## Interpreting the backtest requires knowing how much the loss distribution
## itself changed. If 2021 is much calmer than 2020, too FEW exceedances is the
## expected failure mode, and Kupiec rejects for over-conservatism rather than
## for under-estimating risk. That is still a model-risk finding -- capital held
## against a stale crisis-calibrated model is capital not deployed.

cat("\n--- Loss distribution by window ---\n")
print(data.frame(
  window   = c("In-sample 2020", "Out-of-sample 2021"),
  n        = c(T_in, T_out),
  mean     = round(c(mean(L_in),  mean(L_out)), 3),
  sd       = round(c(sd(L_in),    sd(L_out)), 3),
  worst    = round(c(max(L_in),   max(L_out)), 2),
  q95      = round(c(quantile(L_in, .95), quantile(L_out, .95)), 2),
  row.names = NULL))

cat("\nRatio of out-of-sample to in-sample loss volatility: ",
    round(sd(L_out) / sd(L_in), 3), "\n", sep = "")

## ---- 12. Save ----------------------------------------------------------------

write.csv(portfolio,      "results/portfolio_value_loss.csv", row.names = FALSE)
write.csv(risk_table,     "results/task3_var_es.csv",         row.names = FALSE)
write.csv(backtest_table, "results/task3_backtest.csv",       row.names = FALSE)
write.csv(backtest_rel,   "results/task3_backtest_relative.csv", row.names = FALSE)

saveRDS(list(portfolio = portfolio, shares = shares, w0 = w0,
             risk_table = risk_table, backtest_table = backtest_table,
             backtest_rel = backtest_rel,
             fit_normal = c(mu = mu_n, sigma = sigma_n),
             fit_t = c(mu = mu_t, scale = s_t, nu = nu_t, sd = sd_t),
             L_in = L_in, L_out = L_out),
        "data/task3_results.rds")

message("Task 3 results written to results/ and data/task3_results.rds")