## =============================================================================
## ME317 Project 2026 -- Task 2: Returns, extreme events, and (non-)normality
##
## Assumes `prices_clean` is already in the environment:
##   data.frame, column 1 = date, columns 2:11 = adjusted prices for
##   TPR, PFE, GILD, MSFT, AMZN, MAR, AAL, BA, BP, AXP
##
## Produces:
##   log_returns, simple_returns : data.frames (date + 10 return columns)
##   stats_table                 : moments + normality test results
##   extremes                    : top 5 losses and gains per stock
##   top_moves                   : largest moves overall, matched to event dates
##
## NOTE: not executed in this environment. Run section by section and read the
##       printed output before using any of it in the report.
## =============================================================================

## ---- 0. Setup ----------------------------------------------------------------

set.seed(12345)

## moments  : skewness(), kurtosis()
## tseries  : jarque.bera.test()   -- the course's normality test (Workshop 7)
## QRM      : MardiaTest()         -- joint normality (optional but recommended)
for (p in c("moments", "tseries")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Install ", p, ": install.packages('", p, "')")
}
library(moments)
library(tseries)

has_QRM <- requireNamespace("QRM", quietly = TRUE)
if (has_QRM) library(QRM) else
  message("QRM not installed -- the Mardia joint-normality test will be skipped.")

stopifnot(exists("prices_clean"), "date" %in% names(prices_clean))

## ---- 1. Extract matrices and define the sector map ---------------------------

dates   <- as.Date(prices_clean$date)
P       <- as.matrix(prices_clean[, -1, drop = FALSE])
TICKERS <- colnames(P)
n_obs   <- nrow(P)

## Sector labels drive the within- vs across-sector comparison in Section 8.
sectors <- data.frame(
  ticker = c("TPR", "PFE", "GILD", "MSFT", "AMZN",
             "MAR", "AAL", "BA", "BP", "AXP"),
  sector = c("Consumer/Luxury", "Pharma", "Pharma", "Technology", "Technology",
             "Hotels", "Airlines", "Aerospace", "Energy", "Financials"),
  stringsAsFactors = FALSE
)
print(sectors)

## ---- 2. Returns --------------------------------------------------------------
## Log return:    r_t = log(P_t) - log(P_{t-1})   -- additive over TIME
## Simple return: R_t = P_t / P_{t-1} - 1         -- additive across ASSETS
##
## Built by explicit row offset rather than diff() so the date alignment is
## unambiguous: row t of the return matrix corresponds to ret_dates[t].

R_log <- log(P[-1, , drop = FALSE]) - log(P[-n_obs, , drop = FALSE])
R_sim <- P[-1, , drop = FALSE] / P[-n_obs, , drop = FALSE] - 1
ret_dates <- dates[-1]

stopifnot(nrow(R_log) == n_obs - 1, !anyNA(R_log), all(is.finite(R_log)))

log_returns <- data.frame(date = ret_dates, R_log,
                          check.names = FALSE, row.names = NULL)
simple_returns <- data.frame(date = ret_dates, R_sim,
                             check.names = FALSE, row.names = NULL)

## Sanity check: log and simple returns should agree closely for small moves and
## diverge for large ones (log(1+R) < R). Largest gap tells you where it matters.
message("Max |log - simple| return gap: ",
        round(max(abs(R_log - R_sim)), 4),
        "  (expected to occur on the most extreme day)")

## ---- 3. Price dynamics -------------------------------------------------------

## (a) Normalised index, all ten on one panel. Prices differ by an order of
##     magnitude (MSFT vs AAL), so raw levels cannot be compared on one axis.
##     Rebasing to 100 at the start makes relative pandemic performance legible
##     -- this single figure carries more information per page than ten small
##     panels, which matters given the 8-page limit.
P_idx <- sweep(P, 2, P[1, ], "/") * 100

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 1), mar = c(4, 4, 3, 6))
matplot(dates, P_idx, type = "l", lty = 1, lwd = 1.4,
        col = rainbow(length(TICKERS), v = 0.85),
        xlab = "", ylab = "Index (start = 100)", xaxt = "n",
        main = "Cumulative performance, rebased to 100")
axis.Date(1, at = seq(min(dates), max(dates), by = "3 months"), format = "%b %y")
abline(h = 100, col = "grey60", lty = 3)
legend("topleft", legend = TICKERS, col = rainbow(length(TICKERS), v = 0.85),
       lty = 1, lwd = 1.4, cex = 0.7, ncol = 2, bty = "n")
par(op)

## (b) Multi-panel raw prices.
par(mfrow = c(5, 2), mar = c(2, 3.6, 1.8, 0.8), oma = c(0, 0, 2, 0))
for (tk in TICKERS) {
  plot(dates, P[, tk], type = "l", col = "steelblue4",
       xlab = "", ylab = "USD", main = tk, cex.main = 0.95)
}
mtext("Adjusted closing prices, 2020-2021", outer = TRUE, cex = 1.0)
par(op)

## (c) Multi-panel log returns. Fixing a common y-axis across panels makes the
##     differences in volatility between sectors immediately visible.
ylim_common <- range(R_log)
par(mfrow = c(5, 2), mar = c(2, 3.6, 1.8, 0.8), oma = c(0, 0, 2, 0))
for (tk in TICKERS) {
  plot(ret_dates, R_log[, tk], type = "l", col = "grey25", ylim = ylim_common,
       xlab = "", ylab = "log return", main = tk, cex.main = 0.95)
  abline(h = 0, col = "grey70")
}
mtext("Daily log returns (common scale)", outer = TRUE, cex = 1.0)
par(op)

## ---- 4. Extreme moves: top 5 losses and gains per stock ---------------------

extremes <- do.call(rbind, lapply(TICKERS, function(tk) {
  r  <- R_log[, tk]
  ord <- order(r)                      # ascending: worst first
  worst <- head(ord, 5)
  best  <- head(rev(ord), 5)
  rbind(
    data.frame(ticker = tk, type = "loss", rank = 1:5,
               date = ret_dates[worst], log_ret = r[worst]),
    data.frame(ticker = tk, type = "gain", rank = 1:5,
               date = ret_dates[best],  log_ret = r[best])
  )
}))
extremes$pct_move <- 100 * (exp(extremes$log_ret) - 1)   # log -> simple, in %
extremes$log_ret  <- round(extremes$log_ret, 4)
extremes$pct_move <- round(extremes$pct_move, 2)
rownames(extremes) <- NULL

print(extremes[extremes$type == "loss", ])
print(extremes[extremes$type == "gain", ])

## ---- 5. Largest moves overall, matched to dated events ----------------------

## Flatten the return matrix. as.vector() reads column-major, which is why the
## ticker column uses each = (number of rows).
flat <- data.frame(
  date    = rep(ret_dates, times = length(TICKERS)),
  ticker  = rep(TICKERS, each = length(ret_dates)),
  log_ret = as.vector(R_log)
)
flat$abs_move <- abs(flat$log_ret)
top_moves <- head(flat[order(-flat$abs_move), ], 10)
rownames(top_moves) <- NULL

## CANDIDATE event dates. These are the canonical macro dates for the window --
## treat them as leads to CHECK, not as established fact. Verify each one
## against a news archive (FT, Reuters) before it goes in the report, and cite
## the source. Do not assert a link you have not confirmed.
covid_events <- data.frame(
  date = as.Date(c(
    "2020-02-24", "2020-03-09", "2020-03-11", "2020-03-12", "2020-03-13",
    "2020-03-16", "2020-03-18", "2020-03-23", "2020-03-24", "2020-03-27",
    "2020-04-09", "2020-04-20", "2020-04-29", "2020-06-11", "2020-11-09",
    "2020-11-16", "2020-12-11", "2021-11-26")),
  event = c(
    "Italian outbreak; first major global sell-off",
    "OPEC+ talks collapse, Saudi-Russia oil price war; US circuit breaker",
    "WHO declares COVID-19 a pandemic",
    "US suspends travel from Europe; sharp equity sell-off",
    "US national emergency declared; sharp rebound",
    "Worst single-day S&P 500 fall since 1987",
    "Peak liquidity stress; further heavy losses",
    "Fed announces open-ended asset purchases; approximate market bottom",
    "Large relief rally on US fiscal stimulus expectations",
    "CARES Act signed into law",
    "Fed expands emergency lending facilities",
    "WTI front-month futures settle negative",
    "Gilead remdesivir NIH trial results reported",
    "Second-wave fears; renewed sell-off in travel and reopening names",
    "Pfizer/BioNTech vaccine efficacy announcement; reopening-trade rotation",
    "Moderna efficacy results (also: Pfizer Upjohn/Viatris spin-off completes)",
    "FDA emergency use authorisation for the Pfizer/BioNTech vaccine",
    "Omicron variant identified; travel and aerospace names sold off sharply"),
  stringsAsFactors = FALSE
)

## Match each large move to the nearest candidate event within a few days.
## A window is needed because news released after the close moves the NEXT
## session, and weekend announcements move Monday.
nearest_event <- function(d, ev, window = 3) {
  gap <- as.numeric(ev$date - d)
  i   <- which.min(abs(gap))
  if (abs(gap[i]) > window) return(NA_character_)
  sprintf("%s [%s, %+d day(s)]", ev$event[i], ev$date[i], -gap[i])
}

top_moves$pct_move      <- round(100 * (exp(top_moves$log_ret) - 1), 2)
top_moves$log_ret       <- round(top_moves$log_ret, 4)
top_moves$abs_move      <- NULL
top_moves$candidate_event <- vapply(top_moves$date, nearest_event,
                                    character(1), ev = covid_events)

cat("\n--- Largest absolute moves, all stocks (VERIFY EVENTS) ---\n")
print(top_moves)

cat("\n--- Three largest overall ---\n")
print(top_moves[1:3, ])

## ---- 6. Moments and normality statistics ------------------------------------
## moments::kurtosis() returns RAW kurtosis (3 for a normal), so excess kurtosis
## is kurtosis - 3. Getting this wrong is the single most common error here.

stats_table <- data.frame(
  ticker       = TICKERS,
  mean_daily   = colMeans(R_log),
  sd_daily     = apply(R_log, 2, sd),
  ann_return   = colMeans(R_log) * 252,
  ann_vol      = apply(R_log, 2, sd) * sqrt(252),
  skewness     = apply(R_log, 2, moments::skewness),
  exc_kurtosis = apply(R_log, 2, moments::kurtosis) - 3,
  min          = apply(R_log, 2, min),
  max          = apply(R_log, 2, max),
  row.names    = NULL
)

## ---- 7. Formal normality tests ----------------------------------------------

## Jarque-Bera, following the Workshop 7 wrapper. Null: the margin is normal.
## The statistic is (n/6)(S^2 + (K-3)^2/4), asymptotically chi-squared with 2 df.
JBTest <- function(dataMatrix) {
  out <- data.frame(ticker = colnames(dataMatrix),
                    JB_stat = NA_real_, JB_pvalue = NA_real_)
  for (i in seq_len(ncol(dataMatrix))) {
    tst <- tseries::jarque.bera.test(dataMatrix[, i])
    out$JB_stat[i]   <- as.numeric(tst$statistic)
    out$JB_pvalue[i] <- tst$p.value
  }
  out
}
jb <- JBTest(R_log)

## Shapiro-Wilk as an independent cross-check (valid for n <= 5000).
sw_p <- apply(R_log, 2, function(x) shapiro.test(x)$p.value)

stats_table$JB_stat   <- jb$JB_stat
stats_table$JB_pvalue <- jb$JB_pvalue
stats_table$SW_pvalue <- as.numeric(sw_p)

## Display copy, rounded. Very small p-values are reported as "< 2e-16" rather
## than as a spuriously precise number.
display <- stats_table
num <- sapply(display, is.numeric)
display[num] <- lapply(display[num], function(x) signif(x, 4))
display$JB_pvalue <- ifelse(stats_table$JB_pvalue < 2e-16, "< 2e-16",
                            format(signif(stats_table$JB_pvalue, 3)))
display$SW_pvalue <- ifelse(stats_table$SW_pvalue < 2e-16, "< 2e-16",
                            format(signif(stats_table$SW_pvalue, 3)))

cat("\n--- Descriptive statistics and normality tests ---\n")
print(display, row.names = FALSE)

cat("\nMargins rejecting normality at 5% (Jarque-Bera): ",
    sum(stats_table$JB_pvalue < 0.05), " of ", length(TICKERS), "\n", sep = "")

## Joint normality. Task 3's normal portfolio model assumes MULTIVARIATE
## normality, which testing ten margins separately does not establish.
if (has_QRM) {
  cat("\n--- Mardia test of joint normality ---\n")
  print(MardiaTest(R_log))
}

## ---- 8. QQ-plots -------------------------------------------------------------
## Returns are standardised so all ten panels share a scale and the reference
## line is exactly y = x. Departures at the ends are fat tails; asymmetry
## between the two ends is skew.

par(mfrow = c(5, 2), mar = c(3.4, 3.4, 1.8, 0.8), mgp = c(2, 0.7, 0),
    oma = c(0, 0, 2, 0))
for (tk in TICKERS) {
  z <- as.numeric(scale(R_log[, tk]))
  qqnorm(z, main = tk, cex.main = 0.95, pch = 19, cex = 0.4,
         col = rgb(0, 0, 0, 0.45),
         xlab = "Theoretical quantiles", ylab = "Standardised return")
  abline(0, 1, col = "red", lwd = 1.3)
}
mtext("Normal QQ-plots of standardised daily log returns", outer = TRUE, cex = 1.0)
par(op)

## ---- 9. Pairwise dependence -------------------------------------------------

## Full correlation matrix, ordered by sector so blocks are visible.
ord_by_sector <- sectors$ticker[order(sectors$sector)]
C <- cor(R_log[, ord_by_sector])
cat("\n--- Correlation matrix of daily log returns (sector-ordered) ---\n")
print(round(C, 2))

## Helper for one scatter panel: correlation in the title, OLS fit, zero lines.
scatter_pair <- function(a, b, tag) {
  x <- R_log[, a]; y <- R_log[, b]
  plot(x, y, pch = 19, cex = 0.45, col = rgb(0, 0, 0, 0.35),
       xlab = paste(a, "log return"), ylab = paste(b, "log return"),
       main = sprintf("%s vs %s (%s), rho = %.2f", a, b, tag, cor(x, y)),
       cex.main = 0.9)
  abline(h = 0, v = 0, col = "grey85")
  abline(lm(y ~ x), col = "red", lwd = 1.3)
}

within_pairs <- list(c("PFE", "GILD"), c("MSFT", "AMZN"), c("AAL", "BA"))
across_pairs <- list(c("AAL", "MSFT"), c("MAR", "PFE"), c("BP", "MSFT"))

par(mfrow = c(2, 3), mar = c(4, 4, 2.6, 1))
for (pr in within_pairs) scatter_pair(pr[1], pr[2], "within")
for (pr in across_pairs) scatter_pair(pr[1], pr[2], "across")
par(op)

## Correlation breakdown: average pairwise correlation in the crisis window
## versus in the calmer later period. Diversification is not a fixed property.
crisis <- ret_dates >= as.Date("2020-02-19") & ret_dates <= as.Date("2020-04-30")
calm   <- ret_dates >= as.Date("2021-01-01")

avg_cor <- function(idx) {
  M <- cor(R_log[idx, ])
  mean(M[lower.tri(M)])
}
cat("\nAverage pairwise correlation\n")
cat("  Crisis (19 Feb - 30 Apr 2020, n =", sum(crisis), "): ",
    round(avg_cor(crisis), 3), "\n")
cat("  Calmer (2021, n =", sum(calm), "):              ",
    round(avg_cor(calm), 3), "\n")

## Empirical lower-tail dependence: given stock A is in its worst 5% of days,
## how often is B also? Under independence this would be about 5%. This is the
## quantity the copula in Task 4 has to reproduce, so computing it now gives you
## a benchmark to judge the fitted copulas against.
tail_codep <- function(a, b, q = 0.05) {
  x <- R_log[, a]; y <- R_log[, b]
  ix <- x <= quantile(x, q); iy <- y <= quantile(y, q)
  data.frame(pair = paste(a, b, sep = "-"),
             joint_days = sum(ix & iy),
             expected_if_indep = round(q * q * length(x), 1),
             lambda_L_hat = round(sum(ix & iy) / sum(ix), 3))
}
cat("\n--- Empirical lower-tail co-exceedance (5% level) ---\n")
print(do.call(rbind, lapply(c(within_pairs, across_pairs),
                            function(pr) tail_codep(pr[1], pr[2]))),
      row.names = FALSE)

## ---- 10. Save ----------------------------------------------------------------

saveRDS(list(log_returns = log_returns, simple_returns = simple_returns,
             R_log = R_log, ret_dates = ret_dates,
             stats_table = stats_table, extremes = extremes,
             top_moves = top_moves, sectors = sectors,
             correlation = cor(R_log)),
        "data/task2_results.rds")

write.csv(log_returns,    "data/log_returns.csv",    row.names = FALSE)
write.csv(simple_returns, "data/simple_returns.csv", row.names = FALSE)
write.csv(stats_table,    "data/task2_stats.csv",    row.names = FALSE)

message("Task 2 objects saved to data/task2_results.rds")